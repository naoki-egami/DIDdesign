

#' Function to implement DDID under the Staggered Adption Design
#' @param formula A formula indicating outcome and treatment.
#' @param data A long form panel data.
#' @param id_subject A character variable indicating subject index.
#' @param id_time A character variable indicating time index.
#' @importFrom dplyr %>% as_tibble group_by
#' @importFrom tibble tibble
#' @importFrom future.apply future_lapply
#' @importFrom future plan multiprocess sequential
#' @keywords internal
did_sad <- function(formula, data, id_subject, id_time, option) {

  ## keep track of long-form data with panel class from \code{panelr} package
  # dat_panel <- panel_data(data, id = id_subject, wave = id_time)
  ## assign group index
  ## normalize the time index
  dat_panel <- as_tibble(data) %>%
    rename(id_subject = !!sym(id_subject), id_time = !!sym(id_time)) %>%
    mutate(id_time = as.numeric(as.factor(id_time)))

  ## extract variable informations
  if (!("formula" %in% class(formula))) {
    formula <- as.formula(formula)
  }
  all_vars  <- all.vars(formula)
  outcome   <- all_vars[1]
  treatment <- all_vars[2]

  ## --------------------------------------
  ## obtain point estimates
  ## estimate DID and sDID
  ## --------------------------------------
  est <- sa_double_did(dat_panel, treatment, outcome, option$lead)


  ## --------------------------------------
  ## obtain the weighting matrix via bootstrap
  ## --------------------------------------
  if (isTRUE(option$parallel)) {
    plan(multicore)
  } else {
    plan(sequential)
  }

  est_boot <- do.call(rbind, future_lapply(1:option$n_boot, function(i) {
      dat_boot <- sample_panel(dat_panel)
      sa_double_did(dat_boot, treatment, outcome, option$lead)
  }, future.seed = TRUE))

  ## --------------------------------------
  ## double DID estimator
  ## --------------------------------------

  ## estimate GMM weighting matrix
  W      <- cov(est_boot)
  w_did  <- (W[1,1] - W[1,2]) / (sum(diag(W)) - 2 * W[1,2])
  w_sdid <- (W[2,2] - W[1,2]) / (sum(diag(W)) - 2 * W[1,2])
  w_vec  <- c(w_did, w_sdid)

  ## double did estimate
  double_did <- as.vector(t(w_vec) %*% est)

  ## compute the variance of double did
  var_ddid <- as.vector(t(w_vec^2) %*% apply(est_boot, 2, var))

  ## --------------------------------------
  ## prepare output
  ## --------------------------------------
  estimates <- tibble(
    estimator = c("SA-Double-DID", "SA-DID", "SA-sDID"),
    estimate  = c(est, double_did),
    std.error = c(apply(est_boot, 2, sd), sqrt(var_ddid))
  )

  weights <- list(W = W, weight_did = w_did, weight_sdid = w_sdid)
  return(list(estimates = estimates, weights = weights))
}


#' Compute the time-weighted DID and sDID
#' @keywords internal
#' @return A vector of two elements: Time-weighted DID (the first element) and time-weighted sequential DID (second).
sa_double_did <- function(dat_panel, treatment, outcome, lead) {
  ## create Gmat and index for each design
  Gmat <- create_Gmat(dat_panel, treatment = treatment)
  id_time_use <- get_periods(Gmat)
  id_subj_use <- get_subjects(Gmat, id_time_use)
  time_weight <- get_time_weight(Gmat, id_time_use)

  ## estimate DID and sDID for each periods and weight them by the proportion of treated units
  out <- compute_did(dat_panel, outcome, treatment,
                     id_time_use, id_subj_use, time_weight, lead = lead)
  return(out)
}


#' Estimate DID and sDID
#' @keywords internal
#' @importFrom plm plm
#' @importFrom dplyr lag
compute_did <- function(dat_panel, outcome, treatment,
  id_time_use, id_subj_use, time_weight, min_time = 3, lead) {
  est_did <- est_sdid <- list(); iter <- 1

  ## we need to renormalize weights if past periods is not avaialbe
  time_weight_new <- list()
  max_time <- max(dat_panel$id_time)

  ## compute individual time specific DID
  for (i in 1:length(id_time_use)) {
    if ((id_time_use[i] >= min_time) && ((id_time_use[i] + lead) <= max_time)) {

      ## -------------------------------
      ## compute DID
      ## -------------------------------
      ## subset the data
      dat_use <- dat_panel %>%
        filter(id_subject %in% id_subj_use[[i]]) %>%
        filter(id_time == (id_time_use[i] + lead) |
               id_time == id_time_use[i]-1)

      ## estimate DID using plm
      fm  <- as.formula(paste(outcome, treatment, sep = "~"))
      fit <- plm::plm(fm, data = dat_use, index = c("id_subject", "id_time"),
                      effects = 'twoways', model = 'within')
      est_did[[iter]] <- fit$coef

      ## -------------------------------
      ## compute sDID
      ## -------------------------------
      dat_use <- dat_panel %>%
        group_by(id_subject) %>%
        mutate(outcome_lag = lag(!!sym(outcome), order_by = id_time)) %>%
        mutate(outcome_diff = !!sym(outcome) - outcome_lag) %>%
        filter(id_subject %in% id_subj_use[[i]]) %>%
        filter(id_time == (id_time_use[i] + lead) |
               id_time == id_time_use[i]-1)

      ## estimate DID using plm
      fm  <- as.formula(paste("outcome_diff", treatment, sep = "~"))
      fit <- plm(fm, data = dat_use, index = c("id_subject", "id_time"),
                  effects = 'twoways', model = 'within')

      est_sdid[[iter]] <- fit$coef

      ## -------------------------------
      ## time weights
      ## -------------------------------
      time_weight_new[[iter]] <- time_weight[iter]
      iter <- iter + 1
    }
  }

  ## -------------------------------
  ## combine all did estimates
  ## -------------------------------
  time_weight_new <- unlist(time_weight_new)
  did_vec  <- as.vector(unlist(est_did) %*% (time_weight_new / sum(time_weight_new)))
  sdid_vec <- as.vector(unlist(est_sdid) %*% (time_weight_new / sum(time_weight_new)))

  res <- c(did_vec, sdid_vec); names(res) <- c("DID", "sDID")
  return(res)
}


#' Block Bootstrap used in SA design
#' @importFrom tibble as_tibble
#' @importFrom dplyr %>% filter mutate select pull
#' @keywords internal
sample_panel <- function(panel_dat) {

  ## get subject id for bootstrap
  id_vec     <- unique(pull(panel_dat, .data$id_subject))
  id_boot    <- sample(id_vec, size = length(id_vec), replace = TRUE)

  ## constrcut a panel data for bootstrap
  tmp <- as.data.frame(panel_dat)
  dat_list <- list()
  for (i in 1:length(id_boot)) {
    dat_list[[i]] <- filter(tmp, id_subject == id_boot[i]) %>%
      mutate(id_new = i, id_old = id_boot[i])
  }

  boot_dat <- do.call("rbind", dat_list) %>%
    as_tibble() %>%
    select(-.data$id_subject) %>%
    mutate(id_subject = .data$id_new)

  return(boot_dat)
}
