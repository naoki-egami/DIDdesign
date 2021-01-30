

#' Function to implement DDID under the Staggered Adoption Design
#' @param formula A formula indicating outcome and treatment.
#' @param data A long form panel data.
#' @param id_subject A character variable indicating subject index.
#' @param id_time A character variable indicating time index.
#' @return Double DID estimates.
#' @importFrom dplyr %>% as_tibble group_by arrange
#' @importFrom tibble tibble
#' @importFrom future.apply future_lapply
#' @importFrom future plan multicore sequential
#' @keywords internal
did_sad <- function(formula, data, id_subject, id_time, option) {

  ## keep track of long-form data with panel class from \code{panelr} package
  # dat_panel <- panel_data(data, id = id_subject, wave = id_time)
  ## assign group index
  ## normalize the time index
  dat_panel <- as_tibble(data) %>%
    rename(id_subject = !!sym(id_subject), id_time = !!sym(id_time)) %>%
    mutate(id_time = as.numeric(as.factor(id_time))) %>%
    arrange(id_subject, id_time)

  is_panel <- TRUE
  var_cluster <- option$id_cluster
  if (is.null(var_cluster) && isTRUE(is_panel)) {
    var_cluster <- "id_unit"
    option$var_cluster_pre <- id_subject
  }

  ## extract variable informations
  if (!("formula" %in% class(formula))) {
    formula <- as.formula(formula)
  }

  all_vars  <- all.vars(formula)
  outcome   <- all_vars[1]
  treatment <- all_vars[2]

  ## prepare custom formula
  fm_prep <- did_formula(formula, is_panel)

  ## --------------------------------------
  ## obtain point estimates
  ## estimate DID and sDID
  ## --------------------------------------
  est <- sa_double_did(fm_prep, dat_panel, treatment, outcome, option)

  ## --------------------------------------
  ## obtain the weighting matrix via bootstrap
  ## --------------------------------------
  if (isTRUE(option$parallel)) {
    plan(multicore)
  } else {
    plan(sequential)
  }

  est_boot <- future_lapply(1:option$n_boot, function(i) {
      dat_boot <- sample_panel(dat_panel)
      sa_double_did(fm_prep, dat_boot, treatment, outcome, option)
  }, future.seed = TRUE)

  ## --------------------------------------
  ## double DID estimator
  ## --------------------------------------
  estimates <- sa_did_to_ddid(est, est_boot, lead = option$lead)
  return(estimates)
}


#' Compute the time-weighted DID and sDID
#' @keywords internal
#' @param formula A formula.
#' @param dat_panel A panel data.
#' @param treatment Name of the treatment variable.
#' @param outcome Name of the outcome variable.
#' @param option A list of options.
#' @return A vector of two elements: Time-weighted DID (the first element) and time-weighted sequential DID (second).
sa_double_did <- function(formula, dat_panel, treatment, outcome, option) {
  ## create Gmat and index for each design
  Gmat <- create_Gmat(dat_panel, treatment = treatment)
  id_time_use <- get_periods(Gmat, option$thres)
  id_subj_use <- get_subjects(Gmat, id_time_use)
  time_weight <- get_time_weight(Gmat, id_time_use)

  ## estimate DID and sDID for each periods and weight them by the proportion of treated units
  out <- compute_did(formula, dat_panel, outcome, treatment,
                     id_time_use, id_subj_use, time_weight,
                     lead = option$lead, var_cluster_pre = option$var_cluster_pre)
  return(out)
}


#' Compute lead specific variance covariance matrix
#' @param obj A list of bootstrap outputs.
#' @param lead A vector of the lead parameter.
#' @return Estimate of the variance covariance matrix.
#' @keywords internal
#' @importFrom purrr map
sa_calc_cov <- function(obj, lead) {
  tmp <- map(obj, ~cbind(.x$DID[[lead]], .x$sDID[[lead]]))
  cov_est <- cov(do.call(rbind, tmp))
  return(cov_est)
}


#' Convert estimates into double did
#' @keywords internal
#' @param obj_point Point estimates.
#' @param obj_boot Bootstrap outputs.
#' @param lead A vector of the lead parameter.
#' @return A list of estimates and weights.
#' @importFrom dplyr as_tibble
sa_did_to_ddid <- function(obj_point, obj_boot, lead) {

  estimates <- W_save <- vector("list", length = length(lead))


  for (ll in 1:length(lead)) {
    tmp <- matrix(NA, nrow = 3, ncol = 3)

    ## compute variance covariance matrix
    W <- sa_calc_cov(obj_boot, lead[ll]+1)

    ## compute weight from Vcov
    w_did  <- (W[1,1] - W[1,2]) / (sum(diag(W)) - 2 * W[1,2])
    w_sdid <- (W[2,2] - W[1,2]) / (sum(diag(W)) - 2 * W[1,2])
    w_vec  <- c(w_did, w_sdid)

    ## compute double did
    est  <- c(obj_point$DID[[ll]], obj_point$sDID[[ll]])
    ddid <- t(w_vec) %*% est

    ## variance
    var_ddid <- as.vector(t(w_vec^2) %*% diag(W))

    ## save weights
    estimates[[ll]] <- data.frame(
      estimator = c("SA-Double-DID", "SA-DID", "SA-sDID"),
      lead      = lead[ll],
      estimate  = c(ddid, est),
      std.error = sqrt(c(var_ddid, diag(W))),
      ddid_weight = c(NA, w_vec)
    )

    W_save[[ll]] <- solve(W)
  }

  ## summarise
  estimates <- as_tibble(bind_rows(estimates))
  return(list(estimates = estimates, weights = W_save))

}

#' Estimate DID and sDID
#' @keywords internal
#' @inheritParams sa_double_did
#' @param id_time_use A list of time index.
#' @param id_subj_use A list of unit index.
#' @param time_weight A vector of time weights.
#' @param min_time A threshold value of minimum number of units in the treatment.
#' @param lead A vector of lead parameters.
#' @param var_cluster_pre A variable name used for clustering.
#' @return A list of DID and sDID estimates.
#' @importFrom dplyr lag bind_rows filter select left_join
#' @importFrom rlang !! sym
#' @importFrom purrr map_dbl
compute_did <- function(formula, dat_panel, outcome, treatment,
  id_time_use, id_subj_use, time_weight, min_time = 3, lead, var_cluster_pre) {

  est_did <- list()
  iter <- 1

  ## we need to renormalize weights if past periods is not avaialbe
  time_weight_new <- list()
  max_time <- max(dat_panel$id_time)

  ## compute individual time specific DID
  for (i in 1:length(id_time_use)) {
    if ((id_time_use[i] >= min_time) && ((id_time_use[i] + lead) <= max_time)) {

      est_did[[iter]] <- vector("list", length = length(lead))

      ## -------------------------------
      ## compute DID
      ## -------------------------------
      ## subset the data
      dat_use <- dat_panel %>%
        filter(.data$id_subject %in% id_subj_use[[i]]) %>%
        filter(.data$id_time == (id_time_use[[i]]) |
               .data$id_time == id_time_use[[i]]-1 |
               .data$id_time == id_time_use[[i]]-2)

      if (max(lead) > 0) {
        ## handle lead cases
        treatment_info <- dat_panel %>%
          filter(.data$id_subject %in% id_subj_use[[i]]) %>%
          filter(.data$id_time == id_time_use[[i]]) %>%
          select(.data$id_subject, !!sym(treatment))

        outcome_lead <- dat_panel %>%
          filter(.data$id_subject %in% id_subj_use[[i]]) %>%
          filter(.data$id_time <= (id_time_use[[i]] + max(lead)) &
                 .data$id_time >= (id_time_use[[i]] + max(1, min(lead)))) %>%
          select(-!!sym(treatment)) %>%
          left_join(treatment_info, by = 'id_subject')

        dat_use <- bind_rows(dat_use, outcome_lead)
      }


      ## -------------------------------
      ## estimate DID and sDID
      ## -------------------------------
      dat_did <- did_panel_data(
        formula$var_outcome, formula$var_treat, formula$var_covars,
        var_cluster_pre, id_unit = "id_subject", id_time = "id_time", dat_use
      )

      for (ll in 1:length(lead)) {
        est_did[[iter]][[ll]]  <- ddid_fit(formula$fm_did, dat_did, lead = lead[ll])
      }

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

  did_vec <- sdid_vec <- list()
  for (ll in 1:length(lead)) {
    did_est <- map_dbl(est_did, ~.x[[ll]][1])
    sdid_est <- map_dbl(est_did, ~.x[[ll]][2])
    did_vec[[ll]] <- did_est %*% (time_weight_new / sum(time_weight_new))
    sdid_vec[[ll]] <- sdid_est %*% (time_weight_new / sum(time_weight_new))
  }


  res <- list("DID" = did_vec, "sDID" = sdid_vec)
  return(res)
}


#' Block Bootstrap used in SA design
#' @param panel_dat A panel data.
#' @return A panel data, sampled from the input data.
#' @importFrom tibble as_tibble
#' @importFrom dplyr %>% filter mutate select pull
#' @keywords internal
sample_panel <- function(panel_dat) {

  ## get subject id for bootstrap
  id_vec  <- unique(pull(panel_dat, .data$id_subject))
  id_boot <- sample(id_vec, size = length(id_vec), replace = TRUE)

  ## constrcut a panel data for bootstrap
  tmp <- as.data.frame(panel_dat)
  dat_list <- list()
  for (i in 1:length(id_boot)) {
    dat_list[[i]] <- filter(tmp, .data$id_subject == id_boot[i]) %>%
      mutate(id_new = i, id_old = id_boot[i])
  }

  boot_dat <- do.call("rbind", dat_list) %>%
    as_tibble() %>%
    select(-.data$id_subject) %>%
    mutate(id_subject = .data$id_new)

  return(boot_dat)
}
