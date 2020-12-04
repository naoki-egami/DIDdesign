

#' Function to implement the Staggered Adption Design
#' @param formula A formula indicating outcome and treatment.
#' @param data A long form panel data.
#' @param id_subject A character variable indicating subject index.
#' @param id_time A character variable indicating time index.
#' @examples
#' # simulate data
#' set.seed(1234)
#' dat <- simulate_sad(30, 10)
#'
#' # estimate
#' sa_did_est(outcome ~ treatment, data = dat$dat,
#'   id_subject = "id_subject", id_time = 'id_time'
#' )
#' @importFrom dplyr %>% as_tibble group_by
#' @importFrom future.apply future_lapply
#' @importFrom future plan multiprocess sequential
#' @export
sa_did_est <- function(formula, data, id_subject, id_time, n_boot = 10, lead = 0, parallel = TRUE) {

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
  est <- sa_double_did(dat_panel, treatment, outcome, lead)


  ## --------------------------------------
  ## obtain the weighting matrix via bootstrap
  ## --------------------------------------
  if (isTRUE(parallel)) {
    plan(multiprocess)
  } else {
    plan(sequential)
  }

  # est_boot <- matrix(NA, nrow = n_boot, ncol = 2)
  # for (b in 1:n_boot) {
  #   dat_boot <- sample_panel(dat_panel)
  #   est_boot[b,] <- sa_double_did(dat_boot, treatment, outcome, lead)
  # }

  est_boot <- do.call(rbind, future_lapply(1:n_boot, function(i) {
      dat_boot <- sample_panel(dat_panel)
      sa_double_did(dat_boot, treatment, outcome, lead)
  }, future.seed = TRUE))


  ## estimate GMM weighting matrix
  W <- cov(est_boot)
  w_did  <- (W[1,1] - W[1,2]) / (sum(diag(W)) - 2 * W[1,2])
  w_sdid <- (W[2,2] - W[1,2]) / (sum(diag(W)) - 2 * W[1,2])
  w_vec  <- c(w_did, w_sdid)

  ## double did estimate
  double_did <- as.vector(t(w_vec) %*% est)

  ## compute the variance of double did
  var_ddid <- as.vector(t(w_vec^2) %*% apply(est_boot, 2, var))

  return(list(est = est, res_boot = est_boot,
    ddid = double_did, ddid_var = var_ddid, weights = w_vec))
}


sa_double_did <- function(dat_panel, treatment, outcome, lead) {
  ## create Gmat and index for each design
  Gmat <- create_Gmat(dat_panel, treatment = treatment)
  id_time_use <- get_periods(Gmat)
  id_subj_use <- get_subjects(Gmat, id_time_use)
  time_weight <- get_time_weight(Gmat, id_time_use)

  ## estimate DID and sDID for each periods
  out <- compute_did(dat_panel, outcome, treatment,
                     id_time_use, id_subj_use, time_weight, lead = lead)
  return(out)
}



#' Create a G matrix
#'
#' @param dat_panel A class of \code{panelr} object.
#' @importFrom rlang sym !!
#' @importFrom dplyr %>% select mutate case_when
#' @importFrom tidyr pivot_wider
#' @keywords internal
create_Gmat <- function(dat_panel, treatment) {
  Gmat <- dat_panel %>%
    group_by(id_subject) %>%
    mutate(g_sum = max(id_time) - sum(!!sym(treatment)) + 1) %>%
    mutate(g = case_when(
      g_sum > id_time ~ 0,
      g_sum == id_time ~ 1,
      g_sum < id_time ~ -1
    )) %>%
    select(g, id_subject, id_time) %>%
    pivot_wider(id_cols = id_subject, names_from = id_time, values_from = g) %>%
    ungroup() %>%
    select(-id_subject) %>%
    data.matrix()

  return(Gmat)
}



#' Obtain periods used for the analysis
#' @param Gmat G matrix produced in \code{create_Gmat()}.
#' @param thres A minimum number of treatd units for the period included in the analysis. Default is 2.
#' @keywords internal
get_periods <- function(Gmat, thres = 1) {
  ## check which periods to use
  ## only use periods that are more than "thres" observations treated
  use_id <- which(apply(Gmat, 2, function(x) sum(x == 1) >= thres))

  ## remove the last period if everyone is eventually treated
  n_treated <- apply(Gmat, 2, function(x) sum(x == 1))
  if (sum(n_treated[n_treated >= thres]) == nrow(Gmat)) {
    use_id <- use_id[-length(use_id)]
  }
  return(use_id)
}


#' Obtain subject index for each periods
#' @param Gamt G matrix created by \code{create_Gmat()}.
#' @param id_time_use A vector of time index. Should be normalized.
#' @keywords internal
get_subjects <- function(Gmat, id_time_use) {
  id_use <- list(); iter <- 1
  for (i in 1:length(id_time_use)) {
    tmp <- Gmat[,id_time_use[i]]
    id_use[[iter]] <- which(tmp >= 0)
    iter <- iter + 1
  }

  return(id_use)
}


#' Compute time-sepcific weights
#' @param Gmat G matrix.
#' @param id_time_use A vector of time index (normalized). Output of \code{get_periods}.
#' @keywords internal
get_time_weight <- function(Gmat, id_time_use) {
  n_treated   <- sum(Gmat[,id_time_use] == 1)
  time_weight <- rep(NA, length(id_time_use))
  for (i in 1:length(id_time_use)) {
    tmp <- Gmat[,id_time_use[i]]
    time_weight[i] <- sum(tmp == 1) / n_treated
  }
  return(time_weight)
}


#' Estimate DID
#' @keywords internal
#' @import plm
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


#' Block Bootstrap
#' @importFrom tibble as_tibble
#' @keywords internal
sample_panel <- function(panel_dat) {


  ## get subject id for bootstrap
  id_vec     <- unique(pull(panel_dat, id_subject))
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
    select(-id_subject) %>%
    mutate(id_subject = id_new)

  return(boot_dat)
}
