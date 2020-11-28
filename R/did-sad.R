

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
#' @export
sa_did_est <- function(formula, data, id_subject, id_time, n_boot = 100, lead = 0) {

  ## keep track of long-form data with panel class from \code{panelr} package
  dat_panel <- panel_data(data, id = id_subject, wave = id_time)

  ## extract variable informations
  if (!("formula" %in% class(formula))) {
    formula <- as.formula(formula)
  }
  all_vars  <- all.vars(formula)
  outcome   <- all_vars[1]
  treatment <- all_vars[2]

  ## estimate DID and sDID
  est <- sa_double_did(dat_panel, treatment, outcome, lead)

  ## bootstrap
  est_boot <- matrix(NA, nrow = n_boot, ncol = 2)
  for (b in 1:n_boot) {
    dat_boot <- sample_panel(dat_panel)
    est_boot[b,] <- sa_double_did(dat_boot, treatment, outcome, lead)
  }

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
  did <- compute_did(dat_panel, outcome, treatment,
                     id_time_use, id_subj_use, time_weight, lead = lead )
  sdid <- compute_sdid(dat_panel, outcome, treatment,
                     id_time_use, id_subj_use, time_weight, lead = lead)
  out <- c(did, sdid); names(out) <- c('DID', "sDID")

  attr(out, "did_all") <- attr(did, "all_est")
  return(out)
}



#' Create a G matrix
#'
#' @param dat_panel A class of \code{panelr} object.
#' @importFrom rlang sym !!
#' @importFrom dplyr %>% select mutate
#' @keywords internal
create_Gmat <- function(dat_panel, treatment) {
  Gmat <- dat_panel %>%
    mutate(g_sum = max(!!sym(get_wave(dat_panel))) - sum(!!sym(treatment)) + 1) %>%
    mutate(g = ifelse(g_sum < !!sym(get_wave(dat_panel)), 0,
                    ifelse(g_sum == !!sym(get_wave(dat_panel)), 1, -1))) %>%
    select(g) %>%
    widen_panel() %>%
    select(-id_subject) %>%
    data.matrix()

  return(Gmat)
}



#' Obtain periods used for the analysis
#' @param Gmat G matrix produced in \code{create_Gmat()}.
#' @param thres A minimum number of treatd units for the period included in the analysis. Default is 2.
#' @keywords internal
get_periods <- function(Gmat, thres = 3) {
  ## check which periods to use
  ## only use periods that are more than "thres" observations treated
  use_id <- which(apply(Gmat, 2, function(x) sum(x == 1) >= thres))
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
    id_use[[iter]] <- which(tmp != 0)
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
compute_did <- function(dat_panel, outcome, treatment,
  id_time_use, id_subj_use, time_weight, min_time = 3, lead) {
  est <- list(); iter <- 1

  ## we need to renormalize weights if past periods is not avaialbe
  time_weight_new <- list()
  max_time <- max(dat_panel %>% pull(get_wave(dat_panel)))

  ## compute individual time specific DID
  for (i in 1:length(id_time_use)) {
    if ((id_time_use[i] >= min_time) && ((id_time_use[i] + lead) <= max_time)) {
      ## subset the data
      dat_use <- dat_panel %>%
        filter(!!sym(get_id(dat_panel)) %in% id_subj_use[[i]]) %>%
        filter(!!sym(get_wave(dat_panel)) == (id_time_use[i] + lead) |
               !!sym(get_wave(dat_panel)) == id_time_use[i]-1) %>%
        as_pdata.frame()

      ## estimate DID using plm
      fm  <- as.formula(paste(outcome, treatment, sep = "~"))
      fit <- plm(fm, data = dat_use,
                  effects = 'twoways', model = 'within')
      est[[iter]] <- fit$coef
      time_weight_new[[iter]] <- time_weight[iter]
      iter <- iter + 1
    }
  }


  ## combine all did estimates
  time_weight_new <- unlist(time_weight_new)
  did_vec <- as.vector(unlist(est) %*% (time_weight_new / sum(time_weight_new)))

  attr(did_vec, "all_est") <- unlist(est)
  return(did_vec)
}


#' Estimate sDID
#' @keywords internal
compute_sdid <- function(dat_panel, outcome, treatment,
  id_time_use, id_subj_use, time_weight, min_time = 3, lead) {
  est <- list(); iter <- 1

  ## we need to renormalize weights if past periods is not avaialbe
  time_weight_new <- list()

  max_time <- max(dat_panel %>% pull(get_wave(dat_panel)))

  ## compute individual time specific DID
  for (i in 1:length(id_time_use)) {
    if ((id_time_use[i] >= min_time) && ((id_time_use[i] + lead) <= max_time)) {
      ## subset the data
      dat_use <- dat_panel %>%
        mutate(outcome_diff = c(NA, diff(!!sym(outcome)))) %>%
        filter(!!sym(get_id(dat_panel)) %in% id_subj_use[[i]]) %>%
        filter(!!sym(get_wave(dat_panel)) == (id_time_use[i] + lead) |
               !!sym(get_wave(dat_panel)) == id_time_use[i]-1) %>%
        as_pdata.frame()

      ## estimate DID using plm
      fm  <- as.formula(paste("outcome_diff", treatment, sep = "~"))
      fit <- plm(fm, data = dat_use,
                  effects = 'twoways', model = 'within')
      est[[iter]] <- fit$coef
      time_weight_new[[iter]] <- time_weight[iter]
      iter <- iter + 1
    }
  }


  ## combine all did estimates
  time_weight_new <- unlist(time_weight_new)
  sdid_vec <- as.vector(unlist(est) %*% (time_weight_new / sum(time_weight_new)))

  return(sdid_vec)
}



#' Block Bootstrap
#' @importFrom tibble as_tibble
#' @keywords internal
sample_panel <- function(panel_dat) {


  ## get subject id for bootstrap
  id_subject <- get_id(panel_dat)
  id_vec     <- unique(pull(panel_dat, !!sym(id_subject)))
  id_boot    <- sample(id_vec, size = length(id_vec), replace = TRUE)

  ## constrcut a panel data for bootstrap
  tmp <- as.data.frame(panel_dat)
  dat_list <- list()
  for (i in 1:length(id_boot)) {
    dat_list[[i]] <- filter(tmp, !!sym(id_subject) == id_boot[i]) %>%
      mutate(id_new = i, id_old = id_boot[i])
  }

  boot_dat <- do.call("rbind", dat_list) %>%
    as_tibble() %>%
    select(-id_subject) %>%
    mutate(id_subject = id_new) %>%
    panel_data(id = id_subject, wave = !!sym(get_wave(panel_dat)))

  return(boot_dat)
}
