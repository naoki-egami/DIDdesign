

#' Function to implement the Staggered Adption Design
#' @param formula A formula indicating outcome and treatment.
#' @param data A long form panel data.
#' @param id_subject A character variable indicating subject index.
#' @param id_time A character variable indicating time index.
#' @import panelr
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
sa_did_est <- function(formula, data, id_subject, id_time) {

  ## keep track of long-form data with panel class from \code{panelr} package
  dat_panel <- panel_data(data, id = id_subject, wave = id_time)

  ## extract variable informations
  all_vars  <- all.vars(formula)
  outcome   <- all_vars[1]
  treatment <- all_vars[2]

  ## create Gmat and index for each design
  Gmat <- create_Gmat(dat_panel, treatment = treatment)
  id_time_use <- get_periods(Gmat)
  id_subj_use <- get_subjects(Gmat, id_time_use)
  time_weight <- get_time_weight(Gmat, id_time_use)

  ## estimate DID and sDID for each periods
  did <- compute_did(dat_panel, outcome, treatment,
                     id_time_use, id_subj_use, time_weight)
  sdid <- compute_sdid(dat_panel, outcome, treatment,
                     id_time_use, id_subj_use, time_weight)
  out <- c(did, sdid); names(out) <- c('DID', "sDID")
  return(out)
}





#' Create a G matrix
#'
#' @param dat_panel A class of \code{panelr} object.
#' @import panelr
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
get_periods <- function(Gmat, thres = 2) {
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
  id_time_use, id_subj_use, time_weight) {
  est <- list(); iter <- 1

  ## we need to renormalize weights if past periods is not avaialbe
  time_weight_new <- list()

  ## compute individual time specific DID
  for (i in 1:length(id_time_use)) {
    if (id_time_use[i] >= 2) {
      ## subset the data
      dat_use <- dat_panel %>%
        filter(!!sym(get_id(dat_panel)) %in% id_subj_use[[i]]) %>%
        filter(!!sym(get_wave(dat_panel)) <= id_time_use[i]) %>%
        filter(!!sym(get_wave(dat_panel)) >= id_time_use[i]-1) %>%
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

  return(did_vec)
}


#' Estimate sDID
#' @keywords internal
compute_sdid <- function(dat_panel, outcome, treatment,
  id_time_use, id_subj_use, time_weight) {
  est <- list(); iter <- 1

  ## we need to renormalize weights if past periods is not avaialbe
  time_weight_new <- list()

  ## compute individual time specific DID
  for (i in 1:length(id_time_use)) {
    if (id_time_use[i] >= 3) {
      ## subset the data
      dat_use <- dat_panel %>%
        filter(!!sym(get_id(dat_panel)) %in% id_subj_use[[i]]) %>%
        filter(!!sym(get_wave(dat_panel)) <= id_time_use[i]) %>%
        filter(!!sym(get_wave(dat_panel)) >= id_time_use[i]-2) %>%
        mutate(outcome_diff = c(NA, diff(!!sym(outcome)))) %>%
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
