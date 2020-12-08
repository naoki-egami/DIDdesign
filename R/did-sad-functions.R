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
get_periods <- function(Gmat, thres = 3) {
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
