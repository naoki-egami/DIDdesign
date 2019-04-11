

##
## set of function to simulate data from k-th order parallel trend
##


#' Generate Synthetic Data from k-th Order Parallel Trend
#'
#' @param n_unit the number of units.
#' @param n_period the number of time periods.
#'  The last time point is always used as the post-period. Maximal value is \code{5L}.
#' @param k order of trend.
#' @param dgp either \code{"parametric"} or \code{"nonparametric"}.
#' @param panel boolean. If \code{TRUE} data set are generated in panel format,
#'  otherwise they are in repeated cross-section form.
#' @param prop_treat proportion of treated unit. Default is 0.5.
#' @export
did_simulate <- function(n_unit, n_period = 3L, k = 1L,
  dgp = 'parametric', panel = TRUE, prop_treat = 0.5
) {

  ##
  ## input check
  ##
  if (!is.integer(k)) { stop("k should be integer type") }
  if (!is.integer(n_period)) { stop("n_period should be integer type.") }
  if (k <= 0L) { stop("k should take positive value.") }
  if (n_period <= 1) { stop("n_period should be greater than 2.") }
  if ((n_period - 1) < k) { stop("k cannot exceed n_period - 1") }
  if (n_period > 5L) stop("we currently support n_period up to 5.")

  ##
  ##
  ##
  if (isTRUE(panel)) {

  } else {

  }

}

#' helper function
#' @keywords internal
sim_np_panel <- function(n1, n0, n_period, k) {

}

#' @keywords internal
sim_pm_rcs <- function(n1, n0, n_period, k) {
  t_trend <- runif(n_period, -1, 1)
  mu0 <- t_trend

  Dmat   <- getDtf(n_period, k-1)
  n_base <- ncol(Dmat) - nrow(Dmat)
  delta0 <- runif(nrow(Dmat), -1, 1)
  window <- (nrow(Dmat)-n_base-1):nrow(Dmat)
  for (i in 1:nrow(Dmat)) {
    Dmat[nrow(Dmat)-(i-1), window]
    window <- window - i
  }
}
