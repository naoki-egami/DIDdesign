
#' Generate staggered adoption design
#' @param n_obs Number of units.
#' @param n_time Number of time periods.
#' @return \code{simulate_sad()} returns a list containing the following elements:
#' \itemize{
#'   \item{dat} A long format dataset generated under the staggered adoption.
#'   \item{Amat} A matrix of treatment patterns.
#' }
#' @keywords internal
simulate_sad <- function(n_obs, n_times) {

  ## draw treatment timing from t = 2, ..., n_times + 1
  ## t = 1 is left as a baseline year where everyone is under control
  treatment_timing <- sort(sample.int(n_times, size = n_obs, replace = TRUE) + 1)

  ## create A matrix (time-varying treatment indicator)
  Amat <- matrix(0, nrow = n_obs, ncol = n_times)
  for (i in 1:n_obs) {
    if (treatment_timing[i] <= n_times) {
      Amat[i,treatment_timing[i]:n_times] <- 1
    }
  }

  ## create unit and time index
  id_subject <- rep(1:n_obs, n_times)
  id_time    <- rep(1:n_times, each = n_obs)
  Avec       <- as.vector(Amat)

  ## generate the outcome
  ## 1. time fixed effect
  delta_t <- rnorm(n_times)

  ## 2. unit fixed effect
  alpha_i <- rnorm(n_obs)

  ## 3. random error
  epsilon <- rnorm(n_times * n_obs)

  ## 4. generate Y
  Yvec <- alpha_i[id_subject] + delta_t[id_time] +
  0.2 * Avec + epsilon

  ## form a panel data
  dat_panel <- data.frame(
      outcome = Yvec, treatment = Avec,
      id_subject = id_subject, id_time = id_time
  )

  # dat_panel <- panelr::panel_data(dat_panel, id = id_subject, wave = id_time)
  return(list(dat = dat_panel, Amat = Amat))
}
