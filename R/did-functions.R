#' Set default value of options
#' @param option A list of options.
#' @return  A list of options.
#' @keywords internal
set_option <- function(option) {

  if (!exists('n_boot', option)) option$n_boot         <- 30
  if (!exists('parallel', option)) option$parallel     <- TRUE
  if (!exists('se_boot', option)) option$se_boot       <- FALSE
  if (!exists('id_cluster', option)) option$id_cluster <- NULL
  if (!exists('lead', option)) option$lead             <- 0
  if (!exists('thres', option)) option$thres           <- 2
  if (!exists('lag', option)) option$lag               <- 1
  if (!exists('stdz', option)) option$stdz             <- TRUE
  return(option)
}

#' Transform formula
#' @param formula Formula.
#' @param is_panel A boolean argument to indicate if data is in the panel format.
#' @return A list of working objects.
#' @keywords internal
#' @import Formula
#' @importFrom stats update
did_formula <- function(formula, is_panel) {

  ## obtain variable names
  fm <- as.Formula(formula)
  n_vars_fm <- length(fm)

  ## handle outcome
  if (n_vars_fm[1] != 1) stop("Invalid formula (lhs)")
  var_outcome <- all.vars(formula(fm, lhs = 1, rhs = 0))

  ## handle covariates
  if (n_vars_fm[2] == 2) {
    var_covars <- all.vars(formula(fm, lhs = 0, rhs = 2))
    fm_covar <- formula(fm, lhs = 0, rhs = 2)
  } else if (n_vars_fm[2] == 1){
    var_covars <- fm_covar <- NULL
  } else {
    stop("Invalid formula (rhs)")
  }

  ## handle treatment & post-period indicator
  if (isTRUE(is_panel)) {
    var_treat <- all.vars(formula(fm, lhs = 0, rhs = 1))
    var_post  <- NULL
  } else {
    var_tmp <- all.vars(formula(fm, lhs = 0, rhs = 1))
    var_treat <- var_tmp[1]
    var_post  <- var_tmp[2]
  }

  ##
  ## update formula
  ##
  if (is.null(var_covars)) {
    fm_did <- list(
      as.formula(outcome ~ Gi + It + Gi : It),
      as.formula(outcome_delta ~ Gi + It + Gi : It)
    )
  } else {
    fm_did <- list(
      update(fm_covar, outcome ~ Gi + It + Gi : It + .),
      update(fm_covar, outcome_delta ~ Gi + It + Gi : It + .)
    )
  }

  return(list(
    fm_did = fm_did, fm_covar = fm_covar,
    var_outcome = var_outcome, var_treat = var_treat,
    var_post = var_post, var_covars = var_covars
  ))

}


#' Setup Parallel Workers
#' @keywords internal
#' @param is_parallel A boolean argument. If \code{TRUE} bootstrap will be conducted on multiple cores.
setup_parallel <- function(is_parallel) {
  if (isTRUE(is_parallel) & parallelly::availableCores(constraints = "multicore") > 1) {
    future::plan(future::multicore)
  } else if (isTRUE(is_parallel) & parallelly::availableCores() > 1) {
    future::plan(future::multisession)
  } else {
    future::plan(future::sequential)
  }
}
