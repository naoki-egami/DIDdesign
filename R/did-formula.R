
#' Transform formula 
#' @keywords internal
#' @import Formula
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
