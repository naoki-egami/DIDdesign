

#' Implement the Double DID under the standard design
#' @keywords internal
did_std <- function(
  formula, data, id_unit, id_time, is_panel = TRUE, option
) {

  ## input check
  if (isFALSE(is_panel)) id_unit <- NULL
  if (isTRUE(is_panel) && is.null(id_unit)) {
    stop("A vaiable name should be provided to id_unit.")
  }

  ## prepare formulas
  fm_prep <- did_formula(formula, is_panel)

  ##
  ## handle cluster variable
  ##
  var_cluster <- option$id_cluster
  if (is.null(var_cluster) && isTRUE(is_panel)) {
    var_cluster <- "id_unit"
    var_cluster_pre <- id_unit
  }

  ##
  ## transform data
  ##
  if (isTRUE(is_panel)) {
    dat_did <- did_panel_data(
      fm_prep$var_outcome, fm_prep$var_treat, fm_prep$var_covars,
      var_cluster_pre, id_unit, id_time, data
    )
  } else {
    dat_did <- did_rcs_data(
      fm_prep$var_outcome, fm_prep$var_treat, fm_prep$var_post,
      fm_prep$var_covars, var_cluster,id_time, data
    )
  }

  ## -------------------------------- ##
  ## estimate ATT via DID and sDID    ##
  ## -------------------------------- ##

  ## point estimate
  fit_did  <- ddid_fit(fm_prep$fm_did, dat_did, lead = option$lead)
  weights  <- did_compute_weights(
    fm_prep, dat_did, var_cluster, is_panel, option
  )

  estimates <- double_did_compute(fit_did, weights, option$se_boot)

  ## compute double did estimate
  return(list(estimates = estimates, weights = weights))
}

#' Fit DID and sDID
#'
#' @param data An output of \code{did_panel_data()}.
#' @param formula A formula of the form \code{y ~ Gi + It + Gi * It + x1 + x2}.
#' @importFrom purrr map_dbl
#' @importFrom dplyr filter
#' @importFrom rlang .data
#' @return A numeric vector that containts the DID estimate and the sequential DID estimate.
#' @keywords internal
ddid_fit <- function(formula, data, lead = 0) {
  time_use <- c(-1, lead)
  est <- map_dbl(formula, function(fm) {
    fit <- lm(fm, data = filter(data, .data$id_time_std %in% time_use))
    return(fit$coef['Gi:It'])
  })

  return(est)
}


#' Compute Weighting Matrix via Bootstrap
#' @importFrom future.apply future_lapply
#' @importFrom future plan multicore sequential
#' @keywords internal
did_compute_weights <- function(
  fm_prep, dat_did, var_cluster, is_panel, option
) {


  ## setup cluster ID for bootstrap
  if (is.null(var_cluster)) {
    id_cluster_vec <- 1:nrow(dat_did)
  } else {
    id_cluster_vec <- pull(dat_did, !!sym(var_cluster)) %>% unique()
  }

  ##
  ## bootstrap
  ##

  ## setup worker
  if (isTRUE(option$parallel)) {
    plan(multicore)
  } else {
    plan(sequential)
  }

  ## point estimates
  ## use future_lapply to implement the bootstrap parallel
  est_boot <- do.call(rbind, future_lapply(1:option$n_boot, function(i) {
    ddid_boot(
      fm_prep, dat_did, id_cluster_vec, var_cluster, is_panel, option$lead
    )
  }, future.seed = TRUE))

  ## compute weights
  W <- cov(est_boot)
  w_did  <- (W[1,1] - W[1,2]) / (W[1,1] + W[2,2] - 2 * W[1,2])
  w_sdid <- (W[2,2] - W[1,2]) / (W[1,1] + W[2,2] - 2 * W[1,2])

  return(list(
    W = solve(W), weight_did = w_did, weight_sdid = w_sdid,
    Vcov = W
  ))

}

#' Bootstrap
#' @keywords internal
ddid_boot <- function(
  fm_prep, dat_did, id_cluster_vec, var_cluster, is_panel, lead
) {
  ## sample index
  id_boot <- sample(id_cluster_vec,
    size = length(id_cluster_vec), replace = TRUE
  )

  ## create dataset
  dat_tmp <- list()
  for (j in 1:length(id_boot)) {
    if (is.null(var_cluster)) {
      id_tmp <- id_boot[j]
    } else {
      id_tmp <- which(dat_did[,var_cluster] == id_boot[j])
    }
    dat_tmp[[j]] <- dat_did[id_tmp, ]
    dat_tmp[[j]]$id_unit <- j
  }

  ## create did_data object
  if (isTRUE(is_panel)) {
    dat_boot <- did_panel_data(
      var_outcome = "outcome", var_treat = 'treatment', fm_prep$var_covars,
      var_cluster, id_unit = "id_unit", id_time = 'id_time',
      data = do.call(rbind, dat_tmp)
    )
  } else {
    dat_boot <- did_rcs_data(
      var_outcome = "outcome", var_treat = "Gi", var_post = "It",
      fm_prep$var_covars, var_cluster,
      id_time = "id_time", do.call(rbind, dat_tmp)
    )
  }

  ## fit DID and sDID
  est <- ddid_fit(fm_prep$fm_did, dat_boot, lead)
  return(est)
}


#' Compute Point and Variance Estimates
#' @keywords internal
#' @importFrom dplyr tibble
double_did_compute <- function(fit, weights, se_boot = TRUE) {
  ## point estimate
  ## τ[d-did] = w[did] * τ[did] + w[s-did] * τ[s-did]
  ddid <- fit[1] * weights$weight_did +
          fit[2] * weights$weight_sdid

  ## variance estimate
  var_did  <- weights$Vcov[1,1]
  var_sdid <- weights$Vcov[2,2]
  if (isTRUE(se_boot)) {
    ## V(τ̂) = w^2[did] * Var(τ̂[did]) + w^2[s-did] * Var(τ̂[did]) +
    ##            2 * w[did] * w[s-did] * Cov(τ̂[did], τ̂[s-did])
    cov_did  <- weights$Vcov[1,2]
    ddid_var <- weights$weight_did^2 * var_did +
                weights$weight_sdid^2 * var_sdid +
                2 * cov_did * weights$weight_did * weights$weight_sdid
  } else {
    ## asymptotic variance
    ddid_var <- (1 / sum(weights$W))
  }

  ## organize output
  est <- tibble(
    estimator = c("Double-DID", "DID", "sDID"),
    estimate = c(ddid, fit[1], fit[2]),
    std.error = c(sqrt(ddid_var), sqrt(var_did), sqrt(var_sdid))
  )

  return(est)
}
