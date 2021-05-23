

#' Implement the Double DID under the standard design
#' @param formula Formula.
#' @param data Panel data.
#' @param id_unit A variable name for unit index. Required for panel data.
#' @param id_time A variable name for time index.
#' @param is_panel A boolean argument. Should be \code{FALSE} for the repeated cross-sectional data.
#' @param option A list of option parameters.
#' @return Double DID estimates.
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
      fm_prep$var_covars, var_cluster, id_time, data
    )
  }

  ## -------------------------------- ##
  ## estimate ATT via DID and sDID    ##
  ## -------------------------------- ##

  ## point estimate for DID and sDID
  fit_did <- lapply(option$lead, function(ll) ddid_fit(fm_prep$fm_did, dat_did, lead = ll))

  ## compute weights via bootstrap --> compute W
  boot_out  <- did_compute_weights(fm_prep, dat_did, var_cluster, is_panel, option)

  ## compute double did estimate and variance
  estimates <- double_did_compute(fit_did, boot_out, option$lead, option$se_boot)

  ## compute double did estimate
  return(estimates)
}

#' Fit DID and sDID
#'
#' @param data An output of \code{did_panel_data()}.
#' @param formula A formula of the form \code{y ~ Gi + It + Gi * It + x1 + x2}.
#' @param lead A value of the lead parameter.
#' @importFrom purrr map_dbl
#' @importFrom dplyr filter
#' @importFrom rlang .data
#' @importFrom stats lm
#' @return A numeric vector that contains the DID estimate and the sequential DID estimate.
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
#' @param formula Formula.
#' @param dat_did Panel data.
#' @param var_cluster A variable used for clustering standard errors.
#' @param is_panel A boolean argument to indicate if data is in the panel format or the repeated cross-section format.
#' @param option A list of options.
#' @return A list of double DID weights.
#' @importFrom future.apply future_lapply
#' @importFrom future plan multicore sequential
#' @importFrom purrr map
#' @importFrom dplyr pull
#' @importFrom rlang !! sym
#' @importFrom stats cov
#' @keywords internal
did_compute_weights <- function(
  formula, dat_did, var_cluster, is_panel, option
) {

  ## setup cluster ID for bootstrap
  if (is.null(var_cluster)) {
    id_cluster_vec <- 1:nrow(dat_did)
  } else {
    id_cluster_vec <- unique(pull(dat_did, !!sym(var_cluster)))
  }

  ## --------------------------------------------
  ## Compute weights
  ## --------------------------------------------

  if (isTRUE(option$se_boot)) {
    ##
    ## bootstrap to compute weights
    ##

    ## setup worker
    if (isTRUE(option$parallel)) {
      plan(multicore)
    } else {
      plan(sequential)
    }

    ## use future_lapply to implement the bootstrap parallel
    est_boot <- future_lapply(1:option$n_boot, function(i) {
      tryCatch({
        ddid_boot(formula, dat_did, id_cluster_vec, var_cluster, is_panel, option$lead)
      }, error = function(e) {
        NULL
      })
    }, future.seed = TRUE)
    est_boot <- est_boot[lengths(est_boot) != 0]

    ## convert VCVO to weights
    weights_save <- vector("list", length = length(option$lead))
    for (ll in 1:length(option$lead)) {
      tmp <- do.call(rbind, map(est_boot, ~.x[[ll]]))
      W <- cov(tmp)
      w_did  <- (W[1,1] - W[1,2]) / (W[1,1] + W[2,2] - 2 * W[1,2])
      w_sdid <- (W[2,2] - W[1,2]) / (W[1,1] + W[2,2] - 2 * W[1,2])
      weights_save[[ll]] <- list(W = solve(W), vcov = W, weights = c(w_did, w_sdid))
    }
  } else {
    ##
    ## analytically compute weights
    ##
    weights_save <- vector("list", length = length(option$lead))
    for (ll in 1:length(option$lead)) {
      ## obtain residuals
      resid <- ddid_resid(formula$fm_did, dat_did, lead = option$lead[[ll]])
      ep_vcov  <- get_vcov(resid, dat_did, var_cluster = NULL)

      ## obtain design matrix
      X <- model.matrix(formula$fm_did[[1]], data = dat_did)

      ## compute variance covariance matrix
      XX_inv <- solve( t(X) %*% X )
      tmp <- XX_inv %*% ( t(X) %*% ep_vcov %*% X) %*% XX_inv

      ## variance covariance matrix of τ[DID] and τ[s-DID]
      W <- tmp[c(4, 8), c(4, 8)]
      w_did  <- (W[1,1] - W[1,2]) / (W[1,1] + W[2,2] - 2 * W[1,2])
      w_sdid <- (W[2,2] - W[1,2]) / (W[1,1] + W[2,2] - 2 * W[1,2])
      weights_save[[ll]] <- list(W = solve(W), vcov = W, weights = c(w_did, w_sdid))
    }

    ## compute variance covariance matrix

  }

  return(list(weights = weights_save, estimates = est_boot))
}

#' Bootstrap
#' @inheritParams did_compute_weights
#' @param lead A vector of lead parameters.
#' @return A list of DID and sDID estimates.
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
  est <- lapply(lead, function(ll) ddid_fit(fm_prep$fm_did, dat_boot, ll))
  return(est)
}


#' Compute Point and Variance Estimates
#' @param fit A lift of fitted objects.
#' @param boot An output from bootstrap.
#' @param lead A vector of lead parameters.
#' @param se_boot A boolean argument to indicate if standard errors are computed based on the bootstrap-based (i.e., empirical variance),
#'                 or based on the asymptotic approximation.
#' @return A list of estimates and weights.
#' @keywords internal
#' @importFrom dplyr as_tibble bind_rows
double_did_compute <- function(fit, boot, lead, se_boot) {

  weights  <- boot$weights
  boot_est <- boot$estimates
  estimates <- W_save <- vector("list", length = length(lead))
  for (ll in 1:length(lead)) {
    ## compute the
    ddid <- as.vector(weights[[ll]]$weights %*% fit[[ll]])

    ## variance estimate for DID and sDID
    var_did  <- weights[[ll]]$vcov[1,1]
    var_sdid <- weights[[ll]]$vcov[2,2]

    ## variance estimate for Double DID
    if (isTRUE(se_boot)) {
      ## bootstrap based variance
      ddid_boot <- purrr::map_dbl(boot_est, ~ as.vector(weights[[ll]]$weights %*% .x[[ll]]))
      ddid_var  <- var(ddid_boot, na.rm = TRUE)
    } else {
      ## asymptotic variance
      ddid_var <- (1 / sum(weights[[ll]]$W))
    }

    ## summarize estimates
    estimates[[ll]] <- data.frame(
      estimator   = c("Double-DID", "DID", "sDID"),
      lead        = lead[ll],
      estimate    = c(ddid, fit[[ll]]),
      std.error   = c(sqrt(ddid_var), sqrt(var_did), sqrt(var_sdid)),
      ddid_weight = c(NA, weights[[ll]]$weights)
    )

    W_save[[ll]] <- weights[[ll]]$W
  }

  estimates <- as_tibble(bind_rows(estimates))
  return(list(estimates = estimates, weights = W_save))
}
