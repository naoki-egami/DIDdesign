

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

  ## point estimate for DID and sDID
  fit_did <- lapply(option$lead, function(ll) ddid_fit(fm_prep$fm_did, dat_did, lead = ll))

  ## compute weights via bootstrap --> compute W
  weights  <- did_compute_weights(fm_prep, dat_did, var_cluster, is_panel, option)

  ## compute double did estimate and variance
  estimates <- double_did_compute(fit_did, weights, option$lead, option$se_boot)

  ## compute double did estimate
  return(estimates)
}

#' Fit DID and sDID
#'
#' @param data An output of \code{did_panel_data()}.
#' @param formula A formula of the form \code{y ~ Gi + It + Gi * It + x1 + x2}.
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
#' @importFrom future.apply future_lapply
#' @importFrom future plan multicore sequential
#' @importFrom purrr map
#' @importFrom dplyr pull
#' @importFrom rlang !! sym
#' @importFrom stats cov
#' @keywords internal
did_compute_weights <- function(
  fm_prep, dat_did, var_cluster, is_panel, option
) {

  ## setup cluster ID for bootstrap
  if (is.null(var_cluster)) {
    id_cluster_vec <- 1:nrow(dat_did)
  } else {
    id_cluster_vec <- unique(pull(dat_did, !!sym(var_cluster)))
  }

  ## --------------------------------------------
  ## bootstrap to compute weights
  ## --------------------------------------------

  ## setup worker
  if (isTRUE(option$parallel)) {
    plan(multicore)
  } else {
    plan(sequential)
  }

  ## use future_lapply to implement the bootstrap parallel
  est_boot <- future_lapply(1:option$n_boot, function(i) {
    tryCatch({
      ddid_boot(fm_prep, dat_did, id_cluster_vec, var_cluster, is_panel, option$lead)
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

  return(weights_save)
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
  est <- lapply(lead, function(ll) ddid_fit(fm_prep$fm_did, dat_boot, ll))
  return(est)
}


#' Compute Point and Variance Estimates
#' @keywords internal
#' @importFrom dplyr as_tibble bind_rows
double_did_compute <- function(fit, weights, lead, se_boot = FALSE) {

  estimates <- W_save <- vector("list", length = length(lead))
  for (ll in 1:length(lead)) {
    ## compute the
    ddid <- as.vector(weights[[ll]]$weights %*% fit[[ll]])

    ## variance estimate
    var_did  <- weights[[ll]]$vcov[1,1]
    var_sdid <- weights[[ll]]$vcov[2,2]

    if (isTRUE(se_boot)) {
      ## V(τ̂) = w^2[did] * Var(τ̂[did]) + w^2[s-did] * Var(τ̂[did]) +
      ##            2 * w[did] * w[s-did] * Cov(τ̂[did], τ̂[s-did])
      cov_did  <- weights[[ll]]$vcov[1,2]
      ddid_var <- weights[[ll]]$weights[1]^2 * var_did +
                  weights[[ll]]$weights[2]^2 * var_sdid +
                  2 * cov_did * prod(weights[[ll]]$weights)
    } else {
      ## asymptotic variance
      ddid_var <- (1 / sum(weights[[ll]]$W))
    }

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
