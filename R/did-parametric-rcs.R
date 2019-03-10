

#' Parametric Double Difference-in-Differences Estimator with Repeated Cross-Section Data
#' @param data A \code{diddesign_data} object.
#' @export
did_parametric_rcs <- function(data) {
  ## input checks
  if (!('diddesign_data' %in% class(data))) {
    stop("diddesign_data class object should be provided as data.")
  }

  t_pre <- length(data[[1]]$formula)
  if (t_pre <= 1) {
    stop("We reuqire more than two pre-treatment periods.\n")
  }

  m_vec <- 1:t_pre

  n_post <- length(data)
  result <- list()

  # ********************************************************* #
  #                                                           #
  #         model selection by GMM J-stats or t-test          #
  #         - use only pre-treatment data                     #
  #                                                           #
  # ********************************************************* #
  select_tmp <- rcs_selection(data[[1]]$pdata, data[[1]]$formula, attr(data[[1]], 'post_treat'))
  min_model  <- select_tmp$min_model
  attr(result, 'selection')  <- select_tmp[c('test_theta', 'test_se', 'min_model')]

  for (tt in 1:n_post) {

    dat_use    <- data[[tt]]$pdata
    fm_list    <- data[[tt]]$formula

    # ********************************************************* #
    #                                                           #
    #       run fixed effect model: get demeaned data           #
    #                                                           #
    # ********************************************************* #

    ## fit individual twoway fixed effect model
    ## get demeand matrix and response
    fit <- dat_trans <- list()
    for (ff in 1:length(fm_list)) {
      fit[[ff]]       <- lm(fm_list[[ff]], data = dat_use)
      y_var_name      <- all.vars(fm_list[[ff]])[1]
      dat_trans[[ff]] <- getX_rcs(fit[[ff]], dat_use[,y_var_name])
    }

    tmp <- tmp_min <- list()
    for (m in m_vec) {
      use_moments <- m_vec[m:length(m_vec)]

      ## point estimate
      est <- didgmmT_parametric_rcs(dat_trans[use_moments])
      tmp[[m]] <- est

      ## compute variance
      att_var <- cugmm_var_rcs(par = est$ATT, dat = est$data)

      ## compute Ci
      tmp_ci95     <- c(est$ATT - 1.96 * sqrt(att_var), est$ATT + 1.96 * sqrt(att_var))
      tmp_ci90     <- c(est$ATT - 1.64 * sqrt(att_var), est$ATT + 1.64 * sqrt(att_var))
      tmp_min[[m]] <- list("boot_est" = NULL, 'ci95' = tmp_ci95, 'ci90' = tmp_ci90)

    }

    ci95 <- tmp_min[[min_model]]$ci95
    ci90 <- tmp_min[[min_model]]$ci90

    # ********************************************************* #
    # standard did
    # ********************************************************* #
    # RETURN TWO-WAY FIXED EFFECT ESTIMATE
    did_est <- fit[[1]]$coef['treatment:post']

    ## variance
    did_var <- summary(fit[[use_moments]])$coef['treatment:post', 2]^2

    ## confidence interval
    tmp_se95 <- c(did_est + qnorm(0.025) * sqrt(did_var), did_est + qnorm(1 - 0.025) * sqrt(did_var))
    tmp_se90 <- c(did_est + qnorm(0.050) * sqrt(did_var), did_est + qnorm(1 - 0.050) * sqrt(did_var))

    did_boot_list <- list('boot_est' = NULL, 'ci95' = tmp_se95, 'ci90' = tmp_se90)
    did_save <- list("ATT" = did_est, 'results_variance' = did_boot_list)

    # ********************************************************* #
    #                                                           #
    #                     save results                          #
    #                                                           #
    # ********************************************************* #

    result[[tt]] <- list(
      'results_estimates' = tmp,
      'results_variance' = tmp_min,
      'results_standardDiD' = did_save,
      'min_model' = min_model,
      'ATT' = tmp[[1]]$ATT,
      'ci95' = ci95,
      'ci90' = ci90
    )
    attr(result[[tt]], 'post_treat') <- attr(data[[tt]], 'post_treat')
    attr(result[[tt]], 'method') <- 'parametric'
  }


  class(result) <- "diddesign"
  return(result)
}


#' get design matrix and outcome vector
#' @param lm_fit Fitted output from \code{\link{lm}}.
#' @keywords internal
getX_rcs <- function(lm_fit, response_orginal) {
  is_na <- is.na(response_orginal)
  x_tmp <- model.matrix(lm_fit)
  covariates <- x_tmp[, !(colnames(x_tmp) %in% c("treatment:post"))]
  treatment  <- x_tmp[, "treatment:post"]
  outcome    <- as.vector(model.frame(lm_fit)[,1])

  return(list("X" = covariates, 'y' = outcome, 'D' = treatment, 'is_na' = is_na))
}

#' CU-GMM function for RCS data
#' @param A list of dataset. Output of \code{link{getX_rcs}}.
#' @param Xdesign a list of design matrix. The variable of interest is named 'treatment:post'
#' @keywords internal
didgmmT_parametric_rcs <- function(dat, par_init = NULL) {


  ## 1. residualize outcome and treatment
  ##   - Regress Y on X
  ##   - Regress D on X
  for (d in 1:length(dat)) {
    outcome    <- dat[[d]]$y
    treatment  <- dat[[d]]$D
    covariates <- dat[[d]]$X ## this includes intercept term

    dat[[d]]$y_resid <- lm(outcome ~ covariates - 1)$residuals
    dat[[d]]$d_resid <- lm(treatment ~ covariates - 1)$residuals
  }


  ## 2. use resideualized outcome to estimate AT
  ## initialize parameter (ATT) if not given
  if(is.null(par_init)) par_init <- runif(1)
  ##   - moment conditions: E[D(Y- Db)] = 0
  est <- optim(par = par_init, fn = cugmm_loss_rcs, method = "BFGS", dat = dat)

  return(list('est' = est, 'ATT' = est$par, 'data' = dat))
}


#' CU-GMM loss function for RCS data
#' @param par parameter.
#' @param dat a list of data consists of y_resid and d_reisd.
#' @keywords internal
cugmm_loss_rcs <- function(par, dat) {
  n_moment_condition <- length(dat)

  XXe <- matrix(NA, nrow = length(dat[[1]]$is_na), ncol = n_moment_condition)
  for (j in 1:n_moment_condition) {
    is_na <- dat[[j]]$is_na
    loss  <- dat[[j]]$d_resid * (dat[[j]]$y_resid - dat[[j]]$d_resid* par)
    XXe[!is_na, j] <- as.vector(loss)
  }

  gbar  <- colMeans(XXe, na.rm = TRUE)
  XXe   <- na.omit(XXe)
  Omega <- (t(XXe) %*% XXe)
  loss  <- as.vector(t(gbar) %*% solve(Omega / nrow(XXe), gbar))
  # cat("loss = ", loss, '\n')
  return(loss)
}


#' CU-GMM variance function for RCS data
#' @param par parameter.
#' @param dat a list of data consists of y_resid and d_reisd.
#' @keywords internal
cugmm_var_rcs <- function(par, dat) {
  n_moment_condition <- length(dat)

  G <- matrix(NA, nrow = n_moment_condition, ncol = 1)
  XXe <- matrix(NA, nrow = length(dat[[1]]$is_na), ncol = n_moment_condition)
  for (j in 1:n_moment_condition) {
    is_na <- dat[[j]]$is_na
    loss  <- dat[[j]]$d_resid * (dat[[j]]$y_resid - dat[[j]]$d_resid * par)
    XXe[!is_na, j] <- as.vector(loss)

    G[j,1]   <- mean(dat[[j]]$d_resid^2)
  }

  XXe   <- na.omit(XXe)
  Omega <- (t(XXe) %*% XXe) / nrow(XXe)

  gmm_var <- as.vector(solve(t(G) %*%  solve(Omega, G)))
  att_var <- gmm_var / length(dat[[1]]$is_na)

  return(att_var)
}

#
# #' Generate dummy variables
# #' @param data An element of \code{diddesign_data} object.
# #' @return A data matrix. The first column is ones' if \code{add_intercept = TRUE}.
# #' @keywords internal
# gen_dummies <- function(data, add_intercept = TRUE) {
#   pdata <- data$pdata
#
#   ## generate group dummies
#   ## we can use id_subject but make sure it is at the treatment group levels
#   group_dummy <- model.matrix(~ factor(treatment) - 1, data = pdata)
#   ## generate time dummies
#   time_dummy  <- model.matrix(~ factor(id_time) - 1, data = pdata)
#   ## generate interaction dummies
#   interact_dummy <- model.matrix(~ factor(treatment):factor(id_time) - 1, data = pdata)
#
#   if(isTRUE(add_intercept)) {
#     ones <- matrix(1, nrow = nrow(pdata), ncol = 1)
#     Xmat <- list('ones' = ones, 'group' = data.matrix(group_dummy),
#     'time'  = data.matrix(time_dummy), 'interaction' = data.matrix(interact_dummy))
#   } else {
#     Xmat <- cbind(group_dummy, time_dummy, interact_dummy)
#   }
#
#   return(Xmat)
# }
#
#
#
# get_diffmat <- function(time_pre) {
#   time_total <- time_pre + 1
#   diff_mat <- matrix(0, nrow = time_pre, ncol = 2 * time_total)
#   counter <- 0
#   for (tt in 1:time_pre) {
#     diff_mat[tt,1+counter] <--1
#     diff_mat[tt,time_total+1+counter] <- 1
#     counter <- counter + 1
#   }
#   return(diff_mat)
# }
