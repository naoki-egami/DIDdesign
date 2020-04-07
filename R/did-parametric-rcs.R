

#' Parametric Double Difference-in-Differences Estimator with Repeated Cross-Section Data
#' @param data A \code{diddesign_data} object.
#' @family estimation functions
#' @export
did_parametric_rcs <- function(data, se_boot = TRUE,
  n_boot = 500, boot_min = TRUE, id_cluster = NULL,
  only_last = TRUE, verbose = TRUE, alpha = alpha) {
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
  select_tmp <- rcs_selection(data[[1]]$pdata, data[[1]]$formula, attr(data[[1]], 'post_treat'), alpha = alpha)
  min_model  <- select_tmp$min_model
  attr(result, 'selection')  <- select_tmp[c('test_theta', 'test_se', 'min_model')]
  attr(result, 'sign')       <- select_tmp$sign

  for (tt in 1:n_post) {
    dat_use    <- data.frame(data[[tt]]$pdata)
    fm_list    <- data[[tt]]$formula
    id_time2   <- as.numeric(as.factor(dat_use$id_time))
    max_time2  <- max(id_time2)

    if (isTRUE(only_last)) {
      for (i in 0:(t_pre-2)) {
        dat_use[,paste("yd", i, sep = '')] <- ifelse(id_time2 >= (max_time2-1),
          dat_use[,paste("yd", i, sep = '')], NA)
      }
    }


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
      if (isTRUE(se_boot)) {
        boot_est <- var_boot_parametric_rcs(dat_trans[use_moments], n_boot, id_cluster, verbose)
        att_var  <- var(boot_est)
      } else {
        att_var  <- cugmm_var_rcs(par = est$ATT, dat = est$data)
        boot_est <- NULL
      }


      ## compute Ci
      tmp_ci95     <- c(est$ATT - 1.96 * sqrt(att_var), est$ATT + 1.96 * sqrt(att_var))
      tmp_ci90     <- c(est$ATT - 1.64 * sqrt(att_var), est$ATT + 1.64 * sqrt(att_var))
      tmp_min[[m]] <- list("boot_est" = boot_est, 'ci95' = tmp_ci95, 'ci90' = tmp_ci90, 'se' = sqrt(att_var))

    }

    ci95 <- tmp_min[[min_model]]$ci95
    ci90 <- tmp_min[[min_model]]$ci90
    se   <- tmp_min[[min_model]]$se

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

    did_boot_list <- list('boot_est' = NULL, 'ci95' = tmp_se95, 'ci90' = tmp_se90, 'se' = sqrt(did_var))
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
      'ATT' = tmp[[min_model]]$ATT,
      'ci95' = ci95,
      'ci90' = ci90,
      'se'   = se
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
  covariates <- x_tmp[, !(colnames(x_tmp) %in% c("treatment:post", "(Intercept)"))]
  treatment  <- x_tmp[, "treatment:post"]
  outcome    <- as.vector(model.frame(lm_fit)[,1])

  return(list("X" = covariates, 'y' = outcome, 'D' = treatment, 'is_na' = is_na))
}

#' CU-GMM function for RCS data
#' @param A list of dataset. Output of \code{link{getX_rcs}}.
#' @param Xdesign a list of design matrix. The variable of interest is named 'treatment:post'
#' @keywords internal
didgmmT_parametric_rcs <- function(dat, par_init = NULL, optim = FALSE) {


  ## 1. residualize outcome and treatment
  ##   - Regress Y on X
  ##   - Regress D on X
  for (d in 1:length(dat)) {
    outcome    <- dat[[d]]$y
    treatment  <- dat[[d]]$D
    covariates <- dat[[d]]$X
    # residualize the outcome and the treatment
    dat[[d]]$y_resid <- lm(outcome ~ covariates)$residuals
    dat[[d]]$d_resid <- lm(treatment ~ covariates)$residuals

  }


  ## 2. use resideualized outcome to estimate AT
  ## initialize parameter (ATT) if not given
  if(is.null(par_init)) par_init <- runif(1)
  ##   - moment conditions: E[D(Y- Db)] = 0
  ## two-step GMM
  if(isTRUE(optim)) {
    # estimate by optim
    # cat("1st stage\n")
    est1st <- optim(par = par_init, fn = cugmm_loss_rcs_init, method = "BFGS", dat = dat)

    # cat("2nd stage\n")
    est2nd <- optim(par = est1st$par, fn = cugmm_loss_rcs, method = "BFGS", dat = dat, init = est1st$par)
  } else {
    # analytical solution
    # cat("1st stage") W = I
    b_bar <- mean(dat[[1]]$d_resid^2)
    a_bar <- rep(NA, length(dat))
    for (d in 1:length(dat)) {
      a_bar[d] <- mean((dat[[d]]$y_resid * dat[[d]]$d_resid)^2)
    }

    I      <- diag(length(dat))
    est1st <- solve_att_rcs(a_bar = a_bar, b_bar = b_bar, W = I)

    ## compute  the new weighting matrix
    g <- matrix(NA, nrow = length(dat), ncol = length(dat[[1]]$d_resid))
    for (d in 1:length(dat)) {
      g[d, ] <- (dat[[d]]$y_resid - dat[[d]]$d_resid * est1st) * dat[[d]]$d_resid
    }
    W <- solve(g %*% t(g) / length(dat[[1]]$d_resid))
    est2nd <- list(par = solve_att_rcs(a_bar = a_bar, b_bar = b_bar, W = W))
  }




  return(list('est' = est2nd, 'ATT' = est2nd$par, 'data' = dat))
}


#' Residualize
# residualize <- function(y, X) {
#   X1 <- cbind(rep(1, nrow(X)), X)
#   I  <- diag(nrow(X))
#   resid  <- y - X1 %*% solve(t(X1) %*% X1, t(X1) %*% y)
#   return(as.vector(resid))
# }

#' Analytical solution to ATT
solve_att_rcs <- function(a_bar, b_bar, W) {
  a_bar <- matrix(a_bar)
  num   <- sum(t(a_bar) %*% W)
  dem   <- b_bar * sum(W)

  ## solve
  beta_est <- num / dem
  return(beta_est)
}

#' CU-GMM loss function for RCS data
#' @param par parameter.
#' @param dat a list of data consists of y_resid and d_reisd.
#' @keywords internal
cugmm_loss_rcs <- function(par, dat, init = NULL) {
  if (is.null(init)) init <- par

  n_moment_condition <- length(dat)

  XXe <- XXe0 <- matrix(NA, nrow = length(dat[[1]]$is_na), ncol = n_moment_condition)
  for (j in 1:n_moment_condition) {
    is_na <- dat[[j]]$is_na
    loss  <- dat[[j]]$d_resid * (dat[[j]]$y_resid - dat[[j]]$d_resid * par)
    loss0  <- dat[[j]]$d_resid * (dat[[j]]$y_resid - dat[[j]]$d_resid * init)
    XXe[!is_na, j] <- as.vector(loss)
    XXe0[!is_na, j] <- as.vector(loss0)
  }

  gbar  <- colMeans(XXe, na.rm = TRUE)
  XXe   <- na.omit(XXe); XXe0 <- na.omit(XXe0)
  # Omega <- (t(XXe) %*% XXe)
  Omega <- (t(XXe0) %*% XXe0)
  loss  <- as.vector(t(gbar) %*% solve(Omega / nrow(XXe), gbar))
  return(loss)
}


#' CU-GMM loss function for RCS data
#' @param par parameter.
#' @param dat a list of data consists of y_resid and d_reisd.
#' @keywords internal
cugmm_loss_rcs_init <- function(par, dat) {
  n_moment_condition <- length(dat)

  XXe <- matrix(NA, nrow = length(dat[[1]]$is_na), ncol = n_moment_condition)
  for (j in 1:n_moment_condition) {
    is_na <- dat[[j]]$is_na
    loss  <- dat[[j]]$d_resid * (dat[[j]]$y_resid - dat[[j]]$d_resid * par)
    XXe[!is_na, j] <- as.vector(loss)
  }

  gbar  <- colMeans(XXe, na.rm = TRUE)
  XXe   <- na.omit(XXe)
  Omega <- diag(ncol(XXe)) # (t(XXe) %*% XXe)
  loss  <- as.vector(t(gbar) %*% gbar) ## as.vector(t(gbar) %*% solve(Omega / nrow(XXe), gbar))
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
  attr(att_var, "W") <- solve(Omega)
  return(att_var)
}


#' Varinace for Repated Cross-Section with Block-Bootstrap
#' @param data a list of data
#' @param n_boot number of bootstrap iterations
#' @param id_cluster cluster index at the level of blocking
#' @keywords internal
var_boot_parametric_rcs <- function(dat, n_boot, id_cluster, verbose) {
  iter     <- 1
  is_na    <- is.na(dat[[1]]$y)
  n_obs    <- length(na.omit(dat[[1]]$y))
  id_full  <- 1:n_obs
  boot_out <- rep(NA, n_boot)
  if (!is.null(id_cluster)) id_cluster <- id_cluster[!is_na]

  iter_show <- round(n_boot * 0.1)
  while(iter <= n_boot) {
    tryCatch({
      ## sample bootstrap
      if (is.null(id_cluster)) {
        id_boot <- sample(id_full, size = n_obs, replace = TRUE)
      } else {
        id_cluster_boot <- sample(unique(id_cluster), size = length(unique(id_cluster)), replace = TRUE)
        id_boot <- NULL
        for (i in 1:length(id_cluster_boot)) {
          id_boot <- c(id_boot, id_full[id_cluster == id_cluster_boot[i]])
        }
      }
      ## 1. residualize outcome and treatment
      ##   - Regress Y on X
      ##   - Regress D on X
      dat_boot <- list()
      for (d in 1:length(dat)) {
        ## extract info
        outcome    <- dat[[d]]$y[!is_na]; outcome     <- outcome[id_boot]
        treatment  <- dat[[d]]$D[!is_na]; treatment   <- treatment[id_boot]
        covariates <- dat[[d]]$X[!is_na,]; covariates <- covariates[id_boot,]
        # residualize the outcome and the treatment
        dat_boot[[d]] <- list()
        dat_boot[[d]]$y_resid <- lm(outcome ~ covariates)$residuals
        dat_boot[[d]]$d_resid <- lm(treatment ~ covariates)$residuals
      }

      b_bar <- mean(dat_boot[[1]]$d_resid^2)
      a_bar <- rep(NA, length(dat))
      for (d in 1:length(dat)) {
        a_bar[d] <- mean((dat_boot[[d]]$y_resid * dat_boot[[d]]$d_resid)^2)
      }

      I      <- diag(length(dat))
      est1st <- solve_att_rcs(a_bar = a_bar, b_bar = b_bar, W = I)

      ## compute  the new weighting matrix
      g <- matrix(NA, nrow = length(dat), ncol = length(dat_boot[[1]]$d_resid))
      for (d in 1:length(dat)) {
        g[d, ] <- (dat_boot[[d]]$y_resid - dat_boot[[d]]$d_resid * est1st) * dat_boot[[d]]$d_resid
      }
      W <- solve(g %*% t(g) / length(dat_boot[[1]]$d_resid))
      boot_out[iter] <- solve_att_rcs(a_bar = a_bar, b_bar = b_bar, W = W)

      ## update iter
      iter <- iter + 1

    }, error = function(e) {
      NULL
    })


    if (isTRUE(verbose)) {
      if ((iter %% iter_show) == 0) {
          cat('\r', iter, "out of", n_boot, "bootstrap iterations")
          flush.console()
      }
    }

  }

  cat("\n")

  return(boot_out)
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
