

#' Parametric Double Difference-in-Differences Estimator
#' @param data data
#' @param se_boot A boolean argument.
#'  If set to \code{TRUE}, standard errors are computed by the block bootstrap.
#'  Bootstrap iterations are set by \code{n_boot} argument.
#' @param n_boot The number of bootstrap iterations. Required when \code{se_boot == TRUE}.
#' @param boot_min If \code{TRUE}, bootstrap is used only for the selected model.
#'  This option helps reduce computational burdens.
#' @param select Selection criteria used,
#'  one of "HQIC", "BIC", "tt1" (T-test), "tt2" (T-test with Bonferroni correction).
#'  The selected model is used to estimate bootstrap variance when \code{boot_min = TRUE}.
#' @importFrom plm plm vcovHC
#' @export
did_parametric <- function(data, se_boot = FALSE, n_boot = 1000, boot_min = TRUE,
                           select = "HQIC", est_did = FALSE, is_covariates
) {

  ## input checks
  if (!('diddesign_data' %in% class(data))) {
    stop("diddesign_data class object should be provided as data.")
  }
  t_pre <- ncol(data[[1]]$Y) - 1
  if (t_pre <= 1) {
    stop("We reuqire more than two pre-treatment periods.\n")
  }


  n_post <- length(data)
  result <- list()

  # ********************************************************* #
  #                                                           #
  #         model selection by GMM J-stats or t-test          #
  #         - use only pre-treatment data                     #
  #                                                           #
  # ********************************************************* #
  select_tmp <- fe_selection(data[[1]]$pdata, data[[1]]$formula, attr(data[[1]], 'post_treat'))
  HQIC       <- 0 #select_tmp$HQIC
  BIC        <- 0 #select_tmp$BIC
  min_model  <- select_tmp$min_model

  ## save 
  attr(result, 'selection')  <- select_tmp[c('test_theta', 'test_se', 'min_model')]

  for (tt in 1:n_post) {

    m_vec <- 1:t_pre
    # select_tmp <- gmm_selection(Y = data[[tt]]$Y, D = data[[tt]]$D,
    #                             mvec = m_vec, t_pre = t_pre, select = select, n_boot = n_boot)
    # HQIC       <- select_tmp$HQIC
    # BIC        <- select_tmp$BIC
    # min_model  <- select_tmp$min_model

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
      fit[[ff]]       <- plm(fm_list[[ff]], data = dat_use, model = 'within', effect = 'twoways')
      dat_trans[[ff]] <- getX(fm_list[[ff]], fit[[ff]], dat_use)
    }

    # ********************************************************* #
    #                                                           #
    #             GMM: point estimate and variance              #
    #                                                           #
    # ********************************************************* #
    ## for each methods
    ##  - M1: D0y ~ DTy
    ##  - M2: D1y ~ DTy
    cat("\n... estimating treatment effect for ", attr(data[[tt]], 'post_treat'), " ...\n")

    tmp <- tmp_min <- list()
    for (m in m_vec) {
      use_moments <- m_vec[m:length(m_vec)]

      ## initialize
      coef_vec <- lapply(fit[use_moments], function(x) x$coef)
      if (length(coef_vec[[1]]) == 1) {
        par_init <- mean(unlist(coef_vec))
      } else {
        par_init <- mean(sapply(coef_vec, function(x) x[1]))
      }

      ## estimate
      if (isTRUE(is_covariates)) {
        est <- didgmmT_parametric_resid(dat_trans[use_moments], dat_use$id_subject, par_init)
      } else {
        est <- didgmmT_parametric(dat_trans[use_moments], dat_use$id_subject, par_init)
      }

      tmp[[m]] <- est

      ## compute variance
      if (isTRUE(se_boot)) {
        boot_out <- list()
      } else {
        # cat("... computing GMM variance ...\n")

        ## compute asymptotic variance
        if (isTRUE(is_covariates)) {
          att_var <- cugmm_var_parametric_resid(par = est$ATT, dat = est$data,
                                                id_subject = dat_use$id_subject)
        } else {
          att_var <- cugmm_var_parametric(par = est$ATT, dat = dat_trans[use_moments],
                                          id_subject = dat_use$id_subject)
        }

        ## compute Ci
        tmp_ci95     <- c(est$ATT - 1.96 * sqrt(att_var), est$ATT + 1.96 * sqrt(att_var))
        tmp_ci90     <- c(est$ATT - 1.64 * sqrt(att_var), est$ATT + 1.64 * sqrt(att_var))
        tmp_min[[m]] <- list("boot_est" = NULL, 'ci95' = tmp_ci95, 'ci90' = tmp_ci90)
      }

    }

    ci95 <- tmp_min[[min_model]]$ci95
    ci90 <- tmp_min[[min_model]]$ci90

    # ********************************************************* #
    #                                                           #
    #                 standard did estimates                    #
    #                                                           #
    # ********************************************************* #
    if (isTRUE(est_did)) {
      cat("... computing the standard DiD estimate ...\n")

      ## point estimate
      did_est <- std_did(Y = data[[tt]]$Y, D = data[[tt]]$D)

      ## bootstrap
      # if (isTRUE(se_boot)) {
      tmp_est <- std_did_boot(Y = data[[tt]]$Y, D = data[[tt]]$D, n_boot = n_boot)
      tmp_se95 <- quantile(tmp_est, prob = c(0.025, 0.975))
      tmp_se90 <- quantile(tmp_est, prob = c(0.05, 0.95))
      did_boot_list <- list('boot_est' = tmp_est, 'ci95' = tmp_se95, 'ci90' = tmp_se90)
      # } else {
      #   did_boot_list <- NULL
      # }

      ## save obj
      did_save <- list("ATT" = did_est, 'results_bootstraps' = did_boot_list)
    } else {
      # RETURN TWO-WAY FIXED EFFECT ESTIMATE
      did_est <- fit[[1]]$coef[1]
      if (isTRUE(is_covariates)) {
        did_var <- vcovHC(fit[[1]], cluster = 'group', type = 'HC2')[1,1]
      } else {
        did_var <- vcovHC(fit[[1]], cluster = 'group', type = 'HC2')
      }

      tmp_se95 <- c(did_est + qnorm(0.025) * sqrt(did_var), did_est + qnorm(1 - 0.025) * sqrt(did_var))
      tmp_se90 <- c(did_est + qnorm(0.050) * sqrt(did_var), did_est + qnorm(1 - 0.050) * sqrt(did_var))

      did_boot_list <- list('boot_est' = NULL, 'ci95' = tmp_se95, 'ci90' = tmp_se90)
      did_save <- list("ATT" = did_est, 'results_bootstraps' = did_boot_list)

    }

    # ********************************************************* #
    #                                                           #
    #                     save results                          #
    #                                                           #
    # ********************************************************* #

    result[[tt]] <- list(
      'results_estimates' = tmp,
      'results_bootstraps' = tmp_min,
      'results_standardDiD' = did_save,
      'BIC' = BIC, "HQIC" = HQIC,
      'BIC_min' = min(BIC), 'HQIC_min' = min(HQIC),
      'min_model' = min_model,
      'select' = select,
      'ATT' = tmp[[min_model]]$ATT,
      'ci95' = ci95,
      'ci90' = ci90
    )

    attr(result[[tt]], 'post_treat') <- attr(data[[tt]], 'post_treat')
    attr(result[[tt]], 'boot')       <- TRUE
    attr(result[[tt]], 'method')     <- 'parametric'
  }

  class(result) <- "diddesign"
  return(result)
}


##
##
##  helper functions
##
##

#' Get demeaned X and y
#' @importFrom plm pmodel.response
#' @importFrom utils getFromNamespace
#' @keywords internal
getX <- function(fm, fit, dat) {

  # import function
  model.matrix <- getFromNamespace("model.matrix.pFormula", "plm")

  ## get y_bar
  # y_bar  <- rep(NA, nrow(dat))
  outvar <- all.vars(fm)[1]
  is_na  <- is.na(dat[,outvar])
  # y_bar[!is_na] <- pmodel.response(fit)

  ## this returns only non-na entries
  y_bar         <- pmodel.response(fit)
  names(y_bar)  <- outvar

  ## get X_bar
  X_bar <- model.matrix(fm, data = dat, model = 'within', effect = 'twoways')
  # dat[!is_na, colnames(X_bar)] <- X_bar
  return(list("y" = y_bar, "X" = X_bar, 'is_na' = is_na))
}


##
## use full X
##

#' Estimate GMM
#' @keywords internal
didgmmT_parametric <- function(dat, id_subject, par_init = NULL) {

  ## prepara inputs
  k <- length(dat)            ## number of y types = number of models
  p <- ncol(dat[[1]]$X)       ## number of covariates

  if(is.null(par_init)) par_init <- runif(1 + k * (p - 1))
  est <- optim(par = par_init, fn = cugmm_loss_parametric, method = "BFGS",
               dat = dat, id_subject = id_subject, k = k, p = p)

  return(list('est' = est, 'ATT' = est$par))
}


#' gmm loss function
#' @param par parameter vector.
#' @param dat list of data. List of outputs from \code{getX} function.
#' @param id_subject vector of subject id.
#' @param k number of outcome types (models)
#' @param p number of variables
#' @keywords internal
cugmm_loss_parametric <- function(par, dat, id_subject, k, p) {

  ## format parameter: first element is beta
  par_mat      <- matrix(NA, nrow = p, ncol = k)

  par_mat[1,]  <- par[1]
  par_mat[-1,] <- par[-1]
  par_mat      <- t(par_mat)

  ## comput E[X'(y - Xb)] = 0 condition
  uid  <- unique(id_subject)
  nobs <- length(uid)

  XXe <- list()
  for (j in 1:k) {
    is_na <- dat[[j]]$is_na
    XXe[[j]] <- loss_loop(X = dat[[j]]$X, y = dat[[j]]$y, par = par_mat[j,],
                          id_subject = id_subject, uid = uid,
                          is_na = as.numeric(is_na), nobs = nobs, p = p )
  }

  MXe   <- do.call("cbind", XXe)
  gbar  <- colMeans(MXe)
  Omega <- (t(MXe) %*% MXe)
  diag(Omega) <- diag(Omega)
  loss  <- as.vector(t(gbar) %*% solve(Omega / nobs, gbar))
  # cat("loss = ", loss, '\n')
  return(loss)
}

#' compute asymptotic GMM variance
#' @keywords internal
cugmm_var_parametric <- function(par, dat, id_subject) {
  ## prepara inputs
  k <- length(dat)            ## number of y types = number of models
  p <- ncol(dat[[1]]$X)       ## number of covariates

  ## format parameter: first element is beta
  par_mat      <- matrix(NA, nrow = k, ncol = p)
  par_mat[,1]  <- par[1]
  par_mat[,-1] <- par[-1]

  ## comput E[X'(y - Xb)] = 0 condition
  uid  <- unique(id_subject)
  nobs <- length(uid)

  XXe <- list()
  G   <- matrix(0, nrow = k, ncol = 1)
  for (j in 1:k) {
    is_na <- dat[[j]]$is_na
    XXe[[j]] <- loss_loop(X = dat[[j]]$X, y = dat[[j]]$y, par = par_mat[j,],
                          id_subject = id_subject, uid = uid,
                          is_na = as.numeric(is_na), nobs = nobs, p = p )
    G[j,1]   <- sum(as.vector(dat[[j]]$X[,1])^2) / nobs
    # cat("G[j,i] = ", G[j,1], "\n")

  }

  MXe   <- do.call("cbind", XXe)
  gbar  <- colMeans(MXe)
  Omega <- crossprod(MXe, MXe) / nobs

  # G     <- matrix(0, nrow = k * p, ncol = 1)
  # G[seq(1, k * p, p),1] <- 1

  gmm_var <- as.vector(solve(t(G) %*%  solve(Omega, G)))
  att_var <- gmm_var / nobs
  return(att_var)
}


##
## Residualize X
##

#' Estimate GMM with residualized outcome / treatment
#' @keywords internal
didgmmT_parametric_resid <- function(dat, id_subject, par_init = NULL) {

  ## prepara inputs
  k <- length(dat); p <- 1

  ## 1. residualize outcome and treatment
  ##   - Regress Y on X
  ##   - Regress D on X
  for (d in 1:length(dat)) {
    outcome    <- dat[[d]]$y
    treatment  <- dat[[d]]$X[, 1]
    covariates <- dat[[d]]$X[, -1]

    dat[[d]]$y_resid <- lm(outcome ~ covariates - 1)$residuals
    dat[[d]]$d_resid <- matrix(lm(treatment ~ covariates - 1)$residuals,
                              nrow = length(outcome), ncol = 1)
  }

  ## 2. use resideualized outcome to estimate AT
  ##   - moment conditions: E[D(Y- Db)] = 0
  if(is.null(par_init)) par_init <- runif(1)
  est <- optim(par = par_init, fn = cugmm_loss_parametric_resid,
                method = "BFGS", dat = dat, id_subject = id_subject, k = k, p = p)

  return(list('est' = est, 'ATT' = est$par, 'data' = dat))
}



#' gmm loss function
#' @param par parameter vector.
#' @param dat list of data. List of outputs from \code{getX} function.
#' @param id_subject vector of subject id.
#' @param k number of outcome types (models)
#' @param p number of variables
#' @keywords internal
cugmm_loss_parametric_resid <- function(par, dat, id_subject, k, p = 1) {

  ## format parameter: first element is beta
  par_mat      <- matrix(NA, nrow = k, ncol = 1)
  par_mat[,1]  <- par

  ## comput E[X'(y - Xb)] = 0 condition
  uid  <- unique(id_subject)
  nobs <- length(uid)

  XXe <- list()
  for (j in 1:k) {
    is_na <- dat[[j]]$is_na
    XXe[[j]] <- loss_loop(X = dat[[j]]$d_resid, y = dat[[j]]$y_resid, par = par_mat[j,],
                          id_subject = id_subject, uid = uid,
                          is_na = as.numeric(is_na), nobs = nobs, p = p )
  }

  MXe   <- do.call("cbind", XXe)
  gbar  <- colMeans(MXe)
  Omega <- (t(MXe) %*% MXe)
  diag(Omega) <- diag(Omega)
  loss  <- as.vector(t(gbar) %*% solve(Omega / nobs, gbar))
  # cat("loss = ", loss, '\n')
  return(loss)
}

#' compute asymptotic GMM variance
#' @keywords internal
cugmm_var_parametric_resid <- function(par, dat, id_subject, p = 1) {
  ## prepara inputs
  k <- length(dat)            ## number of y types = number of models


  ## format parameter: first element is beta
  par_mat      <- matrix(NA, nrow = k, ncol = 1)
  par_mat[,1]  <- par

  ## comput E[X'(y - Xb)] = 0 condition
  uid  <- unique(id_subject)
  nobs <- length(uid)

  XXe <- list()
  G   <- matrix(NA, nrow = k, ncol = 1)
  for (j in 1:k) {
    is_na <- dat[[j]]$is_na
    XXe[[j]] <- loss_loop(X = dat[[j]]$d_resid, y = dat[[j]]$y_resid, par = par_mat[j,],
                          id_subject = id_subject, uid = uid,
                          is_na = as.numeric(is_na), nobs = nobs, p = p )
    G[j,1]   <- sum(dat[[j]]$d_resid^2) / nobs
  }

  MXe   <- do.call("cbind", XXe)
  gbar  <- colMeans(MXe)
  Omega <- crossprod(MXe, MXe) / nobs

  ## dgi / d beta = [D^2]
  # G     <- matrix(0, nrow = k * p, ncol = 1)
  # G[seq(1, k * p, p),1] <- 1

  gmm_var <- as.vector(solve(t(G) %*%  solve(Omega, G)))
  att_var <- gmm_var / nobs
  return(att_var)
}




plm.boot <- function(fm, data, n_boot) {
  ##
  ## block bootstrap
  ##

  id_unique <- unique(data$id_subject)
  n_subject <- length(id_unique)
  fit <- list()
  for (bb in 1:n_boot) {
    ## sample id
    id_use <- sample(id_unique, size = n_subject, replace = TRUE)

    ## resample data
    dat_use <- data[data$id_subject %in% id_use, ]

    fit[[bb]] <- plm(fm, data = dat_use, model = 'within', effect = 'twoways')$coef
  }

  return(do.call('rbind', fit))

}
