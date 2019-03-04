


#' Pre-treatment Trend test with GMM
#' @keywords internal
didgmm_test2 <- function(Y0, D0, M, max_trial = 100) {

  T0 <- ncol(Y0); N  <- nrow(Y0)

  if (M >= T0) stop("M should be smaller than T0 for specification test.")
  n_over <- (T0-1) - M
  converge <- 1; iter <- 1

  if (n_over == 0) {
    ## this is simple  "T0-1"-th order difference-in-differences
    YD1 <- apply(Y0[D0 == 1,], 2, mean, na.rm = TRUE);
    YD0 <- apply(Y0[D0 == 0,], 2, mean, na.rm = TRUE);
    fit <- list("par" = diff(YD1, difference = M) - diff(YD0, difference = M), "value" = 0)

    ## variance (simply adding variances )
    varn1 <- apply(Y0[D0 == 1,], 2, function(x) var(x));
    varn0 <- apply(Y0[D0 == 0,], 2, function(x) var(x));
    se_est <- sqrt(sum(varn1) + sum(varn0))
  } else {
    ## pick initial value with Identify W GMM
    init <- optim(par = rnorm(1), fn = hansenT_tretest_over2, method = 'BFGS',
                  Y = Y0, D = D0, M = M, init = TRUE)$par

    ## run until convergence with different initi values
    while(converge == 1 & iter <= max_trial) {
      fit <- optim(par = init, fn = hansenT_tretest_over2, method = 'BFGS',
                   Y = Y0, D = D0, M = M)

      ## convergence
      converge <- ifelse(fit$convergence == 0, 0, 1)
      iter <- iter + 1
      init <- rnorm(1)
    }

    ## GMM variance
    se_est <- sqrt(hansenT_tretest_over2(fit$par, Y = Y0, D = D0, M = M, return_var = TRUE))

  }

  ## compute relevant statistics
  zeta   <- fit$par
  Jstats <- fit$value * N
  BIC    <- Jstats - n_over * log(N)
  HQIC   <- Jstats - 2.01 * n_over * log(log(N))
  return(list("zeta" = zeta, "J" = Jstats, "BIC" = BIC, "HQIC" = HQIC, 'se' = se_est))
}




#' compute loss based on moment condition
#' @keywords internal
hansenT_tretest_over2 <- function(par, Y0, D0, M, init = FALSE, return_var = FALSE) {
  ## params
  zeta <- par

  ## common quantities
  T0 <- ncol(Y0) - 1; N <- nrow(Y0); N1 <- sum(D0)
  pi_i <- N1 / N
  pi_weight <- 1 / pi_i * (D0 - pi_i) / (1 - pi_i)

  ## compute moments for zeta
  t_use <- (M+1):T0
  yT0diff <- matrix(NA, nrow = length(t_use), ncol = N)
  for (z in 1:length(t_use)) {
    ## take (k-1)-th order diff
    tmp <- diff(t(Y0[,1:t_use[z]]), differences = M)
    ## take 'last' time period among subseted
    yT0diff[z,] <- tmp[nrow(tmp),]
  }

  ## combine moments
  psi_out <- t(yT0diff) * pi_weight - zeta

  ## compute the loss
  x <- colSums(psi_out) / N
  if(isTRUE(init)) {
    W <- diag(ncol(psi_out))
  } else {
    W <- crossprod(psi_out, psi_out) / N
  }


  ## return gWg
  if (isTRUE(return_var)) {
    gmm_var <- 1 / sum(solve(W))
    loss    <- gmm_var / N
  } else {
    loss <- t(x) %*% solve(W,  x)
  }

  return(loss)
}



#' Trend test with mean
#' @keywords internal
trendT_test <- function(Y, D, M, boot = FALSE, n_boot = 1000, alpha = 0.05) {
  T0 <- ncol(Y) - 1; N <- nrow(Y)

  ## input checks
  if (M >= T0) stop("M should be smaller than to T0.")

  ## compute quantities
  YD1 <- Y[D == 1, 1:T0]
  YD0 <- Y[D == 0, 1:T0]

  ## group mean
  YD1M <- colMeans(YD1)
  YD0M <- colMeans(YD0)
  deltaT <- YD1M - YD0M

  ## take diff
  delta <- diff(deltaT, difference = M)

  ## variance　by boot
  delta_boot <- list()
  if (isTRUE(boot)) {
    tmp <- list()
    for (bb in 1:n_boot) {
      use_id <- sample(1:N, size = N, replace = TRUE)
      Yboot <- Y[use_id,]; Dboot <- D[use_id]

      ## group mean
      YD1M <- colMeans(Yboot[Dboot == 1, 1:T0])
      YD0M <- colMeans(Yboot[Dboot == 0, 1:T0])
      deltaT <- YD1M - YD0M

      ## take diff
      tmp[[bb]] <- diff(deltaT, difference = M)
    }
    delta_boot_res <- do.call("rbind", tmp)
    delta_zscore  <- delta / apply(delta_boot_res, 2, sd)
    delta_pval <- sapply(1:length(delta), function(i) {
      2 * min(pnorm(delta_zscore[i]), 1 - pnorm(delta_zscore[i]))
    })

    reject <- ifelse(min(delta_pval) <= alpha, TRUE, FALSE)
    reject_BF <- ifelse(min(delta_pval) <= alpha / length(delta_pval), TRUE, FALSE)
    # delta_pval <- apply(cbind(pnorm(delta_zscore), 1 - pnorm(delta_zscore)), 1,
    # function(x) 2 * min(x))

    delta_boot <- list(
      "estimates" = delta_boot_res, 'z-score' = delta_zscore, 'pval' = delta_pval,
      'reject' = reject, 'reject_correction' = reject_BF
    )

  }

  return(list("delta" = delta, 'boot' = delta_boot))

}




## *********************************************** ##
##                                                 ##
##  LEGACY CODE BELOW: JUST FOR ARCHIVAL PURPOSES  ##
##                                                 ##
## *********************************************** ##

#' Pre-treatment Trend test with GMM
#' @keywords internal
didgmm_test <- function(Y, D, M, ep = 0.01, max_trial = 100) {

  ## quantities
  T0 <- ncol(Y) - 1; N <- nrow(Y)

  ## input checks
  if (M > T0) stop("M should be smaller than or equal to T0.")

  ## count the number of  over-identification conditions
  n_param <- 1
  n_zeta <- T0 - M + 1    # varphi^{k-1}_{k}, ..., varphi^{k-1}_{T0}
  n_over <- n_zeta - n_param
  n_over <- max(0, n_over)

  ## estimate param
  converge <- 1; iter <- 1
  if (n_over == 0) {
    ## this is simple  "T0-1"-th order difference-in-differences
    YD1 <- apply(Y[D == 1, 1:T0], 2, mean, na.rm = TRUE);
    YD0 <- apply(Y[D == 0, 1:T0], 2, mean, na.rm = TRUE);
    fit <- list("par" = diff(YD1, difference = (T0-1)) - diff(YD0, difference = (T0-1)), "value" = 0)

    ## variance (simply adding variances )
    varn1 <- apply(Y[D == 1, 1:T0], 2, function(x) var(x));
    varn0 <- apply(Y[D == 0, 1:T0], 2, function(x) var(x));
    se_est <- sqrt(sum(varn1) + sum(varn0))
  } else {
    ## pick initial value with Identify W GMM
    init <- optim(
      par = rnorm(1), fn = hansenT_tretest_over_init, method = 'BFGS', Y = Y, D = D, M = M, ep = ep
    )$par

    ## run until convergence with different initi values
    while(converge == 1 & iter <= max_trial) {
      fit <- optim(
        par = init, fn = hansenT_tretest_over, method = 'BFGS', Y = Y, D = D, M = M, ep = ep
      )

      ## convergence
      converge <- ifelse(fit$convergence == 0, 0, 1)
      iter <- iter + 1
      init <- rnorm(2)
    }

    ## GMM variance
    se_est <- sqrt(hansenT_tretest_variance(fit$par, Y = Y, D = D, M = M))
  }

  ## compute relevant statistics
  zeta   <- fit$par
  Jstats <- fit$value * N
  BIC    <- Jstats - n_over * log(N)
  HQIC   <- Jstats - 2.01 * n_over * log(log(N))
  return(list("zeta" = zeta, "J" = Jstats, "BIC" = BIC, "HQIC" = HQIC, 'se' = se_est))

}

#' compute loss based on moment condition
#' @keywords internal
hansenT_tretest_over <- function(par, Y, D, M, ep) {
  ## params
  zeta <- par

  ## common quantities
  T0 <- ncol(Y) - 1; N <- nrow(Y); N1 <- sum(D)
  pi_i <- N1 / N
  pi_weight <- 1 / pi_i * (D - pi_i) / (1 - pi_i)

  ## compute moments for zeta
  t_use <- M:T0
  yT0diff <- matrix(NA, nrow = length(t_use), ncol = N)
  for (z in 1:length(t_use)) {
    ## take (k-1)-th order diff
    if (M == 1) {
      tmp <- t(Y[,1:t_use[z]])
    } else {
      tmp <- diff(t(Y[,1:t_use[z]]), differences = M-1) ## fix this
    }
    ## take t time period
    yT0diff[z,] <- tmp[nrow(tmp),]
  }

  ## combine moments
  pi_i <- N1 / N
  pi_weight <- 1 / pi_i * (D - pi_i) / (1 - pi_i)
  psi_out <- t(yT0diff) * pi_weight - zeta

  ## compute the loss
  x <- colSums(psi_out) / N
  W <- crossprod(psi_out, psi_out)
  diag(W) <- diag(W) + ep

  ## return gWg
  loss <- t(x) %*% solve(W / N,  x)
  return(loss)
}


#' function to select intial value
#' @keywords internal
hansenT_tretest_over_init <- function(par, Y, D, M, ep) {
  ## params
  zeta <- par

  ## common quantities
  T0 <- ncol(Y) - 1; N <- nrow(Y); N1 <- sum(D)
  pi_i <- N1 / N
  pi_weight <- 1 / pi_i * (D - pi_i) / (1 - pi_i)

  ## compute moments for zeta
  t_use <- M:T0
  yT0diff <- matrix(NA, nrow = length(t_use), ncol = N)
  for (z in 1:length(t_use)) {
    ## take k-th order diff
    if (M == 1) {
      tmp <- t(Y[,1:t_use[z]])
    } else {
      tmp <- diff(t(Y[,1:t_use[z]]), differences = M-1)
    }

    ## take t time period
    yT0diff[z,] <- tmp[nrow(tmp),]
  }

  ## combine moments
  pi_i <- N1 / N
  pi_weight <- 1 / pi_i * (D - pi_i) / (1 - pi_i)
  psi_out <- t(yT0diff) * pi_weight - zeta

  ## compute the loss
  x <- colSums(psi_out) / N
  # W <- crossprod(psi_out, psi_out)
  # diag(W) <- diag(W) + ep

  ## return gWg
  loss <- crossprod(x, x)
  return(loss)
}

#' GMM variance for the testing moment model
#' @keywords
hansenT_tretest_variance <- function(par, Y, D, M) {
  ## params
  zeta <- par

  ## common quantities
  T0 <- ncol(Y) - 1; N <- nrow(Y); N1 <- sum(D)
  pi_i <- N1 / N
  pi_weight <- 1 / pi_i * (D - pi_i) / (1 - pi_i)

  ## compute moments for zeta
  t_use <- M:T0
  yT0diff <- matrix(NA, nrow = length(t_use), ncol = N)
  for (z in 1:length(t_use)) {
    ## take (k-1)-th order diff
    if (M == 1) {
      tmp <- t(Y[,1:t_use[z]])
    } else {
      tmp <- diff(t(Y[,1:t_use[z]]), differences = M-1) ## fix this
    }
    ## take t time period
    yT0diff[z,] <- tmp[nrow(tmp),]
  }

  ## combine moments
  pi_i <- N1 / N
  pi_weight <- 1 / pi_i * (D - pi_i) / (1 - pi_i)
  psi_out <- t(yT0diff) * pi_weight - zeta

  ## compute the variance
  W <- crossprod(psi_out, psi_out) / N
  ones <- rep(1, nrow(W))
  gmm_var <- 1 / as.vector(ones %*% solve(W, ones))
  par_var <- gmm_var / N

  return(par_var)
}
