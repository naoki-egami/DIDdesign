

#' Pre-treatment Trend test with GMM
#' @export
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
    ## this is simple
    YD1 <- apply(Y[D == 1, 1:T0], 2, mean, na.rm = TRUE);
    YD0 <- apply(Y[D == 0, 1:T0], 2, mean, na.rm = TRUE);
    fit <- list("par" = diff(YD1, difference = (T0-1)) - diff(YD0, difference = (T0-1)), "value" = 0)
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
  }

  ## compute relevant statistics
  zeta   <- fit$par
  Jstats <- fit$value * N
  BIC    <- Jstats - n_over * log(N)
  HQIC   <- Jstats - 2.01 * n_over * log(log(N))
  return(list("zeta" = zeta, "J" = Jstats, "BIC" = BIC, "HQIC" = HQIC))

}



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


##
## choose initial value
##
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



#' Trend test with mean
#' @export
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
