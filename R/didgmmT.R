##
## general version
##

#' DID-GMM with mth order parallel treand assumption
#'
#' @param Y Outcome matrix.
#' @param D Treatment vector.
#' @param M Order specification.
#' @param ep small constant added to diagonal elements of the weighting matrix in GMM. Default is 0.01.
#' @param only_beta A boolean argument.
#'  If \code{TRUE} ATT is estimated only using moment conditions related to the last time-period.
#'  Default is \code{FALSE}.
#' @export
didgmmT <- function(Y, D, M, ep = 0.01, max_trial = 100, only_beta = FALSE, only_oneM = FALSE) {

  ## quantities
  T0 <- ncol(Y) - 1; N <- nrow(Y)

  ## input checks
  if (M > T0) stop("M should be smaller than or equal to T0.")

  ## count the number of  over-identification conditions
  if (isTRUE(only_beta)) {
    n_param <- 1
    n_beta  <- T0 - M + 1
    n_over  <- n_beta - n_param
  } else {
    n_param <- 2
    n_beta <- T0 - M + 1    # psi^k_{T0+1},..., psi^T_{T0+1}
    n_zeta <- T0 - M + 1    # varphi^{k-1}_{k}, ..., varphi^{k-1}_{T0}
    n_over <- n_beta + n_zeta - n_param
  }

  n_over <- max(0, n_over)

  if (isTRUE(only_beta) & isTRUE(only_oneM)) {
    ## this is simply Double DiD
    method_use  <- 1
  } else if (!isTRUE(only_beta) & isTRUE(only_oneM)) {
    ## use zeta
    method_use <- 2
  } else if (isTRUE(only_beta) & !isTRUE(only_oneM)) {
    method_use <- 3
  } else {
    method_use <- 4
  }


  ## estimate param
  converge <- 1; iter <- 1
  while(converge == 1 & iter <= max_trial) {
    if (n_over == 0) {
      ## random initialization
      fit <- optim(par = rnorm(1), fn = hansenT_just, method = "BFGS", Y = Y, D = D, M = M)

      ## convergence check
      converge <- ifelse(fit$convergence == 0, 0, 1)
      iter <- iter + 1
    } else if (n_over > 0) {
      if (method_use == 1) {
        ## only double did with one moment
        fit <- optim(par = rnorm(1), fn = hansenT_just, method = "BFGS", Y = Y, D = D, M = M)

        ## convergence check
        converge <- ifelse(fit$convergence == 0, 0, 1)
        iter <- iter + 1

      } else if (method_use == 2) {
        ## use zeta but only top M
        ## pick initial value with Identify W GMM
        if (iter == 1) {
          init <- optim(par = rnorm(2), fn = hansenT_over_init, method = 'BFGS', Y = Y, D = D, M = M)$par
        }

        ## run until convergence with different initi values
        fit <- optim(par = init, fn = hansenT_over2, method = 'BFGS', Y = Y, D = D, M = M, ep = ep)

        ## convergence
        converge <- ifelse(fit$convergence == 0, 0, 1)
        iter <- iter + 1
        init <- rnorm(2)
      } else if (method_use == 3) {
        ## pick initial value with Identify W GMM
        if (iter == 1) {
          init <- optim(par = rnorm(1), fn = hansenT_over_init_onlyB, method = 'BFGS', Y = Y, D = D, M = M)$par
        }

        ## run until convergence with different initi values
        fit <- optim(par = init, fn = hansenT_over_onlyB, method = 'BFGS', Y = Y, D = D, M = M)

        ## convergence
        converge <- ifelse(fit$convergence == 0, 0, 1)
        iter <- iter + 1
        init <- rnorm(1)
      } else {
        ## pick initial value with Identify W GMM
        if (iter == 1) {
          init <- optim(par = rnorm(2), fn = hansenT_over_init, method = 'BFGS', Y = Y, D = D, M = M, ep = ep)$par
        }

        ## run until convergence with different initi values
        fit <- optim(
          par = init, fn = hansenT_over, method = 'BFGS', Y = Y, D = D, M = M, ep = ep
        )

        ## convergence
        converge <- ifelse(fit$convergence == 0, 0, 1)
        iter <- iter + 1
        init <- rnorm(2)
      }

    }
  }

  ## compute relevant statistics
  ATT    <- fit$par[1]
  zeta   <- fit$par[-1]
  Jstats <- fit$value * N
  BIC    <- Jstats - n_over * log(N)
  HQIC   <- Jstats - 2.01 * n_over * log(log(N))

  return(list("ATT" = ATT, "zeta" = zeta, "J" = Jstats, "BIC" = BIC, "HQIC" = HQIC))
}




didgmmT.boot <- function(Y, D, M, n_boot, only_beta = FALSE, only_oneM = FALSE) {
  ## give warning when proportion of treated / control is small
  # if (isTRUE(se.boot)) {
  #   tr_prop <- mean(D, na.rm = TRUE)
  # }
  N <- nrow(Y)
  boot_att <- rep(NA, n_boot)
  for (b in 1:n_boot) {
    use_id <- sample(1:N, size = N, replace = TRUE)
    Yboot <- Y[use_id,]; Dboot <- D[use_id]
    tmp <- didgmmT(Y = Yboot, D = Dboot, M = M, only_beta = only_beta, only_oneM = only_oneM)
    boot_att[b] <- tmp$ATT
  }

  return(boot_att)
}


#' Objective function
hansenT_just <- function(par, Y, D, M) {
  m_use <- M
  tmax <- 1
  N1    <- sum(D); N  <- length(D)
  Dmat  <- matrix(D, nrow = length(D), ncol = tmax)
  ydiff <- matrix(NA, nrow = tmax, ncol = N)

  for (z in 1:length(m_use)) {
    tmp <- diff(t(Y), differences = m_use[z])
    ydiff[z,] <- tmp[nrow(tmp),]
  }

  ## compute the loss
  pi_i <- N1 / N
  tmp  <- t(ydiff) / pi_i * (D - pi_i) / (1 - pi_i) - par
  x    <- colSums(tmp) / N
  W    <- crossprod(tmp, tmp)
  # return(t(x) %*% solve(W / N,  x))
  return(t(x) %*% x)
}


hansenT_over <- function(par, Y, D, M, ep) {
  ## params
  beta <- par[1]; zeta <- par[2]

  ## common quantities
  T0 <- ncol(Y) - 1; N <- nrow(Y); N1 <- sum(D)
  pi_i <- N1 / N
  pi_weight <- 1 / pi_i * (D - pi_i) / (1 - pi_i)

  ## compute moments for beta
  m_use <- M:T0
  yTdiff <- matrix(NA, nrow = length(m_use), ncol = N)
  for (z in 1:length(m_use)) {
    ## take m-th order diff
    tmp <- diff(t(Y), differences = m_use[z])
    ## take T0+1 time period
    yTdiff[z,] <- tmp[nrow(tmp),]
  }

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
  psi_out <- cbind(
    t(yTdiff) * pi_weight - beta,
    t(yT0diff) * pi_weight - zeta
  )

  ## compute the loss
  x <- colSums(psi_out) / N
  W <- crossprod(psi_out, psi_out)
  diag(W) <- diag(W) + ep

  ## return gWg
  loss <- t(x) %*% solve(W / N,  x)
  return(loss)
}



hansenT_over2 <- function(par, Y, D, M, ep) {
  ## params
  beta <- par[1]; zeta <- par[2]

  ## common quantities
  T0 <- ncol(Y) - 1; N <- nrow(Y); N1 <- sum(D)
  pi_i <- N1 / N
  pi_weight <- 1 / pi_i * (D - pi_i) / (1 - pi_i)

  ## compute moments for beta
  # m_use <- M:T0
  m_use <- M ## only use the highest order
  yTdiff <- matrix(NA, nrow = length(m_use), ncol = N)
  for (z in 1:length(m_use)) {
    ## take m-th order diff
    tmp <- diff(t(Y), differences = m_use[z])
    ## take T0+1 time period
    yTdiff[z,] <- tmp[nrow(tmp),]
  }

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
  psi_out <- cbind(
    t(yTdiff) * pi_weight - beta,
    t(yT0diff) * pi_weight - zeta
  )

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
hansenT_over_init <- function(par, Y, D, M, ep) {
  ## params
  beta <- par[1]; zeta <- par[2]

  ## common quantities
  T0 <- ncol(Y) - 1; N <- nrow(Y); N1 <- sum(D)
  pi_i <- N1 / N
  pi_weight <- 1 / pi_i * (D - pi_i) / (1 - pi_i)

  ## compute moments for beta
  m_use <- M:T0
  yTdiff <- matrix(NA, nrow = length(m_use), ncol = N)
  for (z in 1:length(m_use)) {
    ## take m-th order diff
    tmp <- diff(t(Y), differences = m_use[z])
    ## take T0+1 time period
    yTdiff[z,] <- tmp[nrow(tmp),]
  }

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
  psi_out <- cbind(
    t(yTdiff) * pi_weight - beta,
    t(yT0diff) * pi_weight - zeta
  )

  ## compute the loss
  x <- colSums(psi_out) / N
  # W <- crossprod(psi_out, psi_out)
  # diag(W) <- diag(W) + ep

  ## return gWg
  loss <- crossprod(x, x)
  return(loss)
}





hansenT_over_onlyB <- function(par, Y, D, M) {
  ## params
  beta <- par[1]

  ## common quantities
  T0 <- ncol(Y) - 1; N <- nrow(Y); N1 <- sum(D)
  pi_i <- N1 / N
  pi_weight <- 1 / pi_i * (D - pi_i) / (1 - pi_i)

  ## compute moments for beta
  m_use <- M:T0
  yTdiff <- matrix(NA, nrow = length(m_use), ncol = N)
  for (z in 1:length(m_use)) {
    ## take m-th order diff
    tmp <- diff(t(Y), differences = m_use[z])
    ## take T0+1 time period
    yTdiff[z,] <- tmp[nrow(tmp),]
  }

  ## combine moments
  psi_out <- t(yTdiff) * pi_weight - beta

  ## compute the loss
  x <- colSums(psi_out) / N
  W <- crossprod(psi_out, psi_out)

  ## return gWg
  loss <- t(x) %*% solve(W / N,  x)

  return(loss)
}


##
## choose initial value
##
hansenT_over_init_onlyB <- function(par, Y, D, M) {
  ## params
  beta <- par[1]

  ## common quantities
  T0 <- ncol(Y) - 1; N <- nrow(Y); N1 <- sum(D)
  pi_i <- N1 / N
  pi_weight <- 1 / pi_i * (D - pi_i) / (1 - pi_i)

  ## compute moments for beta
  m_use <- M:T0
  yTdiff <- matrix(NA, nrow = length(m_use), ncol = N)
  for (z in 1:length(m_use)) {
    ## take m-th order diff
    tmp <- diff(t(Y), differences = m_use[z])
    ## take T0+1 time period
    yTdiff[z,] <- tmp[nrow(tmp),]
  }

  psi_out <- t(yTdiff) * pi_weight - beta

  ## compute the loss
  x <- colSums(psi_out) / N

  ## return gWg
  loss <- crossprod(x, x)
  return(loss)
}
