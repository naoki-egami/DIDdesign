

#' Double DiD with Repeated Cross-Section
#' @param data A \code{diddesign_data} object.
#' @param x_colnames A vector of covariate variable names.
#'  Parsed from formula in \code{\link{did}} function.
#' @importFrom CVXR Variable Minimize Problem
#' @importFrom genlasso getDtf
#' @importFrom utils getFromNamespace
#' @importMethodsFrom CVXR %*% psolve
#' @export
did_parametric_repeatedCS <- function(data, x_colnames) {

  ## export function
  # solve.CVXR <- getFromNamespace("psolve", "CVXR")

  ## get data info
  time_pre <- ncol(data[[1]]$Y) - 1

  ##
  ## loop is over the post-treatment periods
  ##
  for (tt in 1:length(data)) {
    pdata  <- data[[tt]]$pdata
    Xdummy <- gen_dummies(data[[tt]], add_intercept = TRUE)
    Xcovar <- pdata[,colnames(pdata) %in% x_colnames]
    # Xbind  <- cbind(Xdummy, Xcovar)

    ## prepare for CVX input
    gamma_coef      <- Variable(ncol(Xcovar))
    beta_treat_coef <- Variable(2)
    beta_time_coef  <- Variable(time_pre+1)
    beta_inter_coef <- Variable(2 * (time_pre+1))
    intercept       <- Variable(1)
    # beta_coef <- c(intercept, beta_treat_coef, beta_time_coef,
    #                beta_inter_coef, gamma_coef)
    objective <- Minimize(sum((as.vector(pdata$outcome) -
      Xdummy$ones %*% intercept -
      Xdummy$group %*% beta_treat_coef - Xdummy$time %*% beta_time_coef -
      Xdummy$interaction %*% beta_inter_coef)^2))
    # objective <- Minimize(sum((data[[tt]]$Y - Xdummy$ones %*% intercept)^2))


    ## prepare constraints
    ## 0. useful quantities
    n_times <- length(unique(pdata$id_time))

    ## 1. common constraints: need these for identification
    const_treat  <- sum(beta_treat_coef) == 0
    const_time   <- sum(beta_time_coef) == 0
    const_inter  <- sum(beta_inter_coef) == 0

    moment_vec <- 1:time_pre

    diff_mat <- get_diffmat(time_pre)

    for (m in 1:length(moment_vec)) {

      if (moment_vec[m] == time_pre) {
        ## cvx problem
        problem <- Problem(objective,
          constraints = list(const_treat, const_time, const_inter))
      } else {
        ## get differencing operator mat (see genlasso)
        Dmat <- getDtf(time_pre, moment_vec[m]-1)
        ## impose additional constraints
        const_add <- (Dmat %*% (diff_mat %*% beta_inter_coef)) == 0

        ## cvx problem
        problem <- Problem(objective,
          constraints = list(const_treat, const_time, const_inter, const_add))
      }

      ## solve the problem
      result <- psolve(problem)

    }

  }








}


# betaHat   <- Variable(ncol(XX) + ncol(Xcov))
# objective <- Minimize(sum((Y - Xuse %*% betaHat)^2))
# const1    <- sum(betaHat[2:4]) == 0
# const2    <- sum(betaHat[5:6]) == 0
# const3    <- sum(betaHat[7:12]) == 0
# problem <- Problem(objective, constraints = list(const1, const2))
# result <- solve(problem)
#
#
# m <- matrix(result$getValue(betaHat), ncol = 1)
# rownames(m) <- paste0("$\\beta_{", 1:ncol(Xuse), "}$")


#' Generate dummy variables
#' @param data An element of \code{diddesign_data} object.
#' @return A data matrix. The first column is ones' if \code{add_intercept = TRUE}.
#' @keywords internal
gen_dummies <- function(data, add_intercept = TRUE) {
  pdata <- data$pdata

  ## generate group dummies
  ## we can use id_subject but make sure it is at the treatment group levels
  group_dummy <- model.matrix(~ factor(treatment) - 1, data = pdata)
  ## generate time dummies
  time_dummy  <- model.matrix(~ factor(id_time) - 1, data = pdata)
  ## generate interaction dummies
  interact_dummy <- model.matrix(~ factor(treatment):factor(id_time) - 1, data = pdata)

  if(isTRUE(add_intercept)) {
    ones <- matrix(1, nrow = nrow(pdata), ncol = 1)
    Xmat <- list('ones' = ones, 'group' = data.matrix(group_dummy),
    'time'  = data.matrix(time_dummy), 'interaction' = data.matrix(interact_dummy))
  } else {
    Xmat <- cbind(group_dummy, time_dummy, interact_dummy)
  }

  return(Xmat)
}



get_diffmat <- function(time_pre) {
  time_total <- time_pre + 1
  diff_mat <- matrix(0, nrow = time_pre, ncol = 2 * time_total)
  counter <- 0
  for (tt in 1:time_pre) {
    diff_mat[tt,1+counter] <--1
    diff_mat[tt,time_total+1+counter] <- 1
    counter <- counter + 1
  }
  return(diff_mat)
}
