


## standard difference-in-differeces estimator 
std_did <- function(Y, D) {
  ## ==== input check ==== ## 
  t_total <- ncol(Y)
  if (t_total == 2) {
    Y_use 
  } else if (t_total > 2) {
    Y_use <- Y[, c(t_total-1, t_total)]
  } else {
    stop("Data should have T >= 2")
  }
  
  ## DiD 
  est1 <- tapply(Y_use[,2], D, mean, na.rm = TRUE)
  est0 <- tapply(Y_use[,1], D, mean, na.rm = TRUE)
  return(diff(est1) - diff(est0))
  
}

#' Block Bootstrap Function for Standard DiD Estimator
#'
std_did_boot <- function(Y, D, n_boot) {
  N <- nrow(Y)
  boot_att <- rep(NA, n_boot)
  for (b in 1:n_boot) {
    use_id <- sample(1:N, size = N, replace = TRUE)
    Yboot <- Y[use_id,]; Dboot <- D[use_id]
    boot_att[b] <- std_did(Y = Yboot, D = Dboot)
  }  
  
  return(boot_att)
}
