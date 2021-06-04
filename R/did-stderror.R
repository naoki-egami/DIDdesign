

##
## Functions to compute standard errors
##


ddid_resid <- function(formula, data, lead = 0) {
  time_use <- c(-1, lead)
  est <- map(formula, function(fm) {
    fit <- lm(fm, data = filter(data, .data$id_time_std %in% time_use))$resid
    return(fit)
  })

  ## create a matrix of residuals
  resid <- do.call(cbind, est)

  return(resid)
}


#' Obtain variance covariance matrix
get_vcov <- function(resid, data, var_cluster) {
  if (is.null(var_cluster)) {
    ## variance covariance matrix of residuals
    vcov <- cov(resid)
  } else {
    # id_cluster_vec <- pull(dat_did, !!sym(var_cluster)))
    # id_unique <- unique(id_cluster_vec)
    # for (i in 1:length(id_unique)) {
    #
    # }

  }


}
