

##
## ouputs
##


#' Summarize DIDdesign output
#' @param object An object of \code{DIDdesign} class.
#' @export
summary.DIDdesign <- function(object, ...) {

  if (!("DIDdesign" %in% class(object))) stop("object should the output of did funciton.")
  out <- object$estimates
  out$statistic <- out$estimate / out$std.error
  out$p_value   <- 2 * pnorm(abs(out$statistic), lower.tail = FALSE)

  out <- out[, c("estimator", "lead", "estimate", "std.error", "statistic", "p_value")]

  # ## add weights
  # weights <- c(NA, object$weights$weight_did, object$weights$weight_sdid)
  # out$ddid_weights <- weights

  class(out) <- c("summary.DIDdesign", class(out))
  return(out)
}


#' Print
#' @export
#' @importFrom cli cat_rule
#' @param x An object of \code{summary.DIDdesign} class. This is typically an output of \code{summary.DIDdesign()} function.
print.summary.DIDdesign <- function(x, ...) {
  cat_rule(left = crayon::bold("ATT Estimates"))
  print(as.data.frame(x), digits = 2)
  # x_out <- data.frame(x)
  # rownames(x_out) <- 1:nrow(x_out)
  # print.default(x_out, digits = 3)
  invisible(x)
  x
}



#' Summary did_check Output
#' @export
summary.DIDdesign_check <- function(object, ...) {
  if (!("DIDdesign_check" %in% class(object))) stop("object should the output of did_check funciton.")
  tmp <- object$estimate
  tmp <- tmp[,c("estimate", "lag", "std.error", "EqCI95_LB", "EqCI95_UB")]
  class(tmp) <- c("summary.DIDdesign_check", class(tmp))
  return(tmp)
}

#' Print
#' @export
#' @importFrom cli cat_rule
print.summary.DIDdesign_check <- function(x, ...) {
  cat_rule(left = crayon::bold("Standardized Estimates"))
  x_out <- data.matrix(x)
  rownames(x_out) <- 1:nrow(x_out)
  print.default(x_out, quote = FALSE, right = TRUE, digits = 3)
  invisible(x)
}
