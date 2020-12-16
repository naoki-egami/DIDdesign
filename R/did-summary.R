

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

  class(out) <- c(class(out), "summary.DIDdesign.")
  return(out)
}


#' Print
#' @export
#' @importFrom cli cat_rule
#' @param x An object of \code{summary.DIDdesign} class. This is typically an output of \code{summary.DIDdesign()} function.
print.summary.DIDdesign <- function(x, ...) {
  cat_rule(left = crayon::bold("Estimates"))
  print(x)
  invisible(x)
}
