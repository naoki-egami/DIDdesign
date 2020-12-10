

##
## ouputs
##


#' Summarize DIDdesign output
#' @export
summary.DIDdesign <- function(obj) {

  if (!("DIDdesign" %in% class(obj))) stop("obj should the output of did funciton.")
  out <- obj$estimates
  out$statistic <- out$estimate / out$std.error
  out$p_value   <- 2 * pnorm(abs(out$statistic), lower.tail = FALSE)

  ## add weights
  weights <- c(NA, obj$weights$weight_did, obj$weights$weight_sdid)
  out$ddid_weights <- weights

  class(out) <- c(class(out), "summary.DIDdesign.")
  return(out)
}


#' Print
#' @export
#' @importFrom cli cat_rule
print.summary.DIDdesign <- function(obj) {
  cat_rule(left = crayon::bold("Estimates"))
  print(obj)
  invisible(obj)
}
