#' Summarize DIDdesign output
#' @param object An object of \code{DIDdesign} class.
#' @param estimator A vector of estimators to print. For example, \code{estimator = c("DID", "sDID")}.
#' @param ... Other parameters. Currently not supported.
#' @importFrom stats pnorm
#' @export
summary.DIDdesign <- function(object, estimator = NULL, ...) {

  if (!("DIDdesign" %in% class(object))) stop("object should the output of did funciton.")
  out <- object$estimates
  out$statistic <- out$estimate / out$std.error
  out$p_value   <- 2 * pnorm(abs(out$statistic), lower.tail = FALSE)

  out <- out[, c("estimator", "lead", "estimate", "std.error", "statistic", "p_value")]

  if (!is.null(estimator)) {
    ## check
    name_check <- estimator %in% c("Double-DID", "DID", "sDID", "SA-Double-DID", "SA-DID", "SA-sDID")
    if (!all(name_check)) stop("Not supported estimator(s): ", estimator[!name_check])
    out <- out[out$estimator %in% estimator, ]
  }
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
#' @param ... Other parameters. Currently not supported.
print.summary.DIDdesign <- function(x, ...) {
  cat_rule(left = "ATT Estimates")
  print(as.data.frame(x), digits = 2)
  # x_out <- data.frame(x)
  # rownames(x_out) <- 1:nrow(x_out)
  # print.default(x_out, digits = 3)
  invisible(x)
  x
}



#' Summary did_check Output
#' @importFrom dplyr %>% select
#' @param object An output of \code{did_check} function.
#' @param ... Other parameters. Currently not supported.
#' @export
summary.DIDdesign_check <- function(object, ...) {
  if (!("DIDdesign_check" %in% class(object))) stop("object should the output of did_check funciton.")
  tmp <- object$estimate
  tmp <- tmp %>% select(
    estimate = .data$estimate_orig, .data$lag, std.error = .data$std.error_orig,
    .data$EqCI95_LB, .data$EqCI95_UB
  )

  class(tmp) <- c("summary.DIDdesign_check", class(tmp))
  return(tmp)
}

#' Print
#' @export
#' @param x An output of \code{summary.DIDdesign_check} function.
#' @param ... Other parameters. Currently not supported.
#' @importFrom cli cat_rule
print.summary.DIDdesign_check <- function(x, ...) {
  cat_rule(left = "Estimates for assessing parallel trends assumption")
  x_out <- data.matrix(x)
  rownames(x_out) <- 1:nrow(x_out)
  print.default(x_out, quote = FALSE, right = TRUE, digits = 3)
  invisible(x)
}
