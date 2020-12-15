
#' Assessing Assumptions for Difference-in-differences
#' @inheritParams did
#' @export
#' @importFrom dplyr mutate select %>%
did_check <- function(
  formula, data, id_unit, id_time, design = "did",
  is_panel = TRUE, option = list()
) {

  ## set option
  option <- set_option(option)

  ## -----------------------------------
  ## implement DID estimator
  ## -----------------------------------
  if (design == "did") {
    ## standard design
    fit <- did_check_std(formula, data, id_unit, id_time, is_panel, option)
  } else if (design == "sa"){
    ## staggered adoption design
    if (isFALSE(is_panel)) stop("Only panel data is supported in the SA design.")
    fit <- did_check_sad(formula, data, id_unit, id_time, option)
  }

  ## -----------------------------------
  ## equivalence CI
  ## -----------------------------------
  estimate <- fit$est %>%
  mutate(
    CI90_UB_ab = abs(.data$estimate + qnorm(0.95) * .data$std.error),
    CI90_LB_ab = abs(.data$estimate - qnorm(0.95) * .data$std.error)
  ) %>%
  mutate(
    EqCI95_LB = -pmax(.data$CI90_UB_ab, .data$CI90_LB_ab),
    EqCI95_UB = pmax(.data$CI90_UB_ab, .data$CI90_LB_ab)
  ) %>%
  select(-.data$CI90_UB_ab, -.data$CI90_LB_ab)

  out <- list(estimate = estimate, plot = fit$plot)
  class(out) <- c(class(out), "DIDdesign_check")
  attr(out, 'design') <- design
  return(out)
}



#' Plot
#' @export
#' @import patchwork
#' @importFrom ggplot2 theme
plot.DIDdesign_check <- function(obj) {

  if (attr(obj, "design") == "sa") {

  } else {
    p1 <- obj$plot[[1]]$plot + theme(aspect.ratio=1)
    p2 <- obj$plot[[2]]$plot + theme(aspect.ratio=1)
    pp <- p2 + p1
  }
  return(pp)
}
