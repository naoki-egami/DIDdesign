#' Difference-in-Differences Design
#'
#' This packages provides tools for estimating treatment effects
#'  and diagnosing difference-in-differences designs.
#'
#' The package has the following main functions:
#' \itemize{
#'   \item \code{\link{did_plot}}: create a DiD plot to visualize difference-in-differences design.
#'   \item \code{\link{did}}: estimate treatment effects.
#'   \item \code{\link{plot.diddesign}}: visualize estimates.
#'   \item \code{\link{summary.diddesign}}: summarize estimates.
#'   \item \code{\link{did_plot_selection}}: visualize statistics used for model selection.
#' }
"_PACKAGE"


#' @useDynLib DIDdesign
#' @importFrom Rcpp sourceCpp
NULL
