#' Dataset from Anzia (2012)
#'
#'
#' @name anzia2012
#' @docType data
#' @format A \code{tbl} format data.
#' @references Anzia, Sarah F.
#'  "The Election Timing Effect: Evidence from a Policy Intervention in Texas."
#'  \emph{Quarterly Journal of Political Science} 7.3 (2012):209-248.
#'  \doi{10.1561/100.00011056}
#' @keywords dataset
"anzia2012"


#' Dataset from de Kadt and Larreguy (2018)
#'
#' @name dekadt2018
#' @docType data
#' @format A \code{tbl} format data.
#' @references  De Kadt, Daniel, and Horacio A. Larreguy.
#'  "Agents of the regime? Traditional leaders and electoral politics in South Africa."
#'  \emph{The Journal of Politics} 80.2 (2018): 382-399.
#'  \doi{10.1086/694540}
#' @examples
#' # load data
#' data(dekadt2018)
#'
#' # did plot
#' did_plot(anc_vs_na ~ treatment, data = dekadt2018,
#'   id_subject = 'ward_id', id_time = 'year',
#'   post_treatment = c(2009, 2011, 2014)
#' )
#'
#' # fit
#' fit <- did(anc_vs_na ~ treatment, data = dekadt2018,
#'   id_subject = 'ward_id', id_time = 'year',
#'   post_treatment = c(2009, 2011, 2014)
#' )
#'
#' # view results
#' summary(fit)
#'
#' par(mar = c(4, 2.5, 3.5, 1))
#' plot(fit, full = TRUE)
#'
#' did_plot_selection(fit, equivalence = FALSE, ylim = c(0, 0.1))
#' @keywords dataset
"dekadt2018"
