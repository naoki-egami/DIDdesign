#' Dataset from Anzia (2012)
#'
#'
#' @name anzia2012
#' @usage data(anzia2012)
#' @docType data
#' @format A panel data in \code{tbl} and \code{data.frame} format
#'   with 6,377 rows and 11 variables:
#'  \describe{
#'    \item{district}{district id}
#'    \item{year}{time index}
#'    \item{lnavgsalary_cpi}{main outcome variable}
#'    \item{oncycle}{main treatment variable}
#'  }
#' @references Anzia, Sarah F.
#'  "The Election Timing Effect: Evidence from a Policy Intervention in Texas."
#'  \emph{Quarterly Journal of Political Science} 7.3 (2012):209-248.
#'  \doi{10.1561/100.00011056}
#' @keywords dataset
"anzia2012"


#' Dataset from de Kadt and Larreguy (2018)
#'
#' @name dekadt2018
#' @usage data(dekadt2018)
#' @docType data
#' @format A panel data in \code{tbl} and \code{data.frame} format.
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
#'   post_treatment = c(2009, 2011, 2014),
#'   method = 'parametric'
#' )
#'
#' # view results
#' summary(fit)
#'
#' plot(fit, full = TRUE)
#'
#' did_plot_selection(fit, equivalence = FALSE, ylim = c(0, 0.1))
#' @keywords dataset
"dekadt2018"



#' Dataset from Malesky et al. (2014)
#'
#' @name malesky2014
#' @usage data(malesky2014)
#' @docType data
#' @format A \code{tbl} format data.
#' @references Malesky, Edmund J., Cuong Viet Nguyen, and Anh Tran.
#' "The impact of recentralization on public services:
#'   A difference-in-differences analysis of the abolition of elected councils in Vietnam."
#'   \emph{American Political Science Review} 108.1 (2014): 144-168.
#'   \doi{10.1017/S0003055413000580}
#' @keywords dataset
"malesky2014"
