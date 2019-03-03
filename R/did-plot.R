
#' DiD plot
#' @examples
#' # load package
#' require(DIDdesign)
#'
#' # load data
#' data(anzia2012)
#'
#' # plot
#' did_plot(
#'   outcome = anzia2012$lnavgsalary_cpi,
#'   treatment = anzia2012$oncycle,
#'   post_treatment = c(2007, 2008, 2009),
#'   id_subject = anzia2012$district,
#'   id_time = anzia2012$year,
#'   ylim = c(10.5, 10.8)
#' )
#' @importFrom dplyr %>% pull tbl_df
#' @importFrom Formula as.Formula
#' @importFrom utils getFromNamespace
#' @export
did_plot <- function(formula, data, post_treatment, id_subject = NULL, id_time,
  xlim = NULL, ylim = NULL, col = NULL, ...
) {

  ## import function
  getFormula <- getFromNamespace("formula.Formula", "Formula")

  ## extract variable infor

  ## convert formula
  formula <- as.Formula(formula)
  f1 <- getFormula(formula, rhs = 1)

  ## extract variable names
  outcome    <- all.vars(f1)[1]
  treatment  <- all.vars(f1)[2]

  ## gether data
  if (is.null(id_subject)) {
    warning("treat data as repeated cross-section data")
    data %>% group_by(treatment, id_time) %>%
      summarise(ymean = mean(outcome, na.rm = TRUE)) -> y_summary
    y1mean <- y_summary %>% filter(treatment == 1) %>% pul(ymean)
    y0mean <- y_summary %>% filter(treatment == 0) %>% pul(ymean)
    # if (is.null(xlim))
    # if (is.null(xlim))
    # plot(1, 1, type = 'n', xlim = xlim, ylim = ylim, col = col, ...)

  } else {
    dat_use <- did_data(
      outcome        = data %>% pull(outcome),
      treatment      = data %>% pull(treatment),
      id_subject     = data %>% pull(id_subject),
      id_time        = data %>% pull(id_time),
      post_treatment = post_treatment,
      long           = TRUE,
      Xcov           = NULL
    )

    ## CALL plot.diddesign
    # par(mfrow = c(1, 3))
    plot(dat_use, xlim = xlim, ylim = ylim, col = col, ...)

  }

}


#' Plot function
#' @param data Input object. This should be either diddesign object or diddesign_data object.
#'  If diddesign object is provided, ATT estimates with 95% confidence intervals are plotted.
#'  If did_double.data objec is provided, mean outcomes for the treated and the control are plotted.
#' @param full A boolean argument, if \code{TRUE} ATT estimates for all models are plotted,
#'  otherwise only estimates from the selected model will be shown. Default is \code{FALSE}.
#' @param xlim xlim in plot function. If left \code{NULL}, it's automatically set.
#' @param ylim ylim in plot function. If left \code{NULL}, it's automatically set.
#' @param ... additional arguments supplied to the plot function.
#' @export
plot.diddesign <- function(data, xlim = NULL, ylim = NULL, col = NULL, lwd = NULL, full = FALSE, ...) {
  if('diddesign_data' %in% class(data)) {
    ##
    ## plot for raw data
    ##

    summary_dat <- summarize.diddesign(data)
    ymean <- summary_dat$ymean
    id_time <- summary_dat$id_time

    ## plot
    if (is.null(xlim)) xlim <- c(min(id_time), max(id_time))
    if (is.null(ylim)) ylim <- c(min(ymean), max(ymean))
    if (is.null(col) | length(col) == 1) col <- c(col, '#006284', 'gray60')
    if (is.null(lwd)) lwd <- 1.5

    time_use <- which(attr(data[[1]], 'post_treat') == id_time)

    ## plot
    plot(1, 1, type = 'n', ylim = ylim, xlim = xlim, xlab = "", ylab = "Mean Outcome", ...)
    abline(v = mean(id_time[(time_use-1):time_use]), col = 'red', lty = 3, lwd = 1.3)
    lines(id_time, ymean[1,], col = col[1], pch = 16, type = 'b', lwd = lwd)
    lines(id_time, ymean[2,], col = col[2],  pch = 17, type = 'b', lwd = lwd)
    legend('topleft', legend = c("treated", 'control'), col = col[1:2],
      lty = 1, pch = c(16, 17), bty = 'n')

  } else if ('diddesign' %in% class(data)) {
    ##
    ## plot for estimation result
    ##
    ATT <- sapply(data, function(x) round(x$ATT, 3))
    xlab <- sapply(data, function(x) attr(x, 'post_treat'))
    if (isTRUE(attr(data[[1]], 'boot') & isTRUE(full)) &
        !is.null(data[[1]]$results_standardDiD)) {
      se_list <- ATT_list <- colors <- list()
      for (i in 1:length(data)) {
        se_mat <- matrix(NA, nrow = length(data[[i]]$results_bootstraps) + 1, ncol = 4)
        att_vec <- rep(NA, length(data[[i]]$results_bootstraps) + 1)

        ## save DiD
        se_mat[1,1:2] <- data[[i]]$results_standardDiD$results_bootstraps$ci90
        se_mat[1,3:4] <- data[[i]]$results_standardDiD$results_bootstraps$ci95
        att_vec[1] <- data[[i]]$results_standardDiD$ATT

        colors[[i]] <- rep('gray60', length(data[[i]]$results_bootstraps))
        for (j in 1:length(data[[i]]$results_bootstraps)) {
          se_mat[j+1,1:2] <- data[[i]]$results_bootstraps[[j]]$ci90
          se_mat[j+1,3:4] <- data[[i]]$results_bootstraps[[j]]$ci95
          att_vec[j+1]    <- data[[i]]$results_estimates[[j]]$ATT
          colors[[i]][j+1]  <- ifelse(data[[i]]$min_model == j, "#006284", 'gray40')
        }
        se_list[[i]] <- se_mat
        ATT_list[[i]] <- att_vec


      }

      ## plot
      points_master <- c(16:17, 15, 2, 4:9)
      locs <- seq(-0.3, 0.3, length.out = length(data[[1]]$results_bootstraps)+1)
      points <- c(1, points_master[1:length(data[[1]]$results_bootstraps)])

      if (is.null(xlim)) xlim <- c(0.5, length(ATT)+0.5)
      if (is.null(ylim)) ylim <- c(min(unlist(se_list)), max(unlist(se_list)))
      plot (1, 1, pch = 16, ylim = ylim, xlim = xlim, xaxt = "n", ylab = "ATT", xlab = "", ...)
      axis(1, at = 1:length(ATT), xlab)
      for (i in 1:length(data)) {
        for (j in 1:(length(data[[i]]$results_bootstraps)+1)) {
          # point
          points(i+locs[j], ATT_list[[i]][j], pch = points[j], col = colors[[i]][j])
          # 90%
          lines(c(i+locs[j],i+locs[j]), c(se_list[[i]][j,1], se_list[[i]][j,2]), col = colors[[i]][j],
                lwd = 2.2)
          # 95%
          arrows(i+locs[j], se_list[[i]][j,3], i+locs[j], se_list[[i]][j,4], col = colors[[i]][j],
            length=0.05, angle=90, code=3)

        }
      }
      lab_method <- ifelse(attr(data[[1]], 'method') == 'parametric', 'TFE', 'DiD')
      legend("topleft", legend = c(lab_method, paste("M", 1:length(data[[1]]$results_bootstraps))),
        pch = points, lty = 1, bty = 'n')

    } else if (isTRUE(attr(data[[1]], 'boot') & isTRUE(full))) {
      se_list <- ATT_list<- list()
      for (i in 1:length(data)) {
        se_mat <- matrix(NA, nrow = length(data[[i]]$results_bootstraps), ncol = 4)
        att_vec <- rep(NA, length(data[[i]]$results_bootstraps))
        for (j in 1:length(data[[i]]$results_bootstraps)) {
          se_mat[j,1:2] <- data[[i]]$results_bootstraps[[j]]$ci90
          se_mat[j,3:4] <- data[[i]]$results_bootstraps[[j]]$ci95
          att_vec[j] <- data[[i]]$results_estimates[[j]]$ATT
        }
        se_list[[i]] <- se_mat
        ATT_list[[i]] <- att_vec
      }

      ## plot
      points_master <- c(16:17, 1:9)
      locs <- seq(-0.3, 0.3, length.out = length(data[[1]]$results_bootstraps))
      points <- points_master[1:length(data[[1]]$results_bootstraps)]

      if (is.null(xlim)) xlim <- c(0.5, length(ATT)+0.5)
      if (is.null(ylim)) ylim <- c(min(unlist(se_list)), max(unlist(se_list)))
      plot (1, 1, pch = 16, ylim = ylim, xlim = xlim, xaxt = "n", ylab = "ATT", xlab = "", ...)
      axis(1, at = 1:length(ATT), xlab)
      for (i in 1:length(data)) {
        for (j in 1:length(data[[i]]$results_bootstraps)) {
          # point
          points(i+locs[j], ATT_list[[i]][j], pch = points[j])
          # 90%
          lines(c(i+locs[j],i+locs[j]), c(se_list[[i]][j,1], se_list[[i]][j,2]), lwd = 2.2)
          # 95%
          arrows(i+locs[j], se_list[[i]][j,3], i+locs[j], se_list[[i]][j,4], length=0.05, angle=90, code=3)

        }
      }
      legend("topleft", legend = paste("M", 1:length(data[[1]]$results_bootstraps)),
        pch = points, lty = 1, bty = 'n')

    } else if (isTRUE(attr(data[[1]], 'boot'))) {
      se_mat <- matrix(NA, nrow = length(data), ncol = 4)
      for (i in 1:length(data)) {
        se_mat[i,1:2] <- data[[i]]$ci90
        se_mat[i,3:4] <- data[[i]]$ci95
      }

      ## plot
      if (is.null(xlim)) xlim <- c(0.5, length(ATT)+0.5)
      if (is.null(ylim)) ylim <- c(min(se_mat), max(se_mat))
      plot (ATT, pch = 16, ylim = ylim, xlim = xlim, xaxt = "n", xlab = "", ...)
      axis(1, at = 1:length(ATT), xlab)
      for (i in 1:length(data)) {
        # 90%
        lines(c(i,i), c(se_mat[i,1], se_mat[i,2]), lwd = 2.2)
        # 95%
        arrows(i, se_mat[i,3], i, se_mat[i,4], length=0.05, angle=90, code=3)
      }

    } else {
      stop("Plot function is not suppored for estimates without bootstrap se's")
    }
  }
}
