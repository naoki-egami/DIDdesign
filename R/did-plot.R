
#' DiD plot
#'
#' Generate a plot for difference-in-differences design.
#' @param formula \code{outcome ~ treatment}.
#' @param data data set
#' @param post_treatment a vector of time index indicating post treatment periods.
#' @param id_subject a variable name of unit index.
#'  This should be left as \code{NULL} when repeated cross-section data is used.
#' @param id_time a variable name of time index.
#' @param diff_order Order of differences. \code{diff_order = 0} generates a plot based on the original outcome.
#'  \code{diff_order = 1} generates a plot based on differenced outcomes.
#' @param xlim xlim
#' @param ylim ylim
#' @param col col
#' @param loc Location of legend. See \code{\link{legend}}.
#' @param lwd Line width.
#'
#' @examples
#' # load package
#' require(DIDdesign)
#'
#' # load data
#' data(anzia2012)
#'
#' # original series
#' did_plot(lnavgsalary_cpi ~ oncycle, data = anzia2012,
#'   post_treatment = c(2007, 2008, 2009),
#'   id_subject = "district",
#'   id_time = "year",
#'   ylim = c(10.55, 10.8),
#'   diff_order = 0,
#'   main = 'Original Series'
#' )
#'
#' # differenced series
#' did_plot(lnavgsalary_cpi ~ oncycle, data = anzia2012,
#'   post_treatment = c(2007, 2008, 2009),
#'   id_subject = "district",
#'   id_time = "year",
#'   ylim = c(-0.05, 0.05),
#'   diff_order = 1,
#'   main = 'Differenced Series'
#' )
#' @importFrom dplyr %>% pull select tbl_df
#' @importFrom Formula as.Formula
#' @importFrom utils getFromNamespace
#' @export
did_plot <- function(formula, data, post_treatment, id_subject = NULL, id_time,
  diff_order = 0, xlim = NULL, ylim = NULL, col = NULL, loc = NULL, lwd = NULL, ...
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
    message("treat data as repeated cross-section data")
    dat_use <- data %>% select_(outcome, treatment, id_time) %>% tbl_df()
    colnames(dat_use) <- c("outcome", "treatment", "id_time")
    dat_use %>% group_by(treatment, id_time) %>%
      summarise(ymean = mean(outcome, na.rm = TRUE)) -> y_summary
    y1mean <- y_summary %>% filter(treatment == 1) %>% pull(ymean)
    y0mean <- y_summary %>% filter(treatment == 0) %>% pull(ymean)
    dat_use <- list("y1mean" = y1mean, 'y0mean' = y0mean,
                    'id_time' = sort(unique(y_summary$id_time)))
    attr(dat_use, 'post_treat') <- min(post_treatment)
    # if (is.null(xlim))
    # if (is.null(xlim))
    # plot(1, 1, type = 'n', xlim = xlim, ylim = ylim, col = col, ...)
    plot_diddesign_data(dat_use, panel = FALSE, diff_order = diff_order,
      xlim = xlim, ylim = ylim, col = col, loc = loc, lwd = lwd, ...)

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
    plot_diddesign_data(dat_use, panel = TRUE, diff_order = diff_order,
      xlim = xlim, ylim = ylim, col = col, loc = loc, lwd = lwd, ...)

  }

}

#' Baseline plot function for did_plot
#' @param data Two possible inputs. If \code{panel = TRUE}, input should be a \code{diddesign_data} class object.
#'    If \code{panel = FALSE}, it should be a list consists of \code{"y1mean"}, \code{"y0mean"} and \code{"id_time"}.
#' @param panel A boolean argument. \code{TRUE} if data are in the panel format and \code{FALSE} if repeated cross-section.
#' @param diff_order Either \code{0} (original series) or \code{1} (differenced series).
#' @keywords internal
plot_diddesign_data <- function(data, panel = TRUE, diff_order = 0,
  xlim = NULL, ylim = NULL, col = NULL, loc = NULL, lwd = NULL, ...) {

  if(isTRUE(panel)) {
    summary_dat <- summarize.diddesign(data)
    ymean <- summary_dat$ymean
    id_time <- summary_dat$id_time
    time_use <- which(attr(data[[1]], 'post_treat') == id_time)

  } else {
    ## data
    ymean <- rbind(data$y1mean, data$y0mean)
    id_time <- data$id_time
    time_use <- which(attr(data, 'post_treat') == id_time)

  }

  ## plot
  if (diff_order == 0) {
    if (is.null(xlim)) xlim <- c(min(id_time), max(id_time))
    if (is.null(ylim)) ylim <- c(min(ymean), max(ymean))
    if (is.null(col) | length(col) == 1) col <- c(col, '#006284', 'gray40')
    if (is.null(lwd)) lwd <- 1.5
    if (is.null(loc)) loc <- 'topleft'

    ## background color setup
    bg_col <- rgb(215, 196, 187, maxColorValue = 255, alpha = 80)

    ub_inc <- ifelse(max(ylim) > 0, 1.10, 1/1.10)
    lb_inc <- ifelse(min(ylim) > 0, 1/1.10, 1.10)
    ## plot
    plot(1, 1, type = 'n', ylim = ylim, xlim = xlim, xlab = "", ylab = "Mean Outcome", ...)
    rect(mean(id_time[(time_use-1):time_use]), min(ylim) * lb_inc,
          max(xlim) * 1.04, max(ylim) * ub_inc, border = NA, col = bg_col)
    abline(v = mean(id_time[(time_use-1):time_use]), col = '#B9887D', lty = 3, lwd = 1.3)
    lines(id_time, ymean[1,], col = col[1], pch = 16, type = 'b', lwd = lwd)
    lines(id_time, ymean[2,], col = col[2],  pch = 17, type = 'b', lwd = lwd)
    legend(loc, legend = c("treated", 'control'), col = col[1:2],
      lty = 1, pch = c(16, 17), bty = 'n')
    box(lwd = 1.2)
  } else {
    ## diff_order = 2

    ## prepare data
    ymean_diff <- t(apply(ymean, 1, diff))

    if (is.null(xlim)) xlim <- c(1, length(id_time)-1) # min(id_time), max(id_time))
    if (is.null(ylim)) ylim <- c(min(ymean_diff), max(ymean_diff))
    if (is.null(col) | length(col) == 1) col <- c(col, '#006284', 'gray40')
    if (is.null(lwd)) lwd <- 1.5
    if (is.null(loc)) loc <- 'topleft'

    time_use <- which(attr(data[[1]], 'post_treat') == id_time)
    time_lab <- sapply(1:(length(id_time)-1), function(i) {
      paste(id_time[i+1], "-", id_time[i], sep = '')
    })
    ## background color setup
    bg_col <- rgb(215, 196, 187, maxColorValue = 255, alpha = 80)

    ub_inc <- ifelse(max(ylim) > 0, 1.10, 1/1.10)
    lb_inc <- ifelse(min(ylim) > 0, 1/1.10, 1.10)

    ## plot
    plot(1, 1, type = 'n', ylim = ylim, xlim = xlim, xlab = "",
      ylab = "Differenced Mean Outcome", xaxt = "n", ...)
    axis(1, at = 1:max(xlim), labels = time_lab, cex.axis = 0.8)
    rect(time_use - 1.5, min(ylim) * lb_inc, max(xlim) * 1.04, max(ylim) * ub_inc,
        border = NA, col = bg_col)
    abline(v = time_use - 1.5, col = '#B9887D', lty = 3, lwd = 1.3)
    lines(ymean_diff[1,], col = col[1], pch = 16, type = 'b', lwd = lwd)
    lines(ymean_diff[2,], col = col[2],  pch = 17, type = 'b', lwd = lwd)
    legend(loc, legend = c("treated", 'control'), col = col[1:2],
      lty = 1, pch = c(16, 17), bty = 'n')
    box(lwd = 1.2)

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
  # if('diddesign_data' %in% class(data)) {
  #   ##
  #   ## plot for raw data
  #   ##
  #
  #   summary_dat <- summarize.diddesign(data)
  #   ymean <- summary_dat$ymean
  #   id_time <- summary_dat$id_time
  #
  #   ## plot
  #   if (is.null(xlim)) xlim <- c(min(id_time), max(id_time))
  #   if (is.null(ylim)) ylim <- c(min(ymean), max(ymean))
  #   if (is.null(col) | length(col) == 1) col <- c(col, '#006284', 'gray60')
  #   if (is.null(lwd)) lwd <- 1.5
  #
  #   time_use <- which(attr(data[[1]], 'post_treat') == id_time)
  #
  #   ## plot
  #   plot(1, 1, type = 'n', ylim = ylim, xlim = xlim, xlab = "", ylab = "Mean Outcome", ...)
  #   abline(v = mean(id_time[(time_use-1):time_use]), col = 'red', lty = 3, lwd = 1.3)
  #   lines(id_time, ymean[1,], col = col[1], pch = 16, type = 'b', lwd = lwd)
  #   lines(id_time, ymean[2,], col = col[2],  pch = 17, type = 'b', lwd = lwd)
  #   legend('topleft', legend = c("treated", 'control'), col = col[1:2],
  #     lty = 1, pch = c(16, 17), bty = 'n')
  #
  # } else

  if ('diddesign' %in% class(data)) {
    ##
    ## plot for estimation result
    ##
    ATT <- sapply(data, function(x) round(x$ATT, 3))
    xlab <- sapply(data, function(x) attr(x, 'post_treat'))
    if (isTRUE(attr(data[[1]], 'boot') & isTRUE(full)) &
        !is.null(data[[1]]$results_standardDiD)) {
      se_list <- ATT_list <- colors <- list()
      for (i in 1:length(data)) {
        se_mat <- matrix(NA, nrow = length(data[[i]]$results_variance) + 1, ncol = 4)
        att_vec <- rep(NA, length(data[[i]]$results_variance) + 1)

        ## save DiD
        se_mat[1,1:2] <- data[[i]]$results_standardDiD$results_variance$ci90
        se_mat[1,3:4] <- data[[i]]$results_standardDiD$results_variance$ci95
        att_vec[1] <- data[[i]]$results_standardDiD$ATT

        colors[[i]] <- rep('gray60', length(data[[i]]$results_variance))
        for (j in 1:length(data[[i]]$results_variance)) {
          se_mat[j+1,1:2] <- data[[i]]$results_variance[[j]]$ci90
          se_mat[j+1,3:4] <- data[[i]]$results_variance[[j]]$ci95
          att_vec[j+1]    <- data[[i]]$results_estimates[[j]]$ATT
          colors[[i]][j+1]  <- ifelse(data[[i]]$min_model == j, "#006284", 'gray40')
        }
        se_list[[i]] <- se_mat
        ATT_list[[i]] <- att_vec


      }

      ## plot
      points_master <- c(16:17, 15, 2, 4:9)
      locs <- seq(-0.3, 0.3, length.out = length(data[[1]]$results_variance)+1)
      points <- c(1, points_master[1:length(data[[1]]$results_variance)])

      if (is.null(xlim)) xlim <- c(0.5, length(ATT)+0.5)
      if (is.null(ylim)) ylim <- c(min(unlist(se_list)), max(unlist(se_list)))
      plot (1, 1, pch = 16, ylim = ylim, xlim = xlim, xaxt = "n", ylab = "ATT", xlab = "", ...)
      axis(1, at = 1:length(ATT), xlab)
      for (i in 1:length(data)) {
        for (j in 1:(length(data[[i]]$results_variance)+1)) {
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
      legend("topleft", legend = c(lab_method, paste("M", 1:length(data[[1]]$results_variance))),
        pch = points, lty = 1, bty = 'n')

    } else if (isTRUE(attr(data[[1]], 'boot') & isTRUE(full))) {
      se_list <- ATT_list<- list()
      for (i in 1:length(data)) {
        se_mat <- matrix(NA, nrow = length(data[[i]]$results_variance), ncol = 4)
        att_vec <- rep(NA, length(data[[i]]$results_variance))
        for (j in 1:length(data[[i]]$results_variance)) {
          se_mat[j,1:2] <- data[[i]]$results_variance[[j]]$ci90
          se_mat[j,3:4] <- data[[i]]$results_variance[[j]]$ci95
          att_vec[j] <- data[[i]]$results_estimates[[j]]$ATT
        }
        se_list[[i]] <- se_mat
        ATT_list[[i]] <- att_vec
      }

      ## plot
      points_master <- c(16:17, 1:9)
      locs <- seq(-0.3, 0.3, length.out = length(data[[1]]$results_variance))
      points <- points_master[1:length(data[[1]]$results_variance)]

      if (is.null(xlim)) xlim <- c(0.5, length(ATT)+0.5)
      if (is.null(ylim)) ylim <- c(min(unlist(se_list)), max(unlist(se_list)))
      plot (1, 1, pch = 16, ylim = ylim, xlim = xlim, xaxt = "n", ylab = "ATT", xlab = "", ...)
      axis(1, at = 1:length(ATT), xlab)
      for (i in 1:length(data)) {
        for (j in 1:length(data[[i]]$results_variance)) {
          # point
          points(i+locs[j], ATT_list[[i]][j], pch = points[j])
          # 90%
          lines(c(i+locs[j],i+locs[j]), c(se_list[[i]][j,1], se_list[[i]][j,2]), lwd = 2.2)
          # 95%
          arrows(i+locs[j], se_list[[i]][j,3], i+locs[j], se_list[[i]][j,4], length=0.05, angle=90, code=3)

        }
      }
      legend("topleft", legend = paste("M", 1:length(data[[1]]$results_variance)),
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




#' Plot function for selection
#' @param data \code{diddesign} object.
#'    Typically, it is an output from \code{\link{did}} function.
#' @param alpha Level of rejection.
#'    Confidence intervals are constructed with the 1-alpha and 1 - alpha/2 levels.
#' @param equivalence boolean; if \code{TRUE}, equivalence region is superimposed on a plot.
#'    See also \code{eL} and \code{eU} arguments to set the upper and the lower bound of a region.
#' @param eU upper bound of equivalence region. If unspecified, the largest value of 1 - alpha level sets is used.
#' @param eL lower bound of equivalence region. If unspecified, the smallest value of 1 - alpha level sets is used.
#' @param xlim xlim. See \code{\link{plot}} function.
#' @param ylim ylim. See \code{\link{plot}} function.
#' @param loc location of legend. See \code{\link{legend}} function.
#' @examples
#' # load package
#' require(DIDdesign)
#'
#' # load data
#' data(anzia2012)
#'
#' # fit models
#' fit1 <- did(lnavgsalary_cpi ~ oncycle, data = anzia2012,
#'             id_subject = "district", id_time = "year",
#'             post_treatment = c(2007, 2008, 2009),
#'             method = "parametric")
#' fit2 <- did(lnavgsalary_cpi ~ oncycle, data = anzia2012,
#'             id_subject = "district", id_time = "year",
#'             post_treatment = c(2007, 2008, 2009),
#'             method = "nonparametric", se_boot = FALSE)
#'
#' # make a plot with equivalence region
#' par(mfrow = c(1,2), mar = c(4, 2.5, 3.5, 1))
#' did_plot_selection(fit1, equivalence = TRUE, eL = -0.005, eU = 0.005,
#'  ylim = c(-0.015, 0.015), main = "Parametric")
#' did_plot_selection(fit2, equivalence = TRUE, eL = -0.005, eU = 0.005,
#'  ylim = c(-0.015, 0.015), main = "Nonparametric")
#' @export
did_plot_selection <- function(
  data, alpha = 0.05, equivalence = TRUE, eL = NULL, eU = NULL,
  xlim = NULL, ylim = NULL, loc = NULL,...) {

  if (!(class(data) %in% c("diddesign"))) {
    stop("Input data should be a 'diddesign' object")
  }

  if(exists("selection", attributes(data))) {
    ## extract info
    selection <- attr(data, "selection")
    theta  <- selection$test_theta
    stderr <- selection$test_se
    model  <- selection$min_model

    ## construct confidence intervals
    CIa <- cbind(theta + qnorm(alpha/2) * stderr,
                 theta + qnorm(1 - alpha/2) * stderr)
    CI2a <- cbind(theta + qnorm(alpha) * stderr,
                  theta + qnorm(1 - alpha) * stderr)

    ## prepare plot parameters
    if (is.null(ylim)) ylim <- c(min(CIa), max(CIa))
    if (is.null(xlim)) xlim <- c(0.7, length(theta) + 0.3)
    if (is.null(loc)) loc <- 'topleft'

    col <- rep('gray30', length(theta))
    if (model <= length(theta)) col[model] <- "#006284"
    col <- rev(col)


    if (isTRUE(equivalence)) {
      if(is.null(eL)) eL <- min(CI2a)
      if(is.null(eU)) eU <- max(CI2a)
      message("Equivalence region: [", round(eL, 3), ",", round(eU, 3),
              "] at the alpha = ", alpha, " level.")
      eq_color <- rgb(123, 162, 63, maxColorValue = 255, alpha = 60)
    }

    ## make plot
    plot(theta, type = 'n', pch = 16, cex = 1.2, ylim = ylim, xlim = xlim,
        xaxt = 'n', xlab  = "", ylab = "Theta", col = col, ...)
    if (isTRUE(equivalence)) {
      rect(min(xlim) / 1.3, eL, max(xlim) * 1.3, eU, border = NA, col = eq_color)
    }
    abline(h = 0, col = 'gray60', lwd = 1.3, lty = 3)
    for (i in 1:nrow(CIa)) {
      lines(c(i,i), CI2a[i,], lwd = 2.2, col = col[i])
      arrows(i, CIa[i,1], i, CIa[i,2], col = col[i], length = 0.05, angle = 90, code = 3)
    }
    points(theta, pch = 16, cex = 1.2, col = col)
    axis(1, at = 1:length(theta), labels = paste("M", rev(1:length(theta)), sep = ''))
    box(lwd = 1.25)
    if (isTRUE(equivalence) & (model <= length(theta))) {
    legend(loc, legend = c('selected model', 'equivalence region'),
      lty = c(1, NA), pch = c(16, NA), col = c("#006284", NA),
      fill = c(NA, eq_color), border = c(NA, 'black'), bty = 'n')
    } else if (!isTRUE(equivalence) & (model <= length(theta))) {
      legend(loc, legend = c('selected model'), lty = 1, pch = 16, col = c("#006284"), bty = 'n')
    } else if (isTRUE(equivalence) & !(model <= length(theta))) {
      legend(loc, legend = c( 'equivalence region'),
        fill = c(eq_color), border = c('black'), bty = 'n')
    } else {
      legend(loc, legend = c('selection statistic'),
        lty = 1, pch = 16, col = c("gray30"), bty = 'n')

    }
  } else {
    stop("selection does not exsit.")
  }
}
