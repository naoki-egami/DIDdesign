
##
## Summary and plot functions
##

summarize.diddesign <- function(data) {
  n_post <- length(data)
  ## get some data
  treat  <- data[[1]]$D
  n_obs  <- length(treat)
  t_time <- ncol(data[[1]]$Y)
  Y_pre  <- data[[1]]$Y[,-t_time]
  id_time <- attr(data[[1]], 'id_time')

  if (n_post > 1) {
    Y_post <- matrix(NA, nrow = n_obs, ncol = n_post)
    for (i in 1:n_post) {
      Y_post[,i] <- data[[i]]$Y[,t_time]
      id_time <- c(id_time, attr(data[[i]], 'post_treat'))
    }
    id_time <- unique(id_time)
  } else {
    Y_post <- data[[1]]$Y[,t_time]
  }

  ## compute mean outcome
  Y <- cbind(Y_pre, Y_post)
  ymean <- apply(Y, 2, function(x) tapply(x, treat, mean, na.rm = TRUE))
  ymean <- apply(ymean, 2, rev) ## first row is TREATED | second row is CONTROL
  return(list("ymean" = ymean, 'id_time' = id_time))

}

#' Summary Function
#' @export
summary.diddesign <- function(obj) {
  if ('diddesign_data' %in% class(obj)) {
    summary_dat <- summarize.diddesign(obj)
    post_first <- attr(obj[[1]], 'post_treat')
    post_bool  <- summary_dat$id_time >= post_first
    status <- ifelse(post_bool, "T", "C")
    res_tab <- data.frame(rbind(
      status, round(summary_dat$ymean[1,], 2), round(summary_dat$ymean[2,], 2),
      round(summary_dat$ymean[1,] - summary_dat$ymean[2,], 2)
    ), row.names =  c("Status", "Y:Treated", "Y:Control", "YT-YC"))

    colnames(res_tab) <- summary_dat$id_time
    # rownames(res_tab) <- c("Status", "Y:Treated", "Y:Control")
  } else if ('diddesign' %in% class(obj)){
    ATT <- sapply(obj, function(x) round(x$ATT, 3))
    BIC <- sapply(obj, function(x) round(x$BIC_min, 3))
    HQIC <- sapply(obj, function(x) round(x$HQIC_min, 3))
    selected <- paste("M", sapply(obj, function(x) as.character(x$min_model)), sep = '')

    if (isTRUE(attr(obj[[1]], 'boot'))) {
      se_save <- rep(NA, length(ATT))
      for (i in 1:length(ATT)) {
        se_save[i] <- paste("[", paste(round(obj[[i]]$ci95, 3), collapse = ', '), "]", sep = '')
      }

      # did estimates
      if (!is.null(obj[[1]]$results_standardDiD)) {
        DiD <- sapply(obj, function(x) round(x$results_standardDiD$ATT, 3))
        DiD_se_save <- rep(NA, length(DiD))
        for (i in 1:length(ATT)) {
          DiD_se_save[i] <- paste("[",
            paste(round(obj[[i]]$results_standardDiD$results_bootstraps$ci95, 3),
            collapse = ', '), "]", sep = '')
        }

        # make a table
        tabs_var <- c("D-DiD", "", "", "", "", "Std-DiD", "")
        labs_var <- c("ATT", "95% CI", "BIC", "HQIC", "Selected", "ATT",  "95% CI")
        res_tab <- data.frame(cbind(tabs_var, labs_var,
          rbind(ATT, se_save, BIC, HQIC, selected, DiD, DiD_se_save))
        )

        # col labels
        colnames(res_tab) <- c("", "", sapply(obj, function(x) attr(x, 'post_treat')))
        rownames(res_tab) <- NULL
      } else {
        res_tab <- data.frame(rbind(
          ATT, se_save, BIC, HQIC, selected
        ), row.names = c("ATT", "95% CI", "BIC", "HQIC", "Selected Model"))
        colnames(res_tab) <- sapply(obj, function(x) attr(x, 'post_treat'))
      }

    } else {
      res_tab <- data.frame(rbind(
        ATT, BIC, HQIC, selected
      ), row.names = c("ATT", "BIC", "HQIC", "Selected Model"))
      colnames(res_tab) <- sapply(obj, function(x) attr(x, 'post_treat'))

    }
  }

  return(res_tab)
}

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
#' @export
did_plot <- function(outcome, treatment, post_treatment, id_subject, id_time,
  xlim = NULL, ylim = NULL, col = NULL, ...
) {
  
  ## FORMAT DATA 
  dat <- did_data(
    outcome = outcome,
    treatment = treatment,
    post_treatment = post_treatment,
    id_subject = id_subject,
    id_time = id_time,
    long = TRUE
  )
  
  ## CALL plot.diddesign
  # par(mfrow = c(1, 3))
  plot(dat, xlim = xlim, ylim = ylim, col = col, ...)

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
      se_list <- ATT_list<- list()
      for (i in 1:length(data)) {
        se_mat <- matrix(NA, nrow = length(data[[i]]$results_bootstraps) + 1, ncol = 4)
        att_vec <- rep(NA, length(data[[i]]$results_bootstraps) + 1)

        ## save DiD
        se_mat[1,1:2] <- data[[i]]$results_standardDiD$results_bootstraps$ci90
        se_mat[1,3:4] <- data[[i]]$results_standardDiD$results_bootstraps$ci95
        att_vec[1] <- data[[i]]$results_standardDiD$ATT

        for (j in 1:length(data[[i]]$results_bootstraps)) {
          se_mat[j+1,1:2] <- data[[i]]$results_bootstraps[[j]]$ci90
          se_mat[j+1,3:4] <- data[[i]]$results_bootstraps[[j]]$ci95
          att_vec[j+1]    <- data[[i]]$results_estimates[[j]]$ATT
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
          points(i+locs[j], ATT_list[[i]][j], pch = points[j])
          # 90%
          lines(c(i+locs[j],i+locs[j]), c(se_list[[i]][j,1], se_list[[i]][j,2]), lwd = 2.2)
          # 95%
          arrows(i+locs[j], se_list[[i]][j,3], i+locs[j], se_list[[i]][j,4], length=0.05, angle=90, code=3)

        }
      }
      legend("topleft", legend = c("DiD", paste("M", 1:length(data[[1]]$results_bootstraps))),
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
