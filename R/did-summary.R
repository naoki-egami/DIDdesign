
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
    selected <- paste("M", sapply(obj, function(x) as.character(x$min_model)), sep = '')

    if (isTRUE(attr(obj[[1]], 'boot'))) {
      se_save <- rep(NA, length(ATT))
      for (i in 1:length(ATT)) {
        se_save[i] <- paste("[", paste(round(obj[[i]]$ci95, 3), collapse = ', '), "]", sep = '')
      }

      # did estimates
      if (!is.null(obj[[1]]$results_standardDiD) & attr(obj[[1]], 'method') == 'nonparametric') {
        # BIC <- sapply(obj, function(x) round(x$BIC_min, 3))
        # HQIC <- sapply(obj, function(x) round(x$HQIC_min, 3))

        DiD <- sapply(obj, function(x) round(x$results_standardDiD$ATT, 3))
        DiD_se_save <- rep(NA, length(DiD))
        for (i in 1:length(ATT)) {
          DiD_se_save[i] <- paste("[",
            paste(round(obj[[i]]$results_standardDiD$results_bootstraps$ci95, 3),
            collapse = ', '), "]", sep = '')
        }

        # make a table
        tabs_var <- c("D-DiD", "", "", "Std-DiD", "")
        labs_var <- c("ATT", "95% CI", "Selected", "ATT",  "95% CI")
        res_tab <- data.frame(cbind(tabs_var, labs_var,
          rbind(ATT, se_save, selected, DiD, DiD_se_save))
        )

        # col labels
        colnames(res_tab) <- c("", "", sapply(obj, function(x) attr(x, 'post_treat')))
        rownames(res_tab) <- NULL
      } else if(attr(obj[[1]], 'method') == 'nonparametric') {
        # BIC <- sapply(obj, function(x) round(x$BIC_min, 3))
        # HQIC <- sapply(obj, function(x) round(x$HQIC_min, 3))

        res_tab <- data.frame(rbind(
          ATT, se_save, selected
        ), row.names = c("ATT", "95% CI", "Selected Model"))
        colnames(res_tab) <- sapply(obj, function(x) attr(x, 'post_treat'))

      } else if (!is.null(obj[[1]]$results_standardDiD) & attr(obj[[1]], 'method') == 'parametric') {
        DiD <- sapply(obj, function(x) round(x$results_standardDiD$ATT, 3))
        DiD_se_save <- rep(NA, length(DiD))
        for (i in 1:length(ATT)) {
          DiD_se_save[i] <- paste("[",
            paste(round(obj[[i]]$results_standardDiD$results_bootstraps$ci95, 3),
            collapse = ', '), "]", sep = '')
        }

        # make a table
        tabs_var <- c("D-DiD", "", "", "2way-FE", "")
        labs_var <- c("ATT", "95% CI", "Selected", "ATT",  "95% CI")
        res_tab <- data.frame(cbind(tabs_var, labs_var,
          rbind(ATT, se_save, selected, DiD, DiD_se_save))
        )

        # col labels
        colnames(res_tab) <- c("", "", sapply(obj, function(x) attr(x, 'post_treat')))
        rownames(res_tab) <- NULL
      } else {
        res_tab <- data.frame(rbind(
          ATT, se_save, selected
        ), row.names = c("ATT", "95% CI", "Selected Model"))
        colnames(res_tab) <- sapply(obj, function(x) attr(x, 'post_treat'))

      }

    } else {
      res_tab <- data.frame(rbind(
        ATT, selected
      ), row.names = c("ATT", "Selected Model"))
      colnames(res_tab) <- sapply(obj, function(x) attr(x, 'post_treat'))

    }
  }

  return(res_tab)
}
