
##
## Summary and plot functions
##


# print function 
print.summary.diddesign <- function(obj) {
  cat("\nCall:\n", paste(deparse(obj$call), sep="\n", collapse = "\n"), "\n\n", sep = "")
  
  cat("\nMain:\n")
  print(obj$main, quote = FALSE, right = TRUE)
  cat("\n")
  invisible(obj)
  
}

print.summary.diddesign_data <- function(obj) {
  cat("\nSummary:\n")
  res_tab <- rbind(
    obj[["Status"]], obj[["Y:Treated"]], obj[["Y:Control"]], obj[["YT-YC"]]
  )

  ##
  colnames(res_tab) <- obj$id_time
  rownames(res_tab) <- c("Status", "Y:Treated", "Y:Control", "YT-YC")
  
  ##
  print.default(res_tab, quote = FALSE, right = TRUE)
  cat("\n")
  invisible(obj)
}

# helper function 
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
summary.diddesign <- function(obj, full = TRUE) {
  if ('diddesign_data' %in% class(obj)) {
    summary_dat <- summarize.diddesign(obj)
    post_first <- attr(obj[[1]], 'post_treat')
    post_bool  <- summary_dat$id_time >= post_first
    status <- ifelse(post_bool, "T", "C")
    
    res_tab <- list() 
    res_tab[["Status"]]    <- status 
    res_tab[["Y:Treated"]] <- round(summary_dat$ymean[1,], 3)
    res_tab[["Y:Control"]] <- round(summary_dat$ymean[2,], 3)
    res_tab[["YT-YC"]]     <- round(summary_dat$ymean[1,] - summary_dat$ymean[2,], 3)
    res_tab[["id_time"]]   <- summary_dat$id_time
    class(res_tab) <- 'summary.diddesign_data'
    
    # rownames(res_tab) <- c("Status", "Y:Treated", "Y:Control")
  } else if ('diddesign' %in% class(obj)){
    
    ## get estimated ATT 
    ATT <- sapply(obj, function(x) round(x$ATT, 3))
    ## get selected models 
    selected <- paste("M", sapply(obj, function(x) as.character(x$min_model)), sep = '')
    ## obtain estiamted se for all years (for selected models )
    se_save <- sapply(1:length(ATT), function(i) {
      tmp <- paste(formatC(round(obj[[i]]$ci95, 3), format = 'f', digits = 3), collapse = ', ')
      paste("[", tmp, "]", sep = '')
    })
    
    ## get information about fitted method 
    is_parametric <- attr(obj[[1]], 'method') == 'parametric'
    
    ## summary function 
    if (isTRUE(is_parametric)) {
      res_tab <- generate_tab_parametric(obj, ATT, se_save, selected, full = full)
      class(res_tab) <- 'summary.diddesign'      
    } else {
      DiD <- sapply(obj, function(x) round(x$results_standardDiD$ATT, 3))
      DiD_se_save <- rep(NA, length(DiD))
      for (i in 1:length(ATT)) {
        DiD_se_save[i] <- paste("[",
          paste(round(obj[[i]]$results_standardDiD$results_variance$ci95, 3),
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
      class(res_tab) <- 'summary.diddesign'  
    }
  }
  
  return(res_tab)
}



#' Generate Table for Parametric Method 
#' @keywords internal
generate_tab_parametric <- function(obj, ATT, se_save, selected, full = TRUE) {
  res_tab <- list()
  
  ## formula 
  res_tab$call <- attr(obj, 'call')
  
  ## get did estimates 
  DiD <- sapply(obj, function(x) round(x$results_standardDiD$ATT, 3))
  DiD_se_save <- sapply(1:length(ATT), function(i) {
    paste("[", 
      paste(round(obj[[i]]$results_standardDiD$results_variance$ci95, 3), collapse = ', '), "]", sep = '')
  })

  ## make a main table
  # tabs_var <- c("D-DiD", "", "", "2way-FE", "")
  # labs_var <- c("ATT", "95% CI", "Selected", "ATT",  "95% CI")
  # main_tab <- data.frame(cbind(tabs_var, labs_var,
  #   rbind(ATT, se_save, selected, DiD, DiD_se_save))
  # )
  # 
  # ## col labels
  # colnames(main_tab) <- c("", "", sapply(obj, function(x) attr(x, 'post_treat')))
  # rownames(main_tab) <- NULL    

  main_tab <- cbind(ATT, se_save, selected)
  
  ## col labels
  rownames(main_tab) <- sapply(obj, function(x) as.character(attr(x, 'post_treat')))
  colnames(main_tab) <- c("ATT", "95% Conf. Int.", "Selected")

  ## add 
  res_tab$main <- main_tab
  
  ## attach full results
  if (isTRUE(full)) {
    for (i in 1:length(obj)) {
      
      ## se 
      ci_tmp <- sapply(1:length(obj[[i]]$results_variance), function(j) {
        tmp <- paste(round(obj[[i]]$results_variance[[j]]$ci95, 3), collapse = ', ')
        paste("[", tmp , "]", sep = '')                
      })
      
      ## ATT 
      ATT <- sapply(obj[[i]]$results_estimates, function(x) round(x$ATT, 3))
      
      ## 
      full_tab <- data.frame(rbind(ATT, ci_tmp))
      rownames(full_tab) <- c("ATT", "95% Conf. Int.")
      colnames(full_tab) <- paste("M", 1:ncol(full_tab), sep = '')
      year_name <- as.character(attr(obj[[i]], 'post_treat'))
      res_tab[[year_name]] <- full_tab
    }
  }
  
  return(res_tab)
}
