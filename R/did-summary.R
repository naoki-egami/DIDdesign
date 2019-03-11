
##
## Summary functions
##


#' Print output from summary function
#' @param obj an object of \code{summary.diddesign} class, typically an ouput from \code{\link{summary.diddesign}}.
#' @export
print.summary.diddesign <- function(obj) {
  cat("\nMethod: ", obj$method, "\n", sep = "")
  cat("\nCall: ", paste(deparse(obj$call), sep="\n", collapse = "\n"), "\n", sep = "")

  cat("\nMain:\n")
  print.default(obj$main, quote = FALSE, right = TRUE)
  cat("\n")

  if (!is.null(obj$results)) {
    cat("\nResults:\n\n")
    for (i in 1:length(obj$results)) {
      cat(" T = ", obj$results[[i]][["post_treat"]], "\n")
      print.default(obj$results[[i]][['results']], quote = FALSE, right = TRUE)
      cat("\n")
    }

    cat("\nSelection:")
    cat(" ", paste("M", obj$selection[['model']], sep = ''), "is selected\n\n")
    print.default(obj$selection[['selection']], quote = FALSE, right = TRUE, digits = 3)
    cat("\n")

  }

  invisible(obj)

}

#' Print output from summary function
#' @param obj an object of \code{summary.diddesign_data} class, typically an ouput from \code{\link{summary.diddesign}}.
#' @export
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

#' helper function to summarize data
#' @param data an object of \code{diddesign_data} class, typically an output from \code{\link{did_data}}.
#' @param return a list consists of the following
#' \itemize{
#'   \item{ymean}{a matrix of mean outcome. The first row is for the treated and the second row is for the control.
#'                Columns correspond to time index (ascending order).}
#'   \item{id_time}{a vector of unique time index (ascending order).}
#' }
#' @keywords internal
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

  ## compute variance
  yvar <- apply(Y, 2, function(x) tapply(x, treat, var, na.rm = TRUE))
  yvar <- apply(yvar, 2, rev) ## first row is TREATED | second row is CONTROL

  return(list("ymean" = ymean, 'yvar' = yvar, 'id_time' = id_time))

}

#' Summarize output from \code{\link{did}} function
#' @param obj an object of \code{diddesign_data} or \code{diddesign} class.
#'  Prints estimated ATT along with 95\% confidence intervals when output of \code{\link{did}} is suppplied.
#' @param full a boolean. If \code{TRUE}, all estimates are returned. Default is \code{TRUE}.
#' @return a list of the following components:
#'   \item{method}{a method used to estimate ATT, either \code{parametric} or \code{nonparametric}}
#'   \item{call}{a formula used in estimation.}
#'   \item{main}{a matrix of selected estimates.}
#'   \item{results}{a list of result matrices.
#'     Each element of list corresponds to a post-period. Returned only when \code{full = TRUE}.}
#'   \item{selection}{a matrix of statistics used to select the model. Returned only when \code{full = TRUE}.}
#' @family main functions
#' @export
summary.diddesign <- function(obj, full = FALSE) {
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

  } else if ('diddesign' %in% class(obj)){
    ## summary function
    res_tab <- generate_tab_parametric(obj, full = full)
    res_tab[['method']] <- attr(obj, 'method')
    class(res_tab) <- 'summary.diddesign'
  }

  return(res_tab)
}



#' Generate Table for Parametric Method
#' @param obj an object of \code{diddesign} class.
#' @param full a boolean; if \code{TRUE} all results are attached.
#' @keywords internal
generate_tab_parametric <- function(obj, full = FALSE) {
  res_tab <- list()

  ## formula
  res_tab$call <- attr(obj, 'call')

  ##
  ## main table
  ##

  ## get estimated ATT
  ATT <- sapply(obj, function(x) round(x$ATT, 3))

  ## obtain estiamted se for all years (for selected models )
  se_save <- sapply(1:length(ATT), function(i) {
    tmp <- paste(formatC(round(obj[[i]]$ci95, 3), format = 'f', digits = 3), collapse = ', ')
    paste("[", tmp, "]", sep = '')
  })

  ## get selected models
  selected <- paste("M", sapply(obj, function(x) as.character(x$min_model)), sep = '')

  ## make a table and add col labels
  main_tab <- cbind(ATT, se_save, selected)
  rownames(main_tab) <- sapply(obj, function(x) as.character(attr(x, 'post_treat')))
  colnames(main_tab) <- c("ATT", "95% Conf. Int.", "Selected")

  ## add to return table
  res_tab$main <- main_tab

  ##
  ## attach full results
  ##
  if (isTRUE(full)) {
    results_list <- list()

    ## get DID estimates
    DiD <- sapply(obj, function(x) round(x$results_standardDiD$ATT, 3))
    DiD_se_save <- sapply(1:length(ATT), function(i) {
      paste("[",
        paste(round(obj[[i]]$results_standardDiD$results_variance$ci95, 3), collapse = ', '), "]", sep = '')
    })

    for (i in 1:length(obj)) {
      ## save obj
      results <- matrix(NA, nrow = length(obj[[1]]$results_estimates)+1, ncol = 3)

      ## save DID
      results[1, 1] <- DiD[i]
      results[1, 2] <- DiD_se_save[i]


      ## se
      ci_tmp <- sapply(1:length(obj[[i]]$results_variance), function(j) {
        tmp <- paste(formatC(round(obj[[i]]$results_variance[[j]]$ci95, 3), format = 'f', digits = 3),
                    collapse = ', ')
        paste("[", tmp , "]", sep = '')
      })

      ## ATT
      ATT <- sapply(obj[[i]]$results_estimates, function(x) round(x$ATT, 3))

      ## selected model
      check_mark <- rep("", length(obj[[1]]$results_estimates)+1)
      check_mark[attr(obj, 'selection')$min_model+1] <- "*"

      ## save
      results[-1, 1] <- ATT
      results[-1, 2] <- ci_tmp
      results[,3]    <- check_mark

      ## add labels
      colnames(results) <- c("ATT", "95% Conf. Int.", "")
      rownames(results) <- c("DiD", paste("M", 1:length(ATT), sep = ''))
      results_list[[i]] <- list('results' = results, 'post_treat' = attr(obj[[i]], 'post_treat'))
    }

    res_tab$results <- results_list

    ##
    ## save selection
    ##
    selection <- matrix(NA, nrow = length(attr(obj, "selection")$test_theta), ncol = 2)
    selection[,1] <- attr(obj, "selection")$test_theta
    selection[,2] <- attr(obj, "selection")$test_se
    colnames(selection) <- c("Theta", 'Std. Error')
    rownames(selection) <- paste("M", length(attr(obj, "selection")$test_theta):1, sep = "")
    res_tab$selection <- list("selection" = selection, "model" = attr(obj, "selection")$min_model)
  }

  return(res_tab)
}
