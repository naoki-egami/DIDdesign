
##
## Summary functions
##


#' Print output from summary function
#' @param obj an object of \code{summary.diddesign} class, typically an ouput from \code{\link{summary.diddesign}}.
#' @importFrom cli cat_rule
#' @export
print.summary.diddesign <- function(obj) {
  cat_rule(left = crayon::bold("Methods"))
  cat("Method: ", obj$method, "\n", sep = "")
  cat("\nCall: ", paste(deparse(obj$call), sep="\n", collapse = "\n"), "\n", sep = "")

  # cat("\nMain:\n")
  cat("\n")
  cat_rule(left = crayon::bold("Main Result"))
  print.default(obj$main, quote = FALSE, right = TRUE, digits = 3)

  if (!is.null(obj$results)) {
    cat("\n")
    cat_rule(left = crayon::bold("Full Results"), line = 2)
    for (i in 1:length(obj$results)) {
      cat(" T = ", obj$results[[i]][["post_treat"]], "\n")
      tab_print  <- round(obj$results[[i]][['estimates']], 3)
      tab_print[,5] <- ifelse(tab_print[,5] == 1, cli::symbol$tick, "")
      print.default(tab_print, quote = FALSE, right = TRUE, digits = 3)
      cat("\n")
    }

    cat_rule(left = crayon::bold("Selection"), line = 2)
    cat(crayon::bold('Bias Test: '), "")
    cat("  ", paste("M", obj$selection[['model']], sep = ''), "is selected\n\n")
    print.default(obj$selection[['selection']], quote = FALSE, right = TRUE, digits = 3)
    cat("\n")
    cat(crayon::bold("Bias test (equivalence based): "), "\n\n")
    
    names(obj$selection[['equivalence']]) <- rev(rownames(obj$selection[['selection']]))
    
    print.default(obj$selection[['equivalence']], quote = FALSE, right = TRUE)
    bias <- round(abs(obj$selection[['bias']]), 3)
    cat("\nEquivalence region is ", "[", -bias, ",", bias, "]\n")
    cat("\n")
    cat(crayon::bold("Sign test for PTT: "), obj$selection[['sign']], "( p-value =", round(obj$selection[['pval']], 3), ")\n")
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
  # ATT <- sapply(obj, function(x) round(x$ATT, 3))
  ATT <- sapply(obj, function(x) x$ATT)

  ## obtain estiamted se for all years (for selected models )
  ci_save <- t(sapply(1:length(ATT), function(i) {
    # tmp <- paste(formatC(round(obj[[i]]$ci95, 3), format = 'f', digits = 3), collapse = ', ')
    # paste("[", tmp, "]", sep = '')
    tmp <- obj[[i]]$ci95
    tmp
  }))

  se_save <- sapply(1:length(ATT), function(i) { obj[[i]]$se } )
  
  ## get selected models
  # selected <- paste("M", sapply(obj, function(x) as.character(x$min_model)), sep = '')
  selected <- sapply(obj, function(x) x$min_model)

  ## make a table and add col labels
  main_tab <- cbind(ATT, se_save, ci_save[,1], ci_save[,2], selected)
  rownames(main_tab) <- sapply(obj, function(x) as.character(attr(x, 'post_treat')))
  colnames(main_tab) <- c("ATT", "SE", "95% CI (LB)", "95% CI (UB)", "Selected")

  ## add to return table
  res_tab$main <- main_tab

  ##
  ## attach full results
  ##
  if (isTRUE(full)) {
    results_list <- list()

    ## get DID estimates
    DiD <- sapply(obj, function(x) round(x$results_standardDiD$ATT, 3))
    DiD_ci_save <- t(sapply(1:length(ATT), function(i) {
      # paste("[",
      #   paste(round(obj[[i]]$results_standardDiD$results_variance$ci95, 3), collapse = ', '), "]", sep = '')
      obj[[i]]$results_standardDiD$results_variance$ci95
    }))
    
    DiD_se_save <- sapply(1:length(ATT), function(i) obj[[i]]$results_standardDiD$results_variance$se)

    for (i in 1:length(obj)) {
      ## save obj
      results <- matrix(0, nrow = length(obj[[1]]$results_estimates)+1, ncol = 5)

      ## save DID
      results[1, 1] <- DiD[i]
      results[1, 2] <- DiD_se_save[i]
      results[1, 3:4] <- DiD_ci_save[i,]


      ## se
      ci_tmp <- t(sapply(1:length(obj[[i]]$results_variance), function(j) {
        obj[[i]]$results_variance[[j]]$ci95
        # tmp <- paste(formatC(round(obj[[i]]$results_variance[[j]]$ci95, 3), format = 'f', digits = 3),
        #             collapse = ', ')
        # paste("[", tmp , "]", sep = '')
      }))

      se_tmp <- sapply(1:length(obj[[i]]$results_variance), function(j) {
        obj[[i]]$results_variance[[j]]$se
      })
      
      ## ATT
      ATT <- sapply(obj[[i]]$results_estimates, function(x) round(x$ATT, 3))

      ## selected model
      # check_mark <- rep("", length(obj[[1]]$results_estimates)+1)
      check_mark <- rep(0, length(obj[[1]]$results_estimates)+1)
      check_mark[attr(obj, 'selection')$min_model+1] <- 1 # cli::symbol$tick

      ## save
      results[-1, 1] <- ATT
      results[-1, 2] <- se_tmp
      results[-1, 3:4] <- ci_tmp
      results[,5]    <- check_mark

      ## add labels
      colnames(results) <- c("ATT", "SE", "95% CI (LB)", "95% CI (UB)", "")
      rownames(results) <- c("DiD", paste("M", 1:length(ATT), sep = ''))
      results_list[[i]] <- list('estimates' = results, 'post_treat' = attr(obj[[i]], 'post_treat'))
    }

    res_tab$results <- results_list

    ##
    ## save selection
    ##
    
    ## equivalence result 
    selection <- matrix(NA, nrow = length(attr(obj, "selection")$test_theta), ncol = 2)
    selection[,1] <- attr(obj, "selection")$test_theta
    selection[,2] <- attr(obj, "selection")$test_se
    colnames(selection) <- c("Theta", 'SE')
    rownames(selection) <- paste("M", length(attr(obj, "selection")$test_theta):1, sep = "")
    res_tab$selection <- list("selection" = selection, "model" = attr(obj, "selection")$min_model, 
                             'sign' = attr(obj, 'sign')$res, 'pval' = attr(obj, 'sign')$pvalue,
                             'equivalence' = attr(obj, 'equivalence'), 'bias' = attr(obj, 'sign')$bias
                           )
  }

  return(res_tab)
}
