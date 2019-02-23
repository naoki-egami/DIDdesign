



#' Parametric Version of Double DiD
#' @param data data
#' @param se_boot A boolean argument.
#'  If set to \code{TRUE}, standard errors are computed by the block bootstrap.
#'  Bootstrap iterations are set by \code{n_boot} argument.
#' @param n_boot The number of bootstrap iterations. Required when \code{se_boot == TRUE}.
#' @param boot_min If \code{TRUE}, bootstrap is used only for the selected model.
#'  This option helps reduce computational burdens.
#' @param select The criteria used to select the best model. The selected model is used to estimate bootstrap variance when \code{boot_min = TRUE}.
#'  Options are "HQIC", "BIC", "tt1" (T-test) and "tt2" (T-test with Bonferroni correction).
#' @importFrom plm plm 
#' @examples
#'
#' # load packages
#' require(didrobust)
#' require(dplyr)
#'
#' # load data
#' data(anzia2012)
#'
#' # prepare for input data
#' dat <- did_double_data(
#'   outcome = anzia2012$lnavgsalary_cpi,
#'   treatment = anzia2012$oncycle,
#'   post_treatment = c(2007, 2008, 2009),
#'   id_subject = anzia2012$district,
#'   id_time = anzia2012$year,
#'   long = TRUE,
#'   Xcov = anzia2012 %>% select(lnenrollment, ami_pc, asian_pc, black_pc)
#' )
#'
#' # fit the model
#' fit <- double_did_parametric(dat)
#'
did_parametric <- function(data, se_boot = FALSE, n_boot = 1000, boot_min = TRUE, select = "HQIC"
) {

  ## input checks
  if (!('diddesign_data' %in% class(data))) {
    stop("diddesign_data class object should be provided as data.")
  }
  t_pre <- ncol(data[[1]]$Y) - 1
  if (t_pre <= 1) {
    stop("We reuqire more than two pre-treatment periods.\n")
  }


  n_post <- length(data)
  out <- list()
  for (tt in 1:n_post) {

    # ********************************************************* #
    #                                                           #
    #         model selection by GMM J-stats or t-test          #
    #         - use only pre-treatment data                     #
    #                                                           #
    # ********************************************************* #
    m_vec <- 1:t_pre
    select_tmp <- gmm_selection(Y = data[[tt]]$Y, D = data[[tt]]$D,
                                mvec = m_vec, t_pre = t_pre, select = select, n_boot = n_boot)
    HQIC       <- select_tmp$HQIC
    BIC        <- select_tmp$BIC
    min_model  <- select_tmp$min_model


    # ********************************************************* #
    #                                                           #
    #       run fixed effect model: get demeaned data           #
    #                                                           #
    # ********************************************************* #
    dat_use  <- data[[tt]]$pdata
    fm_list  <- data[[tt]]$formula

    ## fit individual twoway fixed effect model
    ## get demeand matrix and response 
    fit <- dat_trans <- list()
    for (ff in 1:length(fm_list)) {
      fit[[ff]]       <- plm(fm_list[[ff]], data = dat_use, model = 'within', effect = 'twoways')
      dat_trans[[ff]] <- getX(fm_list[[ff]], fit[[ff]], dat_use)
    }

    # ********************************************************* #
    #                                                           #
    #                          GMM                              #
    #                                                           #
    # ********************************************************* #    





    if (isTRUE(se_boot)) {
      boot_out <- list()
      for (ff in 1:length(fm_list)) {
        boot_out[[ff]] <- plm.boot(fm_list[[ff]], data = dat_use, n_boot)
      }
    }
    out[[tt]] <- fit
  }

  return(out)
}


##
##
##  helper functions
##
##

#' Get demeaned X and y 
#'@importFrom plm model.matrix.pFormula pmodel.response
getX <- function(fm, fit, dat) {
  
  ## get y_bar
  y_bar  <- rep(NA, nrow(dat))
  outvar <- all.vars(fm)[1]
  is_na  <- is.na(dat[,outvar])
  y_bar[!is_na] <- pmodel.response(fit)
  names(y_bar)  <- outvar
  
  ## get X_bar
  X_bar <- model.matrix(fit)
  
  return(list("y" = y_bar, "X" = X_bar))
}


#' Estimate GMM 
#'
didgmmT_parametric <- function(dat) {
  k <- length(dat)            ## number of y types 
  p <- ncol(dat[[1]]$X)       ## number of covariates 
  
  
  
  par <- matrix(0, nrow = k, ncol = p)
  
  
  
  for (j in 1:k) {
    t(dat[[j]]$X) %*% (dat[[j]]$y - t(dat[[j]]$X) %*% par[,k])
  }
}



plm.boot <- function(fm, data, n_boot) {
  ##
  ## block bootstrap
  ##

  id_unique <- unique(data$id_subject)
  n_subject <- length(id_unique)
  fit <- list()
  for (bb in 1:n_boot) {
    ## sample id
    id_use <- sample(id_unique, size = n_subject, replace = TRUE)

    ## resample data
    dat_use <- data[data$id_subject %in% id_use, ]

    fit[[bb]] <- plm(fm, data = dat_use, model = 'within', effect = 'twoways')$coef
  }

  return(do.call('rbind', fit))

}
