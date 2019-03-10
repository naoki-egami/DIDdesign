#' Nonparametric Double Difference-in-Differences Estimator
#'
#' @param data did_double.data object. See \code{\link{did_double_data}}.
#' @param se_boot A boolean argument.
#'  If set to \code{TRUE}, standard errors are computed by the block bootstrap.
#'  Bootstrap iterations are set by \code{n_boot} argument.
#' @param n_boot The number of bootstrap iterations. Required when \code{se_boot == TRUE}.
#' @param boot_min If \code{TRUE}, bootstrap is used only for the selected model.
#'  This option helps reduce computational burdens.
#' @param select The criteria used to select the best model. The selected model is used to estimate bootstrap variance when \code{boot_min = TRUE}.
#'  Options are "GMM", "tt1" (T-test) and "tt2" (T-test with Bonferroni correction).
#' @param est_did If \code{TRUE}, standard difference-in-differences estimates are computed. Default is \code{TRUE}.
#' @examples
#'\donttest{
#'  # load package
#'  require(DIDdesign)
#'  # load data
#'  data(anzia2012)
#'
#'  # convert data into did double format
#'  out <- did_data(
#'    outcome = anzia2012$lnavgsalary_cpi,
#'    treatment = anzia2012$oncycle,
#'    post_treatment = c(2007, 2008, 2009),
#'    id_subject = anzia2012$district,
#'    id_time = anzia2012$year,
#'    long = TRUE, Xcov = NULL
#'  )
#'
#'  # fit nonparametric double-did method
#'  fit <- did_nonparametric(out)
#'
#'  # plot results
#'  par(mar = c(4, 2.5, 3.5, 1))
#'  plot(fit, full = TRUE, ylim = c(-0.05, 0.05))
#'
#'  # show result summary
#'  summary(fit)
#' }
#' @export
did_nonparametric <- function(
  data, se_boot = FALSE, n_boot = 1000, boot_min = TRUE,
  alpha = 0.05, select = "GMM",
  est_did = TRUE
) {
  ## input checks
  if (!('diddesign_data' %in% class(data))) {
    stop("diddesign_data object should be provided as data.")
  }

  t_pre <- ncol(data[[1]]$Y) - 1
  if (t_pre <= 1) {
    stop("We reuqire more than two pre-treatment periods.\n")
  }



  # ********************************************************* #
  #                                                           #
  #    model selection based on pre-treatment observations    #
  #                                                           #
  # ********************************************************* #
  result <- list()

  m_vec <- 1:t_pre
  select_tmp <- gmm_selection(Y = data[[1]]$Y, D = data[[1]]$D,
                              mvec = m_vec, t_pre = t_pre, select = select, alpha = alpha,
                              n_boot = n_boot)

  min_model  <- select_tmp$min_model

  attr(result, 'selection')  <- select_tmp[c('test_theta', 'test_se', 'min_model')]

  # ********************************************************* #
  #                                                           #
  #       estimate effect for each post-treatment period      #
  #                                                           #
  # ********************************************************* #
  for (j in 1:length(data)) {
    cat("\n... estimating treatment effect for ", attr(data[[j]], 'post_treat'), " ...\n")


    ## ==== point estimate ==== ##
    ## set m_vec
    tmp <- list()
    for (m in m_vec) {
      tmp[[m]] <- didgmmT(Y = data[[j]]$Y, D = data[[j]]$D, M = m, only_beta = TRUE, ep = 0)
    }

    ## ==== variance calulations ==== ##
    if (isTRUE(se_boot) & isTRUE(boot_min)) {
      ## do bootstrap if necessary for the selected model
      cat("... bootstraping to compute standard errors ...\n")
      tmp_min <- didgmmT.boot(
        Y = data[[j]]$Y, D = data[[j]]$D, M = min_model, n_boot = n_boot)

      se95 <- quantile(tmp_min, prob = c(0.025, 0.975))
      se90 <- quantile(tmp_min, prob = c(0.05, 0.95))
    } else if (isTRUE(se_boot)) {
      ## do bootstrap for all models
      cat("... bootstraping to compute standard errors ...\n")
      tmp_min <- list()
      for (m in m_vec) {
        tmp_est <- didgmmT.boot(Y = data[[j]]$Y, D = data[[j]]$D, M = m, n_boot = n_boot)
        tmp_se95 <- quantile(tmp_est, prob = c(0.025, 0.975))
        tmp_se90 <- quantile(tmp_est, prob = c(0.05, 0.95))
        tmp_min[[m]] <- list('boot_est' = tmp_est, 'ci95' = tmp_se95, 'ci90' = tmp_se90)
      }

      se95 <- tmp_min[[min_model]]$ci95
      se90 <- tmp_min[[min_model]]$ci90

    } else {
      cat("... computing asymptotic variance ...\n")
      tmp_min <- list()
      for (m in m_vec) {
        var_est <- didgmmT.variance(tmp[[m]], Y = data[[j]]$Y, D = data[[j]]$D, M = m)
        tmp_se95 <- c(tmp[[m]]$ATT + qnorm(0.05/2) * sqrt(var_est),
                      tmp[[m]]$ATT + qnorm(1 - 0.05/2) * sqrt(var_est))
        tmp_se90 <- c(tmp[[m]]$ATT + qnorm(0.05) * sqrt(var_est),
                      tmp[[m]]$ATT + qnorm(1 - 0.05) * sqrt(var_est))
        tmp_min[[m]] <- list('boot_est' = var_est, 'ci95' = tmp_se95, 'ci90' = tmp_se90)
      }

      se95 <- tmp_min[[min_model]]$ci95
      se90 <- tmp_min[[min_model]]$ci90

    }

    ## ==== standard DiD estimate ==== ##
    if (isTRUE(est_did)) {
      cat("... computing the standard DiD estimate ...\n")

      ## point estimate
      did_est <- std_did(Y = data[[j]]$Y, D = data[[j]]$D)

      ## bootstrap
      # if (isTRUE(se_boot)) {
        tmp_est <- std_did_boot(Y = data[[j]]$Y, D = data[[j]]$D, n_boot = n_boot)
        tmp_se95 <- quantile(tmp_est, prob = c(0.025, 0.975))
        tmp_se90 <- quantile(tmp_est, prob = c(0.05, 0.95))
        did_boot_list <- list('boot_est' = tmp_est, 'ci95' = tmp_se95, 'ci90' = tmp_se90)
      # } else {
      #   did_boot_list <- NULL
      # }

      ## save obj
      did_save <- list("ATT" = did_est, 'results_variance' = did_boot_list)
    } else {
      did_save <- NULL
    }

    ## ==== save results ==== ##
    result[[j]] <- list(
      'results_estimates' = tmp,
      'results_variance' = tmp_min,
      'results_standardDiD' = did_save,
      'min_model' = min_model,
      'select' = select,
      'ATT' = tmp[[min_model]]$ATT,
      'ci95' = se95,
      'ci90' = se90
    )

    attr(result[[j]], 'post_treat') <- attr(data[[j]], 'post_treat')
    attr(result[[j]], 'method') <- 'nonparametric'
  }

  class(result) <- "diddesign"
  return(result)
}
