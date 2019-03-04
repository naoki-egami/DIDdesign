#' Model selection by GMM
#' @param Y outcome
#' @param D treatment
#' @param m_vec moment sets
#' @param select selection methods
#' @return a list that consists of
#' \itemize{
#'   \item{min_model}{a index of selected moment.}
#'   \item{HQIC}{HQIC values of all models considered.}
#'   \item{BIC}{BIC values of all models considered.}
#'   \item{test_theta}{test statistics of all models. What value is returned depends on selectin criteria.}
#' }
#' @keywords internal
gmm_selection <- function(Y, D, mvec, t_pre, select, alpha, n_boot) {
  select_list <- list()
  T0 <- ncol(Y) - 1
  Y0 <- Y[, 1:T0];

  ## GMM based selection
  for (m in 1:(length(mvec)-1)) {
    select_list[[m]] <- didgmm_test2(Y = Y0, D0 = D, M = mvec[m])
  }

  HQIC <- rep(0, length(mvec)) # sapply(select_list, function(x) x$HQIC)
  BIC  <- rep(0, length(mvec)) # sapply(select_list, function(x) x$BIC)
  zeta <- sapply(select_list, function(x) x$zeta)

  if (select %in% c("GMM")) {
    ## model selection by BIC and HQIC
    # min_model  <- ifelse(select == "HQIC", which.min(HQIC), which.min(BIC))
    test_theta <- zeta
    test_se    <- sapply(select_list, function(x) x$se)
    min_model  <- check_cover(test_theta, test_se, alpha, mvec)

    names(test_theta) <- names(test_se) <- paste("M", mvec[-length(mvec)], sep = '')
  } else if (select %in% c('tt1', 'tt2')) {
    ## model selection by t-test
    m_vec <- rev(1:(t_pre - 1))
    min_model <- t_pre
    count <- 0
    for (m in 1:(length(m_vec)-1)) {
      m_use <- mvec[m]
      select_list[[m]] <- trendT_test(Y = Y, D = D, M = m_use, boot = TRUE, n_boot = n_boot)
      use_reject <- ifelse(select == "tt1", select_list[[m]]$boot$reject, select_list[[m]]$boot$reject_correction)
      if (use_reject) {
        min_model <- mvec[m+1]
        break;
      } else {
        count <- count + 1
      }
    }
    test_theta <- sapply(select_list, function(x) x$delta)
    test_se    <- sapply(select_list, function(x) sd(x$boot))
    if (count > 0) min_model <- m_vec[count]

  } else {
    stop("invalid input for select")
  }

  out_list <- list('min_model' = min_model, 'HQIC' = HQIC, 'BIC' = BIC,
                   'test_theta' = test_theta, 'test_se' = test_se)
  attr(out_list, 'method') <- 'gmm_selection'
  return(out_list)
}

#' sequentially check coverage of zero
#' @keywords internal
check_cover <- function(theta, se, alpha, mvec) {
  min_model <- max(mvec)
  mvec_rev  <- rev(mvec)[-1]
  theta_rev <- rev(theta)
  se_rev    <- rev(se)
  for (i in 1:length(mvec_rev)) {
    is_cover <- (theta_rev[i] + qnorm(alpha/2) * se_rev[i] <= 0) &
                (theta_rev[i] + qnorm(1 - alpha/2) * se_rev[i] >= 0)
    if (is_cover) {
      min_model <- mvec_rev[i]
    } else {
      break;
    }
  }

  return(min_model)
}

#' Model selectin based on linear FE model
#' @param data panel data from \code{diddesign_data} object.
#' @param fm_list A list of formulae from \code{diddesign_data} object.
#' @param post_time Time index for the post Treatment period.
#' @param alpha Level of rejection. Default is at the 5\% level
#' @importFrom plm plm vcovHC
#' @importFrom utils getFromNamespace
#' @keywords internal
fe_selection <- function(data, fm_list, post_time, alpha = 0.05) {

  # load function
  diff.pseries <- getFromNamespace("diff.pseries", "plm")

  # exclude post-treatment periods
  new_treatment      <- diff.pseries(data$treatment)
  data_use           <- data[data$id_time != post_time, ]
  data_use$treatment <- na.omit(new_treatment)

  # run fe
  min_model  <- length(fm_list)
  test_theta <- test_se <- rep(NA, (length(fm_list)-1))

  # upward model selection (from highest order to lowest order)
  for (ff in (length(fm_list)-1):1) {
    ## estimate "placebo" model
    tmp <- plm(fm_list[[ff]], data = data_use, method = 'within', effect = 'twoways')
    test_theta[ff] <- tmp$coef[1]
    test_se[ff]    <- sqrt(vcovHC(tmp, cluster = 'group', type = 'HC2')[1,1])
    cover_zero     <- ((test_theta[ff] + qnorm(alpha/2) * test_se[ff]) <= 0) &
                      ((test_theta[ff] + qnorm(1 - alpha/2) * test_se[ff]) >= 0)
    if (cover_zero) { min_model <- ff }
  }

  out_list <- list("min_model" = min_model, 'BIC' = NULL, "HQIC" = NULL,
              'test_theta' = test_theta, 'test_se' = test_se)
  attr(out_list, 'method') <- 'fe_selection'
  return(out_list)
}



#' Model selectin based on linear FE model
#' @importFrom dplyr %>% filter
#' @keywords internal
rcs_selection <- function(data, fm_list, post_time, alpha = 0.05) {

  ## shift post indicator
  data_pre <- data %>% filter(id_time < post_time) %>%
    mutate(post = ifelse(id_time == max(id_time), 1, 0))
  min_model  <- length(fm_list)
  test_theta <- test_se <- rep(NA, (length(fm_list)-1))

  for (ff in (length(fm_list)-1):1) {
    tmp <- lm(fm_list[[ff]], data = data_pre)
    test_theta[ff] <- tmp$coef['treatment:post']
    test_se[ff]    <- summary(tmp)$coef['treatment:post', 2]

    if ((test_theta[ff] + qnorm(alpha/2) * test_se[ff]) <= 0 &
        (test_theta[ff] + qnorm(1 - alpha/2) * test_se[ff]) >= 0) {
      min_model <- ff
    }
  }

  return(list("min_model" = min_model, 'BIC' = NULL, "HQIC" = NULL, '
              test_theta' = test_theta, 'test_se' = test_se))

}
