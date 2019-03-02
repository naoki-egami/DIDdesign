#' Model selection by GMM
#' @param Y outcome
#' @param D treatment
#' @param m_vec moment sets
#' @param select selection methods
#' @keywords internal
gmm_selection <- function(Y, D, mvec, t_pre, select, n_boot) {
    select_list <- list()
    for (m in mvec) {
      select_list[[m]] <- didgmm_test(Y = Y, D = D, M = m)
    }

    HQIC <- sapply(select_list, function(x) x$HQIC)
    BIC  <- sapply(select_list, function(x) x$BIC)

    if (select %in% c("BIC", "HQIC")) {
      ## model selection by BIC and HQIC
      min_model <- ifelse(select == "HQIC", which.min(HQIC), which.min(BIC))
    } else if (select %in% c('tt1', 'tt2')) {
      ## model selection by t-test
      m_vec <- rev(1:(t_pre - 1))
      min_model <- t_pre
      count <- 0
      for (m in m_vec) {
        select_list <- trendT_test(Y = Y, D = D, M = m, boot = TRUE, n_boot = n_boot)
        use_reject <- ifelse(select == "tt1", select_list$boot$reject, select_list$boot$reject_correction)
        if (use_reject) {
          min_model <- m+1
          break;
        } else {
          count <- count + 1
        }
      }

      if (count > 0) min_model <- m_vec[count]


    } else {
      stop("invalid input for select")
    }

    return(list('min_model' = min_model, 'HQIC' = HQIC, 'BIC' = BIC))
}



#' Model selectin based on linear FE model
#' @importFrom plm plm vcovHC
#' @keywords internal
fe_selection <- function(data, fm_list, post_time, alpha = 0.05) {

  diff.pseries <- getFromNamespace("diff.pseries", "plm")

  # exclude post-treatment periods
  new_treatment      <- diff.pseries(data$treatment)
  data_use           <- data[data$id_time != post_time, ]
  data_use$treatment <- na.omit(new_treatment)

  # run fe
  min_model  <- length(fm_list)
  test_theta <- test_se <- rep(NA, (length(fm_list)-1))
  for (ff in (length(fm_list)-1):1) {
    tmp <- plm(fm_list[[ff]], data = data_use, method = 'within', effect = 'twoways')
    test_theta[ff] <- tmp$coef[1]
    test_se[ff]    <- sqrt(vcovHC(tmp, cluster = 'group', type = 'HC2')[1,1])

    if ((test_theta[ff] + qnorm(alpha/2) * test_se[ff]) <= 0 &
        (test_theta[ff] + qnorm(1 - alpha/2) * test_se[ff]) >= 0) {
      min_model <- ff
    }
  }

  return(list("min_model" = min_model, 'BIC' = NULL, "HQIC" = NULL, 'test_theta' = test_theta,
              'test_se' = test_se))
}



#' Model selectin based on linear FE model
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

  return(list("min_model" = min_model, 'BIC' = NULL, "HQIC" = NULL, 'test_theta' = test_theta,
              'test_se' = test_se))

}
