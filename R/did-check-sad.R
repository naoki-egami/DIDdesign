


#' Staggered Adoption Design
#' @importFrom dplyr %>% as_tibble arrange group_by bind_rows
#' @keywords internal
did_check_sad <- function(formula, data, id_subject, id_time, option) {


  ## --------------------------------------
  ## Prepare inputs
  ## --------------------------------------
  dat_panel <- as_tibble(data) %>%
    rename(id_subject = !!sym(id_subject), id_time = !!sym(id_time)) %>%
    mutate(id_time = as.numeric(as.factor(id_time))) %>%
    arrange(id_subject, id_time)

  is_panel <- TRUE
  var_cluster <- option$id_cluster
  if (is.null(var_cluster) && isTRUE(is_panel)) {
    var_cluster <- "id_unit"
    option$var_cluster_pre <- id_subject
  }

  ## extract variable informations
  if (!("formula" %in% class(formula))) {
    formula <- as.formula(formula)
  }

  all_vars  <- all.vars(formula)
  outcome   <- all_vars[1]
  treatment <- all_vars[2]

  ## prepare custom formula
  fm_prep <- did_formula(formula, is_panel)

  ## --------------------------------------
  ## obtain point estimates for placebo
  ## --------------------------------------
  did_placebo_est <- did_sad_placebo(fm_prep, dat_panel, treatment, outcome, option)

  ## --------------------------------------
  ## Compute standard error via bootstrap
  ## --------------------------------------
  if (isTRUE(option$parallel)) {
    plan(multicore)
  } else {
    plan(sequential)
  }

  est_boot <- future_lapply(1:option$n_boot, function(i) {
      dat_boot <- sample_panel(dat_panel)
      did_sad_placebo(fm_prep, dat_boot, treatment, outcome, option)
  }, future.seed = TRUE)

  ## --------------------------------------
  ## Summarize results
  ## --------------------------------------
  est_boot <- do.call(rbind, est_boot)
  estimates <- vector("list", length = length(option$lag))
  for (i in 1:length(option$lag)) {
    estimates[[i]] <- data.frame(
      estimate  = did_placebo_est[i],
      lag       = option$lag[i],
      std.error = sd(est_boot[,i])
    )
  }

  estimates <- as_tibble(bind_rows(estimates))


  ## --------------------------------------
  ## plot
  ## --------------------------------------
  gg <- did_sad_plot(estimates)
  return(list(est = estimates, plot = list(gg)))
}


did_sad_placebo <- function(fm_prep, dat_panel, treatment, outcome, option) {
  ## --------------------------------------
  ## Prepare inputs
  ## --------------------------------------
  ## create Gmat and index for each design
  Gmat <- create_Gmat(dat_panel, treatment = treatment)
  id_time_use <- get_periods(Gmat, option$thres)
  id_subj_use <- get_subjects(Gmat, id_time_use)
  time_weight <- get_time_weight(Gmat, id_time_use)

  ## -------------------------------
  ## Run placebo regression
  ## -------------------------------
  est_did <- list()
  for (i in 1:length(id_time_use)) {
    est_did[[i]] <- rep(NA, length(option$lag))

    ## subset the data
    dat_use <- dat_panel %>%
      filter(.data$id_subject %in% id_subj_use[[i]]) %>%
      filter(.data$id_time <= id_time_use[[i]])

    ## create did data
    dat_did <- did_panel_data(
      fm_prep$var_outcome, fm_prep$var_treat, fm_prep$var_covars,
      option$var_cluster_pre, id_unit = "id_subject", id_time = "id_time", dat_use
    )

    ## fit placebo regression
    tmp <- did_std_placebo(fm_prep$fm_did[[1]], dat_did, option$lag)
    est_did[[i]][seq_along(option$lag) %in% as.numeric(names(tmp))] <- tmp
  }

  est_did <- do.call(rbind, est_did)

  ## -------------------------------
  ## Take weighted time average
  ## -------------------------------
  estimates <- rep(NA, length(option$lag))
  for (i in 1:length(option$lag)) {
    tmp <- est_did[,i]
    w_use <- time_weight[!is.na(tmp)]
    estimates[i] <- sum( tmp * (w_use / sum(w_use)) )
  }
  names(estimates) <- option$lag
  return(estimates)
}


#' Create a did plot for standard design
#' @keywords internal
#' @return A ggplot object
#' @importFrom ggplot2 ggplot geom_hline geom_point aes geom_errorbar labs theme_bw scale_x_continuous xlim
#' @importFrom dplyr %>% across group_by summarise mutate ungroup select
did_sad_plot <- function(data) {
  dat_plot <- data %>%
  mutate(
    CI90_UB_ab = abs(.data$estimate + qnorm(0.95) * .data$std.error),
    CI90_LB_ab = abs(.data$estimate - qnorm(0.95) * .data$std.error),
    time_to_treat = -.data$lag
  ) %>%
  mutate(
    EqCI95_LB = -pmax(.data$CI90_UB_ab, .data$CI90_LB_ab),
    EqCI95_UB = pmax(.data$CI90_UB_ab, .data$CI90_LB_ab)
  ) %>%
  select(-.data$CI90_UB_ab, -.data$CI90_LB_ab, -.data$lag)

  gg <- ggplot(dat_plot, aes(x = time_to_treat, y = estimate)) +
    geom_hline(yintercept = 0, color = 'gray50', linetype = 'dotted') +
    geom_errorbar(aes(ymin = EqCI95_LB, ymax = EqCI95_UB), width = 0.05, color = '#1E88A8') +
    geom_point(color = '#1E88A8') +
    theme_bw() +
    labs(x = "Time relative to treatment assignment", y = "Test Statistic (95% Equivalence CI)")

  if (length(unique(dat_plot$time_to_treat)) == 1) {
    tt <- abs(unique(dat_plot$time_to_treat))
    gg <- gg+ xlim(-tt-1, -tt + 1)
    #  + scale_x_continuous(breaks = c(-tt - 1, tt, tt + 1))
  }

  return(list(plot = gg, dat_plot = dat_plot))
}
