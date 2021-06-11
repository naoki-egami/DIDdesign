


#' Staggered Adoption Design
#' @inheritParams did_sad
#' @return A list of placebo estimates and plots.
#' @importFrom dplyr %>% as_tibble arrange group_by bind_rows
#' @importFrom stats as.formula sd
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
  setup_parallel(option$parallel)

  est_boot <- future_lapply(1:option$n_boot, function(i) {
      dat_boot <- sample_panel(dat_panel)
      did_sad_placebo(fm_prep, dat_boot, treatment, outcome, option)
  }, future.seed = TRUE)

  ## --------------------------------------
  ## Summarize results
  ## --------------------------------------
  est_boot_std <- do.call(rbind, purrr::map(est_boot, ~.x[,1]))
  est_boot <- do.call(rbind, purrr::map(est_boot, ~.x[,2]))

  estimates <- vector("list", length = length(option$lag))
  for (i in 1:length(option$lag)) {
    estimates[[i]] <- data.frame(
      estimate      = did_placebo_est[i,1],
      lag           = option$lag[i],
      std.error     = sd(est_boot_std[,i]),
      estimate_orig = did_placebo_est[i,2],
      std.error_orig= sd(est_boot[,i])
    )
  }

  estimates <- as_tibble(bind_rows(estimates))


  ## --------------------------------------
  ## plot
  ## --------------------------------------
  p1 <- did_sad_plot(estimates)
  p2 <- did_sad_pattern(dat_panel, treatment, attr(did_placebo_est, "Gmat"))
  return(list(est = estimates, plot = list(p1, p2)))
}


#' SA Placebo Regression
#' @inheritParams sa_double_did
#' @return A matrix of placebo estimates.
#' @keywords internal
#' @importFrom purrr map
did_sad_placebo <- function(formula, dat_panel, treatment, outcome, option) {
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
  est_did <- est_did_std <- list()
  for (i in 1:length(id_time_use)) {
    est_did[[i]] <- est_did_std[[i]] <- rep(NA, length(option$lag))

    ## subset the data
    dat_use <- dat_panel %>%
      filter(.data$id_subject %in% id_subj_use[[i]]) %>%
      filter(.data$id_time <= id_time_use[[i]])

    ## create did data
    dat_did <- did_panel_data(
      formula$var_outcome, formula$var_treat, formula$var_covars,
      option$var_cluster_pre, id_unit = "id_subject", id_time = "id_time", dat_use
    )

    ## fit placebo regression
    tmp <- did_std_placebo(formula$fm_did[[1]], dat_did, option$lag)
    est_did[[i]][seq_along(option$lag) %in% as.numeric(names(tmp$est))] <- tmp$est
    est_did_std[[i]][seq_along(option$lag) %in% as.numeric(names(tmp$est))] <- tmp$est_std
  }

  est_did <- do.call(rbind, est_did)
  est_did_std <- do.call(rbind, est_did_std)

  ## -------------------------------
  ## Take weighted time average
  ## -------------------------------
  estimates <- matrix(NA, nrow = length(option$lag), ncol = 2)
  for (i in 1:length(option$lag)) {
    tmp <- est_did[,i]
    tmp_std <- est_did_std[,i]
    w_use <- time_weight[!is.na(tmp)]
    estimates[i,1] <- sum( tmp_std * (w_use / sum(w_use)) )
    estimates[i,2] <- sum( tmp * (w_use / sum(w_use)) )
  }
  rownames(estimates) <- option$lag
  attr(estimates, "Gmat") <- Gmat
  return(estimates)
}


#' Create a did plot for standard design
#' @keywords internal
#' @param data A data frame of double did estimates.
#' @return A ggplot2 object
#' @importFrom ggplot2 ggplot geom_hline geom_point aes geom_errorbar labs theme_bw scale_x_continuous xlim
#' @importFrom dplyr %>% across group_by summarise mutate ungroup select
#' @importFrom stats qnorm
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
  select(.data$estimate, .data$time_to_treat, .data$std.error,
         .data$EqCI95_LB, .data$EqCI95_UB)

  gg <- ggplot(dat_plot, aes(x = .data$time_to_treat, y = .data$estimate)) +
    geom_hline(yintercept = 0, color = 'gray50', linetype = 'dotted') +
    geom_errorbar(aes(ymin = .data$EqCI95_LB, ymax = .data$EqCI95_UB),
                  width = 0.05, color = '#1E88A8') +
    theme_bw() +
    labs(x = "Time relative to treatment assignment", y = "95% Standardized Equivalence CI")

  if (length(unique(dat_plot$time_to_treat)) == 1) {
    tt <- abs(unique(dat_plot$time_to_treat))
    gg <- gg+ xlim(-tt-1, -tt + 1)
    #  + scale_x_continuous(breaks = c(-tt - 1, tt, tt + 1))
  }

  return(list(plot = gg, dat_plot = dat_plot))
}


#' Generate a pattern plot
#' @param data Panel data.
#' @param treatment Name of the treatment variable.
#' @param Gmat A matrix of treatment patterns.
#' @return A list of the treatment variation plot and the data used to generate the plot.
#' @importFrom ggplot2 ggplot aes geom_tile scale_fill_manual theme_bw labs theme element_blank element_text
#' @importFrom dplyr %>% mutate
#' @importFrom rlang !! sym
#' @keywords internal
did_sad_pattern <- function(data, treatment, Gmat) {
  treat_timing <- order(apply(Gmat, 1, function(x) ifelse(sum(x == 1) >= 1, which(x == 1), Inf)),
                        decreasing = TRUE)
  data <- data %>%
    mutate(id_subject = factor(.data$id_subject, levels = treat_timing),
           treatment  = ifelse(!!sym(treatment) == 1, "treated", "control"))
  gg <- ggplot(data, aes(x = .data$id_time, y = .data$id_subject, fill = !!sym(treatment))) +
    geom_tile() +
    scale_fill_manual(values = c("lightgray", '#1E88A8')) +
    theme_bw() +
    theme(axis.text.x = element_blank(), axis.text.y = element_blank(),
           axis.title.x = element_text(vjust = -4),  axis.ticks.x = element_blank()) +
    labs(x = "Time", y = "Unit", fill = "Status")

  return(list(plot = gg, dat_plot = data))
}
