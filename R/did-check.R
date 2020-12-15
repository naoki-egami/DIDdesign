
#' Assessing Assumptions for Difference-in-differences
#' @inheritParams did
#' @export
#' @importFrom dplyr mutate select %>%
did_check <- function(
  formula, data, id_unit, id_time, design = "did",
  is_panel = TRUE, option = list()
) {

  ## set option
  option <- set_option(option)

  ## -----------------------------------
  ## implement DID estimator
  ## -----------------------------------
  if (design == "did") {
    ## standard design
    fit <- did_check_std(formula, data, id_unit, id_time, is_panel, option)
  } else if (design == "sa"){
    ## staggered adoption design
    if (isFALSE(is_panel)) stop("Only panel data is supported in the SA design.")
    # fit <- did_sad(formula, data, id_unit, id_time, option)
  }

  ## -----------------------------------
  ## equivalence CI
  ## -----------------------------------
  estimate <- fit$est %>%
  mutate(
    CI90_UB_ab = abs(.data$estimate + qnorm(0.95) * .data$std.error),
    CI90_LB_ab = abs(.data$estimate - qnorm(0.95) * .data$std.error)
  ) %>%
  mutate(
    EqCI95_LB = -pmax(.data$CI90_UB_ab, .data$CI90_LB_ab),
    EqCI95_UB = pmax(.data$CI90_UB_ab, .data$CI90_LB_ab)
  ) %>%
  select(-.data$CI90_UB_ab, -.data$CI90_LB_ab)

  out <- list(estimate = estimate, plot = fit$plot$plot, dat_plot = fit$plot$dat_plot)
  class(out) <- c(class(out), "DIDdesign_check")
  return(out)
}

did_check_std <- function(
  formula, data, id_unit, id_time, is_panel = TRUE, option
) {
  ## input check
  if (isFALSE(is_panel)) id_unit <- NULL
  if (isTRUE(is_panel) && is.null(id_unit)) {
    stop("A vaiable name should be provided to id_unit.")
  }

  ## prepare formulas
  fm_prep <- did_formula(formula, is_panel)

  ##
  ## handle cluster variable
  ##
  var_cluster <- option$id_cluster
  if (is.null(var_cluster) && isTRUE(is_panel)) {
    var_cluster <- "id_unit"
    var_cluster_pre <- id_unit
  }

  ## --------------------------------------------
  ## transform data
  ## --------------------------------------------
  if (isTRUE(is_panel)) {
    dat_did <- did_panel_data(
      fm_prep$var_outcome, fm_prep$var_treat, fm_prep$var_covars,
      var_cluster_pre, id_unit, id_time, data
    )
  } else {
    dat_did <- did_rcs_data(
      fm_prep$var_outcome, fm_prep$var_treat, fm_prep$var_post,
      fm_prep$var_covars, var_cluster,id_time, data
    )
  }

  ## --------------------------------------------
  ## estimate placebo test statistics
  ## --------------------------------------------
  did_placebo_est <- did_std_placebo(fm_prep$fm_did[[1]], dat_did, option$lag)

  ## --------------------------------------------
  ## compute std.error via bootstrap
  ## --------------------------------------------
  if (is.null(var_cluster)) {
    id_cluster_vec <- 1:nrow(dat_did)
  } else {
    id_cluster_vec <- unique(pull(dat_did, !!sym(var_cluster)))
  }

  ## setup worker
  if (isTRUE(option$parallel)) {
    plan(multicore)
  } else {
    plan(sequential)
  }

  ## use future_lapply to implement the bootstrap parallel
  est_boot <- future_lapply(1:option$n_boot, function(i) {
    tryCatch({
      did_std_placebo_boot(fm_prep, dat_did, id_cluster_vec, var_cluster, is_panel, option$lag)
    }, error = function(e) {
      NULL
    })
  }, future.seed = TRUE)
  est_boot <- est_boot[lengths(est_boot) != 0]

  ## --------------------------------------------
  ## summarize results
  ## --------------------------------------------
  est_boot <- do.call(rbind, est_boot)
  estimates <- vector("list", length = length(option$lag))
  for (i in 1:length(option$lag)) {
    estimates[[i]] <- data.frame(
      estimate  = did_placebo_est[i],
      lag       = option$lag[i],
      std.error = sd(est_boot[,i])
    )
  }

  ## --------------------------------------------
  ## generate a DID plot
  ## --------------------------------------------
  gg <- did_std_plot(dat_did)

  return(list(est = as_tibble(bind_rows(estimates)), plot = gg))
}


#' Create a did plot for standard design
#' @keywords internal
#' @param data A data object from \code{did_panel_data} or \code{did_rcs_data}
#' @return A ggplot object
#' @importFrom ggplot2 ggplot geom_line geom_point aes geom_vline labs theme_bw
#' @importFrom dplyr %>% across group_by summarise mutate ungroup select
did_std_plot <- function(data) {
  dat_plot <- data %>% group_by(.data$id_time_std, .data$Gi) %>%
         summarise(across(.data$outcome, list(mean = mean, sd = sd))) %>%
         mutate(group = ifelse(.data$Gi == 1, "Treated", "Control")) %>%
         select(group, time_to_treat = id_time_std, outcome_mean, std.error = outcome_sd) %>%
         mutate(CI90_UB = outcome_mean + qnorm(0.95) * std.error,
                CI90_LB = outcome_mean - qnorm(0.95) * std.error) %>%
         ungroup()
  gg <- ggplot(dat_plot, aes(x = time_to_treat, y = outcome_mean, color = group)) +
          geom_vline(xintercept = 0) +  geom_line() + geom_point() +
          labs(x = "Time relative to treatment assignment", ylab = "Mean of the outcome") +
          theme_bw()
  return(list(plot = gg, dat_plot = dat_plot))
}

#' Run a placebo regression
#' @keywords internal
#' @param formula A formula generated by \code{did_formula} function.
#' @param data An output from \code{did_panel_data} or \code{did_rcs_data} function.
#' @param A vector of non-negative lag parameters.
did_std_placebo <- function(formula, data, lags) {
  ## remove all infeasible lag values
  lags <- abs(lags)
  max_lag <- abs(min(data$id_time_std))
  lags <- lags[lags < max_lag]

  ## run placebo regression
  est <- rep(NA, length(lags))
  for (i in 1:length(lags)) {
    time_use <- c(-lags[i], -lags[i]-1)
    dat_use <- data %>% mutate(It = ifelse(.data$id_time_std >= -lags[i], 1, 0))
    fit <- lm(formula, data = filter(dat_use, .data$id_time_std %in% time_use))
    est[i] <- fit$coef['Gi:It']
  }
  return(est)
}


did_std_placebo_boot <- function(
  fm_prep, dat_did, id_cluster_vec, var_cluster, is_panel, lag
) {
  ## sample index
  id_boot <- sample(id_cluster_vec,
    size = length(id_cluster_vec), replace = TRUE
  )

  ## create dataset
  dat_tmp <- list()
  for (j in 1:length(id_boot)) {
    if (is.null(var_cluster)) {
      id_tmp <- id_boot[j]
    } else {
      id_tmp <- which(dat_did[,var_cluster] == id_boot[j])
    }
    dat_tmp[[j]] <- dat_did[id_tmp, ]
    dat_tmp[[j]]$id_unit <- j
  }

  ## create did_data object
  if (isTRUE(is_panel)) {
    dat_boot <- did_panel_data(
      var_outcome = "outcome", var_treat = 'treatment', fm_prep$var_covars,
      var_cluster, id_unit = "id_unit", id_time = 'id_time',
      data = do.call(rbind, dat_tmp)
    )
  } else {
    dat_boot <- did_rcs_data(
      var_outcome = "outcome", var_treat = "Gi", var_post = "It",
      fm_prep$var_covars, var_cluster,
      id_time = "id_time", do.call(rbind, dat_tmp)
    )
  }

  ## fit DID and sDID
  est <- did_std_placebo(fm_prep$fm_did[[1]], dat_boot, lag)
  return(est)
}
