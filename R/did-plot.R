



#' Plot output from did function
#' @export
#' @importFrom dplyr %>% as_tibble bind_rows mutate select arrange
#' @importFrom ggplot2 ggplot aes geom_hline geom_errorbar geom_point xlim theme_bw labs
plot.DIDdesign <- function(x, check_fit = NULL, ...) {

  if (!is.null(check_fit)) {
    id_use <- ifelse(attr(check_fit, "design") == "sa", 1, 2)
    dat_plot <- bind_rows(
      check_fit$plot[[id_use]]$dat_plot %>%
        mutate(time = time_to_treat) %>%
        select(estimate, std.error, time),
      as_tibble(x$estimate) %>%
        filter(estimator == "SA-Double-DID" | estimator == "Double-DID") %>%
        select(estimate, std.error, time = lead)
    )

  } else {
    dat_plot <- as_tibble(x$estimate)  %>%
      filter(estimator == "SA-Double-DID" | estimator == "Double-DID") %>%
      select(estimate, std.error, time = lead)
  }

  gg <- dat_plot %>% arrange(time) %>%
    mutate(CI90_LB = estimate - qnorm(0.95) * std.error,
           CI90_UB = estimate + qnorm(0.95) * std.error) %>%
    ggplot(aes(x = time, y = estimate)) +
      geom_hline(yintercept = 0, color = 'gray', linetype = 'dashed') +
      geom_errorbar(aes(ymin = CI90_LB, ymax = CI90_UB), width = 0.05) +
      geom_point() +
      labs(x = "Time", y = "Estimates (90% CI)") +
      theme_bw()

  if (length(unique(dat_plot$time)) == 1) {
    tt <- unique(dat_plot$time)
    gg <- gg + xlim(tt - 1, tt + 1)
  }
  return(gg)
}
