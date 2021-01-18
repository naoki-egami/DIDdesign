



#' Plot output from did function
#' @export
#' @importFrom dplyr %>% as_tibble bind_rows mutate select arrange
#' @importFrom ggplot2 ggplot aes geom_hline geom_errorbar geom_ribbon geom_point xlim scale_x_continuous theme_bw labs
#' @param x An output from \code{did} function.
#' @param check_fit An output from \code{did_check} function.
#' @param band A boolean argument. If \code{TRUE}, confidence intervals are shown with \code{geom_ribbon} function.
#' @param ... Other parameters. Currently not supported.
plot.DIDdesign <- function(x, check_fit = NULL, band = FALSE, ...) {

  if (!is.null(check_fit)) {
    id_use <- ifelse(attr(check_fit, "design") == "sa", 1, 2)
    dat_plot <- bind_rows(
      check_fit$estimate %>%
        mutate(time = -.data$lag) %>%
        select(estimate = .data$estimate_orig, std.error = .data$std.error_orig,
               time = .data$time),
      as_tibble(x$estimate) %>%
        filter(.data$estimator == "SA-Double-DID" | .data$estimator == "Double-DID") %>%
        select(.data$estimate, .data$std.error, time = .data$lead)
    ) %>% arrange(.data$time)

  } else {
    dat_plot <- as_tibble(x$estimate)  %>%
      filter(.data$estimator == "SA-Double-DID" | .data$estimator == "Double-DID") %>%
      select(.data$estimate, .data$std.error, time = .data$lead) %>%
      arrange(.data$time)
  }

  if (isTRUE(band)) {
  gg <- dat_plot  %>%
    mutate(CI90_LB = .data$estimate - qnorm(0.95) * .data$std.error,
           CI90_UB = .data$estimate + qnorm(0.95) * .data$std.error) %>%
    ggplot(aes(x = .data$time, y = .data$estimate)) +
      geom_hline(yintercept = 0, color = 'gray', linetype = 'dashed') +
      geom_ribbon(aes(ymin = .data$CI90_LB, ymax = .data$CI90_UB),
                  fill = 'gray', alpha = 0.5) +
      geom_line() +
      geom_point() +
      scale_x_continuous(breaks = unique(dat_plot$time)) +
      labs(x = "Time", y = "Estimates (90% CI)") +
      theme_bw()
  } else {
    gg <- dat_plot  %>%
      mutate(CI90_LB = .data$estimate - qnorm(0.95) * .data$std.error,
             CI90_UB = .data$estimate + qnorm(0.95) * .data$std.error) %>%
      ggplot(aes(x = .data$time, y = .data$estimate)) +
        geom_hline(yintercept = 0, color = 'gray', linetype = 'dashed') +
        geom_errorbar(aes(ymin = .data$CI90_LB, ymax = .data$CI90_UB), width = 0.05) +
        geom_point() +
        labs(x = "Time", y = "Estimates (90% CI)") +
        theme_bw()

  }

  if (length(unique(dat_plot$time)) == 1) {
    tt <- unique(dat_plot$time)
    gg <- gg + xlim(tt - 1, tt + 1)
  }
  return(gg)
}
