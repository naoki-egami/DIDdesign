#' Prepare data for the DID design with panel data
#'
#' @param var_outcome A variable name of the outcome.
#' @param var_treat A variable name of the time-varying treatment.
#' @param var_covars A vector of variables names used as covariates.
#'   This should be left as \code{NULL} is no covariates are specified in formula.
#' @param id_time A variable name of the unit index.
#' @param id_time A variable name of the time index.
#' @param data A data frame.
#' @importFrom dplyr %>% all_of pull group_by select mutate rename summarise ungroup left_join
#' @importFrom rlang !! sym .data
#' @keywords internal
did_panel_data <- function(
  var_outcome, var_treat, var_covars,
  var_cluster, id_unit, id_time, data
) {

  ## time and unit
  # var_unit <- pull(data, !!sym(id_unit))
  # var_year <- pull(data, !!sym(id_time))

  ## create a working dataset
  var_select <- c(var_outcome, var_treat, id_unit, id_time)
  if (!is.null(var_covars)) var_select <- c(var_select, var_covars)
  if (!is.null(var_cluster)) var_select <- c(var_select, var_cluster)

  dat_use <- data %>%
    select(all_of(var_select)) %>%
    rename(outcome = !!sym(var_outcome), treatment = !!sym(var_treat),
           id_unit = !!sym(id_unit), id_time = !!sym(id_time)) %>%
    mutate(id_time = as.numeric(as.factor(.data$id_time)))


  ## treatment info
  treat_info <- dat_use %>% group_by(.data$treatment) %>%
    summarise(min_year = min(.data$id_time),
              max_year = max(.data$id_time))

  ## treat time
  treat_year <- treat_info$min_year[2]
  identical(treat_info$max_year[1], treat_info$max_year[2])

  dat_use <- dat_use %>% group_by(.data$id_unit) %>%
    mutate(min_treat = max(.data$treatment)) %>%
    ungroup()


  dat_use <- dat_use %>%
    rename(Gi = .data$min_treat) %>%
    mutate(It = ifelse(.data$id_time >= treat_year, 1, 0)) %>%
    mutate(id_time_std = .data$id_time - treat_year)

  ## compute ΔY
  lag_y <- dat_use %>%
    group_by(.data$Gi, .data$id_time_std) %>%
    summarise(Ymean = mean(.data$outcome), .groups = 'drop') %>%
    mutate(id_time_std = .data$id_time_std + 1)

  ## data
  dat_use <- dat_use %>% left_join(lag_y, by = c("Gi", "id_time_std")) %>%
      mutate(outcome_delta = .data$outcome - .data$Ymean)

  return(dat_use)
}


#' Prepare data for the DID design with repeated cross-section data
#' @importFrom dplyr %>% pull select all_of rename mutate group_by summarise left_join
#' @importFrom rlang !! sym
#' @keywords internal
did_rcs_data <- function(
  var_outcome, var_treat, var_post,
  var_covars = NULL, var_cluster = NULL, id_time, data
) {


  var_year <- pull(data, !!sym(id_time))

  var_select <- c(var_outcome, var_treat, var_post)
  if (!is.null(var_covars)) var_select <- c(var_select, var_covars)
  if (!is.null(var_cluster)) var_select <- c(var_select, var_cluster)

  dat_use <- data %>%
    select(all_of(var_select)) %>%
    rename(Gi = !!sym(var_treat), It = !!sym(var_post),
           outcome = !!sym(var_outcome)) %>%
    mutate(id_time = as.numeric(as.factor(as.character(var_year))))


  ## treatment info
  treat_info <- dat_use %>% group_by(.data$It) %>%
    summarise(min_year = min(.data$id_time),
              max_year = max(.data$id_time))

  ## treat time
  treat_year <- treat_info$min_year[2]
  dat_use <- dat_use %>%
    mutate(id_time_std = .data$id_time - treat_year)

  ## compute ΔY
  lag_y <- dat_use %>%
    group_by(.data$Gi, .data$id_time_std) %>%
    summarise(Ymean = mean(.data$outcome), .groups = 'drop') %>%
    mutate(id_time_std = .data$id_time_std + 1)

  ## data
  dat_use <- left_join(dat_use, lag_y, by = c("Gi", "id_time_std")) %>%
      mutate(outcome_delta = .data$outcome - .data$Ymean)

  return(dat_use)
}
