


#' DID desing
#' @export
#' @param formula A formula of the form \code{y ~ treatment | x1 + x2}.
#' @param data A data frame.
#' @param id_unit A variable name of unit (e.g., country name, respondent id).
#' @param id_time A variable name of time (e.g., year).
#' @param design The design to be used: either \code{"did"} (the standard difference-in-differences design) or \code{"sa"} (the staggered adoption design).
#' The default is \code{"did"}.
#' @param is_panel A boolean argument. This should be \code{TRUE} when the dataset is panel (i.e., the same units are repeately observed over time); This should be \code{FALSE} when the dataset is the repeated cross-section (RCS) where different sets of units are observed at each time point.
#' @param option A list of option parameters.
#' @import Formula
did_new <- function(
  formula, data, id_unit, id_time,
  design = "did", is_panel = TRUE,
  option = list()) {


  ## obtain variable names
  fm <- as.Formula(formula)
  n_vars_fm <- length(fm)
  if (n_vars_fm[1] != 1) stop("Invalid formula (lhs)")
  var_outcome <- all.vars(formula(fm, lhs = 1, rhs = 0))
  var_treat   <- all.vars(formula(fm, lhs = 0, rhs = 1))

  if (n_vars_fm[2] == 2) {
    var_covars <- all.vars(formula(fm, lhs = 0, rhs = 2))
    fm_covar <- formula(fm, lhs = 0, rhs = 2)
  } else if (n_vars_fm[2] == 1){
    var_covars <- NULL
  } else {
    stop("Invalid formula (rhs)")
  }
  
  ## update formula 
  if (is.null(var_covars)) {
    fm_did <- list(
      as.formula(outcome ~ Gi + It + Gi : It),
      as.formula(outcome_delta ~ Gi + It + Gi : It))
  } else {
    fm_did <- list(
      update(fm_covar, outcome ~ Gi + It + Gi : It + .),
      update(fm_covar, outcome_delta ~ Gi + It + Gi : It + .))
  }
  
  ## transform data 
  dat_did <- did_panel_data(
    var_outcome, var_treat, var_covars,
    id_unit, id_time, data
  )
    
  fit_did  <- ddid_fit(fm_did, dat_did, lead = 1)

  ## bootstrap to compute the weighting matrix W 
  id_cluster <- "id_unit"
  n_boot <- 100
  est_boot <- matrix(NA, nrow = n_boot, ncol = length(fm_did))
  for (i in 1:n_boot) {
    ## sample index 
    id_cluster_vec <- pull(dat_did, !!sym(id_cluster)) %>% unique()
    id_boot <- sample(id_cluster_vec, 
      size = length(id_cluster_vec), replace = TRUE)
    
    dat_tmp <- list()
    for (j in 1:length(id_boot)) {
      id_tmp <- which(dat_did[,id_cluster] == id_boot[j])
      dat_tmp[[j]] <- dat_did[id_tmp, ]
      dat_tmp[[j]]$id_unit <- j
    }


    dat_boot <- did_panel_data(
      var_outcome = "outcome", var_treat = 'treatment', var_covars,
      id_unit = "id_unit", id_time = 'id_time', do.call(rbind, dat_tmp)
    ) 
    est_boot[i,] <- ddid_fit(fm_did, dat_boot, lead = 1)  
  }
  
  
  ## compute weights 
  W <- cov(est_boot)
  w_did  <- (W[1,1] - W[1,2]) / (W[1,1] + W[2,2] - 2 * W[1,2])
  w_sdid <- (W[2,2] - W[1,2]) / (W[1,1] + W[2,2] - 2 * W[1,2])
  
  ddid <- fit_did[1] * w_did + fit_did[2] * w_sdid
  ddid_var <- w_did^2 * W[1,1] + w_sdid^2 * W[2,2] + 
              2 * w_did * w_sdid * W[1,2]
  return(list(
    ddid = ddid, ddid_var = ddid_var,
    W = W, w_did = w_did, w_sdid = w_sdid,
    fit_did = fit_did
  ))
}


#' Prepare data for the did design with panel data 
did_panel_data <- function(
  var_outcome, var_treat, var_covars,
  id_unit, id_time, data 
) {

  ## time and unit 
  var_unit <- pull(data, !!sym(id_unit))
  var_year <- pull(data, !!sym(id_time))
  
  ## create a working dataset 
  if (is.null(var_covars)) {
    var_select <- c(var_outcome, var_treat)    
  } else {
    var_select <- c(var_outcome, var_treat, var_covars)
  }

  dat_use <- data %>% 
    select(all_of(var_select)) %>% 
    mutate(id_unit = var_unit, 
           id_time = as.numeric(as.factor(as.character(var_year)))) %>% 
    rename(outcome = !!sym(var_outcome), treatment = !!sym(var_treat))


  ## treatment info 
  treat_info <- dat_use %>% group_by(treatment) %>% 
    summarise(min_year = min(id_time), 
              max_year = max(id_time))
    
  ## treat time 
  treat_year <- treat_info$min_year[2]
  identical(treat_info$max_year[1], treat_info$max_year[2])

  dat_use <- dat_use %>% group_by(id_unit) %>% 
    mutate(min_treat = max(treatment)) %>% 
    ungroup()


  dat_use <- dat_use %>% 
    rename(Gi = min_treat) %>% 
    mutate(It = ifelse(id_time >= treat_year, 1, 0)) %>% 
    mutate(id_time_std = id_time - treat_year)

  ## compute ΔY 
  lag_y <- dat_use %>% 
    group_by(Gi, id_time_std) %>% 
    summarise(Ymean = mean(outcome), .groups = 'drop') %>% 
    mutate(id_time_std = id_time_std + 1)

  ## data 
  dat_use <- dat_use %>% left_join(lag_y, by = c("Gi", "id_time_std")) %>% 
      drop_na(Ymean) %>% 
      mutate(outcome_delta = outcome - Ymean)
  
  
  return(dat_use)
}




# @param data An output of \code{did_panel_data()}.
# @param formula A formula of the form \code{y ~ Gi + It + Gi * It + x1 + x2}.
ddid_fit <- function(formula, data, lead = 1) {
  time_use <- c(-1, lead - 1)
  est <- map_dbl(formula, function(fm) {
    fit <- lm(fm, data = filter(data, id_time_std %in% time_use))
    return(fit$coef['Gi:It'])
  })
  
  return(est)
}


#' Set default value of options
#' @keywords internal
set_option <- function(option) {

  option$se_boot <- TRUE
  option$n_boot <- 100
  option$id_cluster <- "id_unit"
  option$order <- 2
  return(option)
}
