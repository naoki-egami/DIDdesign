




#' Double DID Estimator for the Staggered Adoption
#' @keywords internal
ddid_sa_new <- function() {
  
}


#' Double DID Data Format for the Staggered Adoption 
#' @keywords internal
ddid_sa_data <- function(var_treat, id_unit, id_time, data) {
  
  ## rename variables 
  dat_use <- data %>% 
    rename(id_unit = !!sym(id_unit), id_time = !!sym(id_time), 
           treatment = !!sym(var_treat)) %>% 
    mutate(id_time = as.numeric(as.factor(as.character(id_time))))


  ## create treatment group 
  treat_summary <- dat_use %>% group_by(id_unit) %>% 
    summarise(treat_prop = mean(treatment)) %>% 
    mutate(treat_group = as.numeric(as.factor(treat_prop)))
      
  ## drop units that is always treated in the sample 
  unit_drop <- treat_summary %>% filter(treat_prop == 1) %>% 
    pull(id_unit)
  
  
  ## merge summary data to the data 
  dat_use <- left_join(dat_use, treat_summary, by = "id_unit")
  
  ## identify the first year of the treatment 
  treat_year <- dat_use %>% filter(treatment == 1) %>% 
    group_by(treat_group) %>% 
    summarise(min_year = min(id_time))
  
  dat_use <- dat_use %>% left_join(treat_year, by = "treat_group") %>% 
    mutate(min_year = ifelse(is.na(min_year), max(id_time) + 1, min_year)) %>%
    mutate(G = case_when(
      id_time < min_year  ~ 0,
      id_time == min_year ~ 1, 
      id_time > min_year  ~ -1
    ))
  
  return(dat_use)

}
