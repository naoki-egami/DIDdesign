

#' Function to implement the Staggered Adption Design
#'
#' @import panelr
sa_did <- function(formula, data, id_subject = NULL, id_time) {

  ## keep track of long-form data with panel class from \code{panelr} package
  dat_panel <- panel_data(data, id = id_subject, wave = id_time)

  ##
  all_vars  <- all.vars(formula)
  outcome   <- all_vars[1]
  treatment <- all_vars[2]

  ##

  
}





#' Create a G matrix
#'
#' @param dat_panel A class of \code{panelr} object.
#' @import panelr
#' @importFrom rlang sym !!
create_Gmat <- function(dat_panel, treatment) {
  ## tranform data into wide
  dat_wide <- dat_panel %>% select(!!sym(treatment)) %>%
    widen_panel()

  ##
}




sd_did_parametric <- function() {

}
