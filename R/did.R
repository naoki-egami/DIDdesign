

#' DiD 
#' @param formula A formula. \code{outcome ~ treatment}. Covariates can be added by \code{outcome ~ treatment | covariates}.
#' @param data \code{\link{did_data}} object
#' @param method Either \code{"parametric"} or \code{"nonparametric"}.
#' @importFrom dplyr %>% pull tbl_df
#' @importFrom Formula as.Formula
#' @importFrom utils getFromNamespace
#' @export 
did <- function(formula, data, id_subject, id_time, 
  post_treatment, method = 'parametric'
) {

  ## import function
  getFormula <- getFromNamespace("formula.Formula", "Formula")  

  ## input checks 
  if (!("tbl_df" %in% class(data))) data <- data %>% tbl_df()
  if (!exists(id_subject, data)) {
    stop("variable specified for id_subject does not exsit in data")
  }
  if (!exists(id_time, data)) {
    stop("variable specified for id_time does not exsit in data")
  }
  if (!is.numeric(post_treatment)) {
    stop("non-numeric variable is supplied to post_treatment")
  }


  # ********************************************************* #
  #                                                           #
  #          get variable infor and transform data            #
  #                                                           #  
  # ********************************************************* #
  
  ## convert formula 
  formula <- as.Formula(formula)
  f1 <- getFormula(formula, rhs = 1)     

  ## extract variable names 
  outcome    <- all.vars(f1)[1]
  treatment  <- all.vars(f1)[2]

  ## check if covariates are specified or not 
  if (length(formula)[2] == 1) {
    ## subset data 
    dat_use <- did_data(
      outcome        = data %>% pull(outcome), 
      treatment      = data %>% pull(treatment),
      id_subject     = data %>% pull(id_subject),
      id_time        = data %>% pull(id_time),
      post_treatment = post_treatment,      
      long           = TRUE,
      Xcov           = NULL
    )
    
  } else {
    ## extract covariates names 
    f2 <- getFormula(formula, lhs = 0, rhs = 2)
    x_colnames <- all.vars(f2)

    ## subset data 
    dat_use <- did_data(
      outcome        = data %>% pull(outcome),
      treatment      = data %>% pull(treatment),
      id_subject     = data %>% pull(id_subject),
      id_time        = data %>% pull(id_time),
      post_treatment = post_treatment,      
      long           = TRUE,
      Xcov           = data %>% select(x_colnames)
    )
    
  }


  
  return(dat_use)
}
