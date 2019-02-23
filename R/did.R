

#' DiD 
#' @param formula A formula. \code{outcome ~ treatment}. Covariates can be added by \code{outcome ~ treatment | covariates}.
#' @param data \code{\link{did_data}} object
#' @param method Either \code{"parametric"} or \code{"nonparametric"}.
#' @param id_subject subject id.
#' @param id_time time id.
#' @param post_treatment time index for post treatment periods.
#' @param se_boot If \code{TRUE} blockbootstrap is used for computing standard errors.
#' @param n_boot Number of bootstrap iterations. 
#' @param boot_min If \code{TRUE} bootstrap is carried out only for the selected model.
#' @param select Selection criteria.
#' @param
#' @importFrom dplyr %>% pull tbl_df
#' @importFrom Formula as.Formula
#' @importFrom utils getFromNamespace
#' @examples
#' # load package
#' require(DIDdesign)
#' 
#' # load  data 
#' data(anzia2012)
#' 
#' # nonparametric version without covariates
#' fit1 <- did(
#'   formula = "lnavgsalary_cpi ~ oncycle",
#'   data = anzia2012,
#'   id_subject = "district", id_time = "year",
#'   post_treatment = c(2007, 2008, 2009),
#'   method = 'nonparametric',
#'   se_boot = FALSE
#' )
#' 
#' @export 
did <- function(formula, data, id_subject, id_time, post_treatment, 
  method = 'parametric', se_boot = FALSE, n_boot = 1000, boot_min = TRUE, 
  select = 'HQIC'
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
  if (!(method %in% c("parametric", "nonparametric"))) {
    stop("Either `parametric' or `nonparametric' is allowed for method option")
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
    
    is_covariates <- FALSE
    
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
    
    is_covariates <- TRUE  
  }

  # ********************************************************* #
  #                                                           #
  #             estimate treatment effect                     #
  #                                                           #
  # ********************************************************* #    
  
  if(method == "nonparametric" & !isTRUE(is_covariates)) {
    # TO USE NONPARAMETRIC METHOD, BOTH CONDITIONS SHOULD BE MET:
    # 1. method should be 'nonparametric', AND
    # 2. is_covariates should be FALSE (NO COVARIATES ALLOWED)
    fit <- did_nonparametric(
      data = dat_use, 
      se_boot = se_boot, n_boot = n_boot, boot_min = boot_min,
      select  = select, est_did = TRUE
    )
  } else if (method == "nonparametric" & isTRUE(is_covariates)) {
    warning("Nonparametric option is not available with covariates. 
             Parametiric method is used instead\n")
  } else {
    
  }
  
  
  
  
  return(fit)
}
