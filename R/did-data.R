

#' Transform data into did_double format
#'
#' A function to transform panel data into did double format.
#'
#' @param outcome Outcome vector when \code{long = TRUE} or outcome matrix when \code{long = FALSE}.
#' @param treatment Treatment vector. When \code{long = TRUE}, the length of this vector should be the same as the length of \code{outcome} vector, while the length of \code{treatment} should be the same as the number of rows (i.e. the number of subjects) in \code{outcome} when the input is in the wide format.
#' @param post_treatment A scalar or a vector of time index that indicats the post treatment period(s).
#'  If left as \code{NULL}, the last period in a variable given as \code{id_time} is used as the post treatment period.
#' @param id_subject A vector of subject index. Required if \code{long = TRUE}.
#' @param id_time A vector of time index. Required if \code{long = TRUE}.
#' @param long A boolean input indicating whether the input is in the long format or not.
#' Default is \code{long = TRUE}.
#' @param x_formula a formula of the right hand side (covariates). If left as \code{NULL}, original variables supplied in \code{Xcov} will be used as is.
#' @return A \code{diddesign} class object which is a list.
#'  Each element corresponds to each post-treatment period.
#'  Each element of the returned \code{diddesign} object is also a list consists of the following:
#' \itemize{
#'  \item \code{Y}: outcome vector.
#'  \item \code{D}: treatment vector.
#'  \item \code{formula}: A list of formula for \code{\link{lm}}.
#'  \item \code{pdata}: A data frame appended transformed outcomes (e.g., yd0, yd1, etc).
#'}
#' @examples
#'  # load data from didrobust package
#'  data(anzia2012)
#'
#'  # convert data into did double format
#'  out <- did_data(
#'    outcome = anzia2012$lnavgsalary_cpi,
#'    treatment = anzia2012$oncycle,
#'    post_treatment = c(2007, 2008, 2009),
#'    id_subject = anzia2012$district,
#'    id_time = anzia2012$year,
#'    long = TRUE
#'  )
#'
#'  # view a summary
#'  summary(out)
#' @export
did_data <- function(
  outcome, treatment,
  post_treatment = NULL, id_subject = NULL, id_time = NULL, long = TRUE,
  Xcov = NULL, x_formula = NULL
) {

  is_rcs <- FALSE
  if (is.null(id_subject)) { is_rcs <- TRUE }

  ##
  ## data transformation
  ##
  if (isTRUE(is_rcs)) {

    ## repeated cross-section + long
    out <- did_data_rcs(outcome, treatment, post_treatment, id_time, Xcov, x_formula)
  } else if (!isTRUE(is_rcs) & isTRUE(long)) {

    ## panel + long
    out <- did_data_panelL(outcome, treatment, post_treatment, id_subject, id_time, Xcov, x_formula)
  } else {

    ## panel + wide
    if (!is.null(Xcov)) {
      warning("To use covariates, please transform data into the long format.")
    }
    out <- did_data_panelW(outcome, treatment, post_treatment)
  }

  class(out) <- c("diddesign", "diddesign_data")
  return(out)
}


#' Data processing function for repeated cross-section data
#' @param outcome a vector of outcome variable.
#' @param treatment a vector of treatment indicator variable.
#' @param post_treatment a vector of numeric time index for post treatment periods.
#' @param id_time a vector of time index variable.
#' @param Xcov a matrix of covariates.
#' @return \code{\link{did_data_rcs}} returns a list. Each element corresponds to each post-treatment period.
#'  Each element of the returned list is also a list consists of the following:
#' \itemize{
#'  \item \code{Y}: outcome vector.
#'  \item \code{D}: treatment vector.
#'  \item \code{formula}: A list of formula for \code{\link{lm}}.
#'  \item \code{pdata}: A data frame appended transformed outcomes (e.g., yd0, yd1, etc).
#'}
#' @importFrom zoo zoo
#' @importFrom dplyr %>% filter summarise group_by pull mutate tbl_df
#' @importFrom utils getFromNamespace
#' @keywords internal
did_data_rcs <- function(outcome, treatment, post_treatment, id_time, Xcov, x_formula = NULL) {

  diff.zoo <- getFromNamespace("diff.zoo", "zoo")
  out <- list()

  if (!is.null(Xcov)) {  ## use this when Xcov is supplied
    if (is.null(colnames(Xcov))) {
      colnames(Xcov) <- paste("XR", 1:ncol(Xcov), sep = '')
    }
    if (any(colnames(Xcov) %in% c("", " "))) {
      ## rename colnames if they contain empty name
      rename_idx <- which(colnames(Xcov) %in% c("", " "))
      colnames(Xcov)[rename_idx] <- paste('XR', length(rename_idx), sep = '')
    }
    x_colname <- colnames(Xcov)
  }


  for (tt in 1:length(post_treatment)) {
    ## make data
    if (is.null(Xcov)) {
      dat2 <- data.frame(outcome, treatment, id_time) %>%
        tbl_df() %>%
        dplyr::filter(id_time < post_treatment[1] | id_time == post_treatment[tt]) %>%
        mutate(id_time = ifelse(id_time == post_treatment[tt], post_treatment[1], id_time)) %>%
        na.omit()
    } else {
      dat2 <- data.frame(outcome, treatment, id_time, Xcov) %>%
        tbl_df() %>%
        dplyr::filter(id_time < post_treatment[1] | id_time == post_treatment[tt]) %>%
        mutate(id_time = ifelse(id_time == post_treatment[tt], post_treatment[1], id_time)) %>%
        na.omit()
    }

    ## take diff
    ymean  <- dat2 %>% group_by(treatment, id_time) %>% summarise(ymean = mean(outcome))
    y0mean <- ymean %>% filter(treatment == 0) %>% pull(ymean)
    y1mean <- ymean %>% filter(treatment == 1) %>% pull(ymean)

    pre_time <- length(unique(id_time)) - length(post_treatment)
    fm <- list()

    ## k = 1
    dat2 <- dat2 %>% mutate(yd0 = outcome, post = ifelse(id_time == post_treatment[1], 1, 0))
    id_time_unique <- sort(unique(dat2$id_time))

    if (is.null(Xcov)) {
      fm[[1]] <- as.formula(paste("yd",0, " ~ treatment + post + treatment * post", sep = ''))
    } else if (is.null(x_formula)){
      fm[[1]] <- as.formula(paste("yd",0, " ~ treatment + post + treatment * post + ",
                            paste(x_colname, collapse = "+"), sep = ''))
    } else {
      fm[[1]] <- as.formula(paste("yd",0, " ~ treatment + post + treatment * post + ",
                            paste(attr(terms(x_formula), 'term.labels'), collapse = "+"), sep = ''))
    }




    ## k = 2 ~ T
    for (k in 2:pre_time) {
      ydiff <- rep(NA, nrow(dat2))

      tmp <- rcs_mean_fill(
        outcome = dat2$outcome, treatment = dat2$treatment,
        time_index = dat2$id_time, time_unique = id_time_unique,
        y1mean = y1mean, y0mean = y0mean
      )

      ytmp <- tmp[[1]]
      ji   <- tmp[[2]][,1] + 1
      ydiff_mat <- t(diff.zoo(zoo(t(ytmp)), differences = k-1, na.pad = TRUE))
      ydiff <- sapply(1:nrow(ydiff_mat), function(i) ydiff_mat[i,ji[i]])
      ## need to subset ydif based on ji

      dat2[paste("yd", k-1, sep='')] <- ydiff

      if (is.null(Xcov)) {
        fm[[k]] <- as.formula(paste("yd", k-1, " ~ treatment + post + treatment * post", sep = ''))
      } else if (is.null(x_formula)){
        fm[[k]] <- as.formula(paste("yd", k-1, " ~ treatment + post + treatment * post + ",
                              paste(x_colname, collapse = "+"), sep = ''))
      } else {
        fm[[k]] <- as.formula(paste("yd", k-1, " ~ treatment + post + treatment * post + ",
                              paste(attr(terms(x_formula), 'term.labels'), collapse = "+"), sep = ''))

      }



    }

    ## output
    out[[tt]] <- list(
      "Y" = dat2$outcome,
      "D" = dat2$treatment,
      'formula' = fm,
      'pdata' = dat2
    )

    attr(out[[tt]], 'post_treat') <- post_treatment[tt]
    attr(out[[tt]], 'id_time') <- sort(unique(id_time))
  }

  return(out)
}


#' Data processing function for panel data (long format)
#' @param outcome A vector of outcome variable.
#' @param treatment A vector of treatment indicator variable.
#' @param post_treatment A vector of numeric time index for post treatment periods.
#' @param id_subject A vector of subject index variable.
#' @param id_time A vector of time index variable.
#' @param Xcov A matrix of covariates.
#' @return A list. Each element corresponds to each post-treatment period.
#'  Each element of the returned list is also a list consists of the following:
#' \itemize{
#'  \item \code{Y}: outcome vector.
#'  \item \code{D}: treatment vector.
#'  \item \code{formula}: A list of formula for \code{\link{lm}}.
#'  \item \code{pdata}: A data frame appended transformed outcomes (e.g., yd0, yd1, etc).
#'}
#' @importFrom plm pdata.frame
#' @importFrom utils getFromNamespace
#' @importFrom dplyr %>% select filter tbl_df
#' @keywords internal
did_data_panelL <- function(outcome, treatment, post_treatment, id_subject, id_time,
  Xcov, x_formula = NULL
) {

  ## na omit
  outcome <- na.omit(outcome); treatment <- na.omit(treatment)

  ## export function
  diff.pseries <- getFromNamespace("diff.pseries", "plm")

  ### ==== intput check ==== ###
  if (length(outcome) != length(treatment)) {
    stop("Length of outcome and length of treatment do not match.")
  }
  if (any(c(min(treatment) != 0, max(treatment) != 1))) {
    stop("We currently support the binary treatment. treatment vector should take either 0 or 1.")
  }
  if (is.null(id_time)) {
    stop("`id_time` should be provided when long = TRUE.")
  } else {
    id_subject <- na.omit(id_subject)
    id_time <- na.omit(id_time)
  }
  if (!is.numeric(id_time)) {
    ## this is a little sketchy conversion
    ## maybe better to ask users to make the input numeric beforehand
    id_time <- as.numeric(as.character(id_time))
  }
  if (is.null(post_treatment)) {
    post_treatment <- max(id_time)
    cat('We treat', post_treatment, 'as the post treatment period\n')
  } else {
    post_treatment <- sort(post_treatment)
  }
  if (min(length(outcome), length(treatment), length(id_subject), length(id_time))  !=
      max(length(outcome), length(treatment), length(id_subject), length(id_time))) {
    stop("outcome, treatment, id_subject and id_time should have the same length.")
  }

  if (!is.null(Xcov)) {  ## use this when Xcov is supplied
    if (is.null(colnames(Xcov))) {
      colnames(Xcov) <- paste("XR", 1:ncol(Xcov), sep = '')
    }
    if (any(colnames(Xcov) %in% c('', "", " ", ' '))) {
      ## rename colnames if they contain empty name
      rename_idx <- which(colnames(Xcov) %in% c('', "", " ", ' '))
      colnames(Xcov)[rename_idx] <- paste('XR', length(rename_idx), sep = '')
    }

    x_colname <- colnames(Xcov)
  }

  ### ==== transform data ==== ###
  dat <- data.frame(outcome, treatment, id_subject, id_time) %>% tbl_df()
  out <- list()
  for (tt in 1:length(post_treatment)) {
    ## prep for nonparametric version
    dat_tmp <- dat %>%
      dplyr::filter(id_time < post_treatment[1] | id_time == post_treatment[tt]) %>%
      na.omit()
    y_tmp <- dat_tmp %>% dplyr::select(-treatment) %>% tidyr::spread(id_time, outcome) %>%
      dplyr::select(-id_subject) %>% data.matrix()
    d_tmp <- dat_tmp %>% select(-outcome) %>% tidyr::spread(id_time, treatment) %>%
      dplyr::select(-id_subject) %>% data.matrix()

    ## prep for parametric model

    # create dynamic treatment variable
    ## overwrite time index so that there's no time gap
    ## this is necesary to take diff in plm
    if (is.null(Xcov)) {
      dat2 <- dat %>%
        dplyr::filter(id_time < post_treatment[1] | id_time == post_treatment[tt]) %>%
        mutate(id_time = ifelse(id_time == post_treatment[tt], post_treatment[1], id_time)) %>%
        mutate(id_time2 = as.numeric(as.factor(id_time))) %>%
        na.omit()
    } else {
      dat2 <- data.frame(outcome, treatment, id_subject, id_time, Xcov) %>%
        tbl_df() %>%
        dplyr::filter(id_time < post_treatment[1] | id_time == post_treatment[tt]) %>%
        mutate(id_time = ifelse(id_time == post_treatment[tt], post_treatment[1], id_time)) %>%
        mutate(id_time2 = as.numeric(as.factor(id_time))) %>%
        na.omit()
    }
    max_time2 <- max(dat2$id_time2)

    dat_plm <- pdata.frame(dat2, index = c("id_subject", "id_time2"))

    # make formula
    pre_time <- ncol(y_tmp) - 1;
    max_diff <- pre_time - 1 ## make sure this is positive
    fm <- list()

    for (i in 0:max_diff) {
      if (is.null(Xcov)) {
        fm[[i+1]] <- as.formula(paste("yd",i, "~ treatment", sep = ''))
      } else if (is.null(x_formula)){
        fm[[i+1]] <- as.formula(paste("yd",i, " ~ treatment + ",
                        paste(x_colname, collapse = "+"), sep = ""))
      } else {
        fm[[i+1]] <- as.formula(paste("yd",i, " ~ treatment + ",
                        paste(attr(terms(x_formula), 'term.labels'), collapse = "+"), sep = ""))
      }
      if (i > 0) {
        dat_plm[paste("yd", i, sep = '')] <- diff.pseries(dat_plm$outcome, lag = i)
      } else {
        dat_plm[paste("yd", i, sep = '')] <- dat_plm$outcome
      }

    }

    ## change id for time to be the original scale
    # dat_plm_out <- pdata.frame(data.frame(dat_plm), index = c("id_subject", "id_time"))

    ## output
    out[[tt]] <- list(
      "Y" = y_tmp,
      "D" = apply(d_tmp, 1, max),
      'formula' = fm,
      'pdata' = dat_plm
    )

    attr(out[[tt]], 'post_treat') <- post_treatment[tt]
    attr(out[[tt]], 'id_time') <- sort(unique(id_time))
  }

  return(out)
}





#' Data processing function for panel data (wide format)
#' @param outcome A vector of outcome variable.
#' @param treatment A vector of treatment indicator variable.
#' @param post_treatment A vector of numeric time index for post treatment periods.
#' @return A list. Each element corresponds to each post-treatment period.
#'  Each element of the returned list is also a list consists of the following:
#' \itemize{
#'  \item \code{Y}: outcome vector.
#'  \item \code{D}: treatment vector.
#'}
#' @importFrom dplyr %>% select tbl_df
#' @keywords internal
did_data_panelW <- function(outcome, treatment, post_treat) {

  ## if the data is already in the wide format,
  ## we will use the input as is

  ## ==== input checks ==== ##
  t_time <- ncol(outcome); n_obs  <- nrow(outcome)
  if ((length(treatment) - n_obs) != 0) {
    stop('Length of the treatment vector and the number of rows in Y do not match.')
  }

  ### ==== make data into did_double obj ==== ###
  if (!is.null(post_treatment)) {
    matched_var <- as.character(post_treatment) %in% colnames(outcome)
    outcome     <- outcome %>% tbl_df()
    Ypre        <- outcome %>% select(-post_treatment) %>% data.matrix()

    if(any(matched_var)) {
      out   <- list()
      Ypost <- outcome %>% select(post_treatment[matched_var]) %>% data.matrix()
      for (i in 1:sum(matched_var)) {
        out[[i]] <- list(
          "Y" = cbind(Ypre, Ypost[,i]),
          "D" = treatment
        )

        attr(out[[i]], 'post_treat') <- post_treatment[i]
        attr(out[[i]], 'id_time') <- NULL
      }
    } else {
      stop("colnames of outcome and variables in post_treatment do not match.")
    }
  } else {
    stop("Please specify post_treatment by providing variable names.")
  }

  return(out)
}
