

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
#' @return A \code{diddesign} class object.
#' @importFrom dplyr %>% tbl_df select filter mutate
#' @importFrom tidyr spread
#' @importFrom plm pdata.frame
#' @importFrom utils getFromNamespace
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
#'  # make a plot
#'  plot(out)
#'
#'  # view a summary
#'  summary(out)
#' @export
did_data <- function(
  outcome, treatment,
  post_treatment = NULL, id_subject = NULL, id_time = NULL, long = TRUE,
  Xcov = NULL
) {

  if (is.null(id_subject)) {
    is_rcs <- TRUE
    # message("treat data as repeated cross-section.\n")
  } else {
    is_rcs <- FALSE
  }

  ## data transformation
  if (isTRUE(is_rcs)) {

    out <- did_data_rcs(outcome, treatment, post_treatment, id_time, Xcov)

  } else if (!isTRUE(is_tcs) & isTRUE(long)) {
    ## if the input is in the long format,
    ## we transform it to the wide format

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
          na.omit()
      } else {
        dat2 <- data.frame(outcome, treatment, id_subject, id_time, Xcov) %>%
          tbl_df() %>%
          dplyr::filter(id_time < post_treatment[1] | id_time == post_treatment[tt]) %>%
          mutate(id_time = ifelse(id_time == post_treatment[tt], post_treatment[1], id_time)) %>%
          na.omit()
      }
      dat_plm <- pdata.frame(dat2, index = c("id_subject", "id_time"))

      # make formula
      pre_time <- ncol(y_tmp) - 1;
      max_diff <- pre_time - 1 ## make sure this is positive
      fm <- list()

      for (i in 0:max_diff) {
        if (is.null(Xcov)) {
          fm[[i+1]] <- as.formula(paste("yd",i, "~ treatment", sep = ''))
        } else {
          fm[[i+1]] <- as.formula(paste("yd",i, " ~ treatment + ", paste(x_colname, collapse = "+"), sep = ""))
        }
        if (i > 0) {
          dat_plm[paste("yd", i, sep = '')] <- diff.pseries(dat_plm$outcome, lag = i)
        } else {
          dat_plm[paste("yd", i, sep = '')] <- dat_plm$outcome
        }

      }

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

  } else {
    ## if the data is already in the wide format,
    ## we will use the input as is

    ## ==== input checks ==== ##
    t_time <- ncol(Y); n_obs  <- nrow(Y)
    if ((length(D) - n_obs) != 0) {
      stop('Length of the treatment vector and the number of rows in Y do not match.')
    }

    ### ==== make data into did_double obj ==== ###
    if (!is.null(post_treatment)) {
      matched_var <- as.character(post_treatment) %in% colnames(outcome)
      outcome <- outcome %>% tbl_df()
      Ypre <- outcome %>% select(-post_treatment) %>% data.matrix()

      if(any(matched_var)) {
        out <- list()
        Ypost <- outcome %>% select(post_treatment[matched_var]) %>% data.matrix()
        for (i in 1:sum(matched_var)) {
          out[[i]] <- list(
            "Y" = cbind(Ypre, Ypost[,i]),
            "D" = treatment
          )
        }
      } else {
        stop("colnames of outcome and variables in post_treatment do not match.")
      }
    } else {
      stop("Please specify post_treatment by providing variable names.")
    }
  }

  class(out) <- c("diddesign", "diddesign_data")
  return(out)
}




#' data processing function for repeated cross-section data
#' @param outcome a vector of outcome observations.
#' @importFrom zoo zoo
#' @importFrom dplyr %>% filter summarise group_by pull tbl_df
#' @importFrom utils getFromNamespace
#' @keywords internal
did_data_rcs <- function(outcome, treatment, post_treatment, id_time, Xcov) {

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
    } else {
      fm[[1]] <- as.formula(paste("yd",0, " ~ treatment + post + treatment * post + ",
                            paste(x_colname, collapse = "+"), sep = ''))
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
      ji    <- tmp[[2]][,1] + 1
      # for (i in 1:nrow(dat2)) {
      #   di <- dat2$treatment[i]
      #   ji <- which(id_time_unique == dat2$id_time[i])
      #   ytmp <- di * y1mean + (1 - di) * y0mean
      #   ytmp[ji] <- dat2$outcome[i]
      #   ydiff[i] <- as.vector(diff.zoo(zoo(ytmp), differences = k-1, na.pad=TRUE))[ji]
      # }


      ydiff_mat <- t(diff.zoo(zoo(t(ytmp)), differences = k-1, na.pad = TRUE))
      ydiff <- sapply(1:nrow(ydiff_mat), function(i) ydiff_mat[i,ji[i]])
      ## need to subset ydif based on ji

      dat2[paste("yd", k-1, sep='')] <- ydiff

      if (is.null(Xcov)) {
        fm[[k]] <- as.formula(paste("yd", k-1, " ~ treatment + post + treatment * post", sep = ''))
      } else {
        fm[[k]] <- as.formula(paste("yd", k-1, " ~ treatment + post + treatment * post + ",
                              paste(x_colname, collapse = "+"), sep = ''))
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


make_basis <- function(k, total_time) replace(numeric(total_time), k, 1)
