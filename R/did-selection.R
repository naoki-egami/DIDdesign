
#' Function to test directin of trends. 
#'
#' This function test if pre-treatment trends are in the same direction. 
#' Currently we support when there are only two pre-treatment time points.
#' @param dat A data object of a class \code{diddesign_data}. Typically an output from \code{did_data}.
#' @param level level of rejection. Default is 0.1.
#' @examples
#' # load data from didrobust package
#' data(anzia2012)
#'
#' # create did_data object 
#' out <- did_data(
#'   outcome = anzia2012$lnavgsalary_cpi,
#'   treatment = anzia2012$oncycle,
#'   post_treatment = c(2007, 2008, 2009),
#'   id_subject = anzia2012$district,
#'   id_time = anzia2012$year,
#'   long = TRUE
#' )
#'
#' # run the test 
#' tt <- sign_test(out)
#' tt$pvalue 
#' @export
sign_test <- function(dat, level = 0.1) {
  dat_type <- attr(dat, 'data_type')
  if (dat_type == "panel") {
    fit <- sign_test_panel(dat, level = level)
  } else {
    ## RCS 
    fit <- sign_test_rcs(dat, level = level)
  }
  
  return(fit)
}



#' Function to test directin of trends for panel 
#' @keywords internal
sign_test_panel <- function(dat, level = 0.1) {


  ## get pretreatment outcomes 
  Ypre <- dat[[1]]$Y[,-ncol(dat[[1]]$Y)]
  Tpre <- ncol(Ypre)  
  Dvec <- dat[[1]]$D
  
  ## get data info
  n01 <- table(Dvec)
  
  if (Tpre >= 3) {
    cat("... conducting a test to assess PTT ... \n")
    Ypre <- Ypre[,c(Tpre-1, Tpre)]
  }
  
  
  ## get the "true" sign 
  Ymean <- apply(Ypre, 2, function(x) tapply(x, Dvec, mean, na.rm = TRUE))
  trend <- apply(Ymean, 1, diff)
  sign_true <- trend[which.max(n01)]
  
  ## compute the bias 
  bias <- diff(Ymean[,2])
  
  ## test 
  # we pick the opposite of "true" sign as H1 (so that two failture to reject H0 is supportive to sign_true)
  h1 <- ifelse(sign_true > 0, 'less', 'greater')
  # 1. control group 
  D0_res <- t.test(apply(Ypre[Dvec==0,], 1, diff), alternative = h1)
  # 2. treatment group 
  D1_res <- t.test(apply(Ypre[Dvec==1,], 1, diff), alternative = h1)
  
  # 3. combine the two tests 
  pval <- min(D0_res$p.value, D1_res$p.value)
  res  <- ifelse(pval > level, 'pass', 'fail to pass')
  
  ## output
  return(list(res = res, pvalue = pval, pval0 = D0_res$p.value, 
              pval1 = D1_res$p.value, T0 = D0_res, T1 = D1_res, bias = bias))
  
}


#' Function to test direction of trends for RCS 
#' @importFrom dplyr %>% ungroup select summarise group_by filter 
#' @keywords internal
sign_test_rcs <- function(dat, level = 0.1) {
  post_time <- attr(dat[[1]], 'post_treat')
  full_time <- attr(dat[[1]], 'id_time')
  pre_time  <- full_time[full_time != post_time]
  dat[[1]]$pdata %>% select(-post_time) %>% 
    group_by(id_time, treatment) %>% 
    summarise(Ymean = mean(outcome), Ysd = sd(outcome), n = n()) %>% 
    ungroup() -> tmp 
    
  if (length(unique(tmp$id_time)) > 3) {
    cat("... Using only the last two periods to assess PTT ... \n")
    tmp <- tmp %>% filter(id_time %in% tail(pre_time))
    pre_time <- unique(tmp$id_time)
  }
  
  ## compute 'true' sign 
  trend <- c(tmp$Ymean[3] - tmp$Ymean[1], tmp$Ymean[4] - tmp$Ymean[2])
  sign_true <- trend[which.max(c(tmp$n[3]+tmp$n[1], tmp$n[4] +tmp$n[2]))]
  
  # compute bias 
  bias <- tmp$Ymean[4] - tmp$Ymean[3]
  
  # we pick the opposite of "true" sign as H1 (so that two failture to reject H0 is supportive to sign_true)
  h1 <- ifelse(sign_true > 0, 1, 2)

  T0 <- trend[1] / sqrt(tmp$Ysd[3]^2+tmp$Ysd[1]^2)
  T1 <- trend[2] / sqrt(tmp$Ysd[4]^2+tmp$Ysd[2]^2)
  
  pval0 <- c(pnorm(T0), 1 - pnorm(T0))
  pval1 <- c(pnorm(T1), 1 - pnorm(T1))
  pvalue <- min(pval0[h1], pval1[h1])
  res  <- ifelse(pvalue > level, 'pass', 'fail to pass')

  return(list(res = res, pvalue = pvalue, pval0 = pval0[h1], pval1 = pval1[h1], T0 = T0, T1 = T1, bias = bias))

  
}



sign_test_parametric_rcs <- function(coefs, vcov, level = 0.1) {
  ## trends 
  trend0 <- coefs[1]; trend1 <- trend0 + coefs[2]
  T0 <- trend0 / sqrt(vcov[1,1])
  T1 <- trend1 / sqrt(vcov[1,1] + vcov[1,2] + 2 * vcov[2,2])
  
  sign_true <- trend0
  h1 <- ifelse(sign_true > 0, 1, 2)
  pval0 <- c(pnorm(T0), 1 - pnorm(T0))
  pval1 <- c(pnorm(T1), 1 - pnorm(T1))
  pvalue <- min(pval0[h1], pval1[h1])
  res  <- ifelse(pvalue > level, 'pass', 'fail to pass')
  return(list(res = res, pvalue = pvalue, pval0 = pval0[h1], pval1 = pval1[h1], T0 = T0, T1 = T1, bias = coefs[2]))  
}


#' equivalence check 
equivalence_test <- function(theta, se, eq, level = 0.1) {
  res <- rep(NA, length(theta))
  for (i in 1:length(theta)) {
    UB <- theta[i] - qnorm(level) * se[i] 
    LB <- theta[i] + qnorm(level) * se[i] 
    res[i] <- ifelse((UB <= abs(eq)) && (LB >= -abs(eq)), "pass", "fail to pass")    
  }
  return(res)
}
