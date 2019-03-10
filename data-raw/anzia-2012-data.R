
# xi: xtreg lnavgsalary_cpi oncycle i.year if commondist==0 & typeunknown==0, fe cluster(district)
require(dplyr)
require(readstata13)

## load data 
dat <- read.dta13("panel data 3.22.12.dta") %>% tbl_df()

## prep for saving 
anzia2012 <- dat %>% 
  filter(commondist==0 & typeunknown==0) %>%
  select(district, year, lnavgsalary_cpi, oncycle,
  lnenrollment, lnt1_cpi, teachers_avg_yrs_exper, ami_pc, asian_pc, black_pc, hisp_pc
  ) %>% na.omit()

## save 
usethis::use_data(anzia2012)
