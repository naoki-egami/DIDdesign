

##
## Sances (2015) "The Distributional Impact of Greater Responsiveness" _Journal of Politics_
## Original dataset:
##  - main.dta
##

require(haven)
require(dplyr)


## load data
dat <- read_dta("main.dta")




## subset dataset
sances2015 <- dat %>%
  select(swis_code, year, muni_name, county_name, reassessed, elected) %>%
  mutate(treatment = if_else(elected == 1, 0, 1))%>%
  select(-elected)


## save
usethis::use_data(sances2015)
