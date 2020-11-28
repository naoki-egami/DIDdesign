# #########################################
# Replication of Malesky et al (2014; APSR)
# #########################################

rm(list=ls())

## load other packages
require(tidyverse)
require(plm)
require(haven)
require(estimatr)


##
## load data 
##

dat_2006_2008 <- read_dta("panel_commune_2006_2008.dta")
dat_2008_2010 <- read_dta("panel_commune_2008_2010.dta")

dat_2006 <- filter(dat_2006_2008, year == 2006)
dat_2008 <- filter(dat_2006_2008, year == 2008)
dat_2010 <- filter(dat_2008_2010, year == 2010)


share <- intersect(colnames(dat_2008), colnames(dat_2010))
dat_main <- bind_rows(
  dat_2006 %>% select(all_of(share)),
  dat_2008 %>% select(all_of(share)),
  dat_2010 %>% select(all_of(share))
)


# 30 main outcomes analyzed in their paper
index1 <- c("goodroadv", "transport", "pro3", "tapwater", "roadv")
index2 <- c("rm2c7d", "rm2c7e", "rm2c7g", "animal_s", "agrvisit", "plant_s", "agrext", "irrigation")
index3 <- c("rm2c7c", "pro5")
index4 <- c("pro4", "rm2c7b", "useschool", "kgarten", "v_prischool")
index5 <- c("broadcast", "post", "vpost")
index6 <- c("rm2c7a", "rm2c7f", "market", "nonfarm", "vmarket1", "vmarket2", "vmarket3")
outcome_list <- c(index1, index2, index3, index4, index5, index6)

##
## Data prepration (following the original replication file)
##
dat_main$agrvisit <- dat_main$agrvisit/100

dat_main$city <- NA
for(i in 1:nrow(dat_main)){
  if(dat_main$year[i] <=2008){
    dat_main$city[i] <- ifelse(dat_main$tinh[i] ==101 | dat_main$tinh[i]==103 | dat_main$tinh[i]==501 |
                                 dat_main$tinh[i]==815 | dat_main$tinh[i]==701, 1, 0)
  }else if(dat_main$year[i] > 2008){
    dat_main$city[i] <- ifelse(dat_main$tinh[i] ==1 | dat_main$tinh[i]==31 | dat_main$tinh[i]== 48 |
                                 dat_main$tinh[i]== 92 | dat_main$tinh[i]== 79, 1, 0)
  }
}
# remove some region (following replication file)
dat_main <- dat_main[dat_main$reg8 !=6, ]


##
## subset data for saving 
##
malesky2014 <- dat_main %>% 
  select(all_of(outcome_list), treatment, 
          lnarea, lnpopden, city, reg8, rnongnghiep, rcongnghiep, rdichvu, 
          year, tinh) %>% 
  mutate(post_treat = ifelse(year == 2010, 1, 0))
usethis::use_data(malesky2014, overwrite = TRUE)
