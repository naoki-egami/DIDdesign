

## test


 paglayan2019 <- paglayan2019 %>%
 mutate(id_time = year,
       id_subject = as.numeric(as.factor(state)),
       log_expenditure = log(pupil_expenditure + 1),
       log_salary      = log(teacher_salary + 1)) %>%
 filter(!(state %in% c("WI", "DC")))


## test create_Gmat
require(tidyr)
tmp <- paglayan2019 %>% group_by(id_subject) %>% nest()
tmp$data[[1]] %>% filter(treatment == 1)
tmp$data[[2]] %>% filter(treatment == 1)
tmp$data[[5]] %>% filter(treatment == 1)


GG <- create_Gmat(paglayan2019, "treatment")
GG[c(1,2,5),]


## test get_periods
get_periods(GG, thres = 1)
which(apply(GG, 2, function(x) sum(x == 1)) >= 1)

get_periods(GG, thres = 2)
which(apply(GG, 2, function(x) sum(x == 1)) >= 2)


## test get subjects
GG[,which(apply(GG, 2, function(x) sum(x == 1)) >= 4)]

apply(
  GG[,which(apply(GG, 2, function(x) sum(x == 1)) >= 4)],
  2, function(x) which(x == 1 | x == 0)
)
get_subjects(GG, get_periods(GG, thres = 4))
