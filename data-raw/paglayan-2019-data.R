require(haven)
require(tidyverse)

## load data
dat <- read_dta("Paglayan Dataset.dta")

dat %>%
  select(CBeverrequired, CBrequired_SY, year, State) %>%
  pivot_wider(names_from = "year",
              values_from = "CBrequired_SY",
              names_prefix = "year")


## full data
full_hist <- ggplot(dat, aes(x = year, y = State, fill = CBrequired_SY)) +
  geom_tile()

after60_hist <- dat %>%
  filter(year >= 1959) %>%
  ggplot(aes(x = year, y = State, fill = CBrequired_SY)) +
    geom_tile() +
    ggtitle("1959 - 2000")

ggsave(after60_hist, file = '~/Desktop/hist_map.pdf')



paglayan2019 <- dat %>%
  filter(year >= 1959) %>%
  select(year, State, CBrequired_SY, ppexpend, avgteachsal) %>%
  rename(state = State, treatment = CBrequired_SY,
        pupil_expenditure = ppexpend, teacher_salary = avgteachsal)

usethis::use_data(paglayan2019)
