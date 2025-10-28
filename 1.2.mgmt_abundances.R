### Elk Abundance - Data management
### Last updated: Oct. 27, 2025
### Contact: xprockox@gmail.com

############################################################
### ------------------ PACKAGES ------------------------ ###
############################################################

library(tidyverse)

############################################################
### ------------------ DATA LOADING -------------------- ###
############################################################

abundances <- read.csv('data/master/elk_counts_2023-07-09.csv')
classification <- read.csv('data/master/elk_classification_2025-10-17.csv')
harvest <- read.csv('data/master/elk_harvest_2021.csv')

############################################################
### ------------------- MANAGEMENT --------------------- ###
############################################################

abundances <- abundances %>%
  select(winter..Jan.., mean, lwr.CL, upr.CL) %>%
  rename(year = winter..Jan..,
         n_total = mean,
         n_total_L95 = lwr.CL,
         n_total_U95 = upr.CL)

classification <- classification %>%
  rename(year = Year,
         n_cow = X.Cows,
         n_calf = X.Calves) %>%
  mutate(total_elk_classified = as.integer(Total..Elk.Classified),
         n_spike = round(n_cow * (X.Yearling.bulls..100.cows/100)),
         n_btb = round(n_cow * (X.BTB..100.cows/100)),
         n_bull = round(n_cow * (X.Bulls.100.cows/100)),
         percent_cow = round(n_cow/total_elk_classified,2),
         percent_calf = round(n_calf/total_elk_classified,2),
         percent_spike = round(n_spike/total_elk_classified,2),
         percent_btb = round(n_btb/total_elk_classified,2),
         percent_bull = round(n_bull/total_elk_classified,2),
         check = percent_cow + percent_calf + percent_spike + percent_btb,
         check2 = n_cow + n_calf + n_bull,
         check2 = total_elk_classified - check2) %>%
  select(year, total_elk_classified, 
         n_cow, n_calf, n_spike, n_btb, n_bull,
         percent_cow, percent_calf, percent_spike, percent_btb, percent_bull, 
         check, check2)

class_percentages <- classification %>%
  select(year, percent_cow, percent_calf, percent_spike, percent_btb, percent_bull)

abundances <- left_join(abundances, class_percentages)

abundances <- abundances %>%
  mutate(n_cow = n_total * percent_cow,
         n_calf = n_total * percent_calf,
         n_spike = n_total * percent_spike,
         n_btb = n_total * percent_btb,
         n_bull = n_total * percent_bull,
         sigma_tot_log = round((log(n_total_U95) - log(n_total_L95)) / (2*1.96),2))

harvest <- harvest %>%
  rename(age = 1) %>%
  pivot_longer(
    cols = -age,
    names_to = "year_raw",
    values_to = "harvest"
  ) %>%
  mutate(
    year = parse_number(year_raw),
    group = case_when(
      age == 1 ~ "age1",
      age >= 2 & age <= 13 ~ "age2_13",
      age >= 14 ~ "age14plus"
    )
  ) %>%
  group_by(year, group) %>%
  summarise(n = sum(harvest, na.rm = TRUE), .groups = "drop_last") %>%
  mutate(total_year = sum(n), prop = n / total_year) %>%
  ungroup() %>%
  select(year, group, prop) %>%
  pivot_wider(names_from = group, values_from = prop) %>%
  arrange(year) %>%
  rename(prop_age1 = age1,
         prop_age2_13 = age2_13,
         prop_age14plus = age14plus)

abundances <- left_join(abundances, harvest)

############################################################
### ------------------ WRITE DATA ---------------------- ###
############################################################

stop('[1.2.mgmt_abundances.R] \n
All required matrices constructed. Code stopped to prevent overwriting data.\n
Continue running code beyond this line to overwrite data exports.')

write.csv(abundances, 'data/intermediate/abundanceEstimates_stages.csv')
