### Elk Abundance - Data management
### Last updated: Oct. 27, 2025
### Contact: xprockox@gmail.com

############################################################
### ------------------ PACKAGES ------------------------ ###
############################################################

library(tidyverse)
library(zoo)

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
  mutate(n_cow = round(n_total * percent_cow),
         n_calf = round(n_total * percent_calf),
         n_spike = round(n_total * percent_spike),
         n_btb = round(n_total * percent_btb),
         n_bull = round(n_total * percent_bull),
         sigma_tot_log = round((log(n_total_U95) - log(n_total_L95)) / (2*1.96),2))

############################################################
########--------------- DATA VIZ ---------------############
############################################################

# visualize timeseries of young vs. old adults harvested
harvest_long <- harvest %>%
  pivot_longer(
    cols = starts_with("Y"),
    names_to = "Year",
    values_to = "Harvest"
  ) %>%
  mutate(
    Year = as.numeric(sub("Y", "", Year)),
    AgeClass = case_when(
      Age >= 1 & Age <= 13 ~ "Young Adult",
      Age >= 14 ~ "Old Adult"
    )
  ) %>%
  group_by(Year, AgeClass) %>%
  summarise(TotalHarvest = sum(Harvest, na.rm = TRUE), .groups = "drop")

p1 <- ggplot(harvest_long, aes(x = Year, y = TotalHarvest)) +
  geom_line(size = 1.1) +
  facet_wrap(~ AgeClass, ncol = 1, scales = "free_y") +
  labs(
    x = "Year",
    y = "Number Harvested"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    strip.text = element_text(face = "bold"),
    panel.spacing = unit(1.2, "lines")
  )
p1 

# filter abundance data down to 1995:2021 because that's the time range of the harvest data
dat_n_2021 <- abundances %>%
  filter(year %in% 1995:2021) %>%
  arrange(year)

# estimate the percent cows in the total elk population using:
# 1. % cow estimates from classification flights
# 2. interpolated % cow estimates for years where flights did not happen
# 3. total elk numbers multiplied by % cow estimates
dat_n_2021$percent_cow_filled <- zoo::na.approx(dat_n_2021$percent_cow)
dat_n_2021$n_cow <- dat_n_2021$percent_cow_filled * dat_n_2021$n_total

harvest_long <- harvest_long %>%
  left_join(
    dat_n_2021 %>% select(year, n_cow),
    by = c("Year" = "year")
  )

# plot proportion of age classes harvested each year
harvest_prop <- harvest_long %>%
  group_by(Year) %>%
  mutate(
    YearTotal = sum(TotalHarvest, na.rm = TRUE),
    Prop = TotalHarvest / YearTotal
  ) %>%
  ungroup()

p2 <- ggplot(harvest_prop, aes(x = Year, y = Prop)) +
  geom_line(size = 1.1) +
  facet_wrap(~ AgeClass, ncol = 1, scales = "free_y") +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  labs(
    x = "Year",
    y = "Proportion of Total Harvest"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    strip.text = element_text(face = "bold"),
    panel.spacing = unit(1.2, "lines")
  )
p2

# estimate the number of individuals in old vs. young classes for the NR female elk population 
# by multiplying n_cow by proportions of each class harvested
harvest_prop <- harvest_prop %>%
  mutate(
    pop_structure = Prop * n_cow
  )

p3 <- ggplot(harvest_prop, aes(x = Year, y = pop_structure)) +
  geom_line(size = 1.1) +
  facet_wrap(~ AgeClass, ncol = 1, scales = "free_y") +
  labs(
    x = "Year",
    y = "Stage-Specific Female Population Estimates"
  ) +
  theme_minimal(base_size = 14)
p3

# Plot all together
cowplot::plot_grid(p1, p2, p3, ncol=3)

############################################################
#####----------- (MORE) DATA MANAGEMENT ---------------#####
############################################################

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
      age >= 1 & age <= 13 ~ "young_ad",
      age >= 14 ~ "old_ad"
    )
  ) %>%
  group_by(year, group) %>%
  summarise(n = sum(harvest, na.rm = TRUE), .groups = "drop_last") %>%
  mutate(total_year = sum(n), prop = n / total_year) %>%
  ungroup() %>%
  select(year, group, prop) %>%
  pivot_wider(names_from = group, values_from = prop) %>%
  arrange(year) %>%
  rename(prop_young_ad = young_ad,
         prop_old_ad = old_ad)

abundances <- left_join(abundances, harvest)

abundances <- abundances %>%
  mutate(n_cow_youngadult = round(n_cow * prop_young_ad),
         n_cow_oldadult = round(n_cow * prop_old_ad))

############################################################
### ------------------ WRITE DATA ---------------------- ###
############################################################

stop('[1.2.mgmt_abundances.R] \n
All required matrices constructed. Code stopped to prevent overwriting data.\n
Continue running code beyond this line to overwrite data exports.')

write.csv(abundances, 'data/intermediate/abundanceEstimates_stages.csv')
