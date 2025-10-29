### Calf Survival and Reproduction Data Management
### Last updated: Oct. 27, 2025
### Contact: xprockox@gmail.com

############################################################
### ------------------- PACKAGES ----------------------- ###
############################################################

library(dplyr)
library(lubridate)

############################################################
### ------------------ DATA IMPORT --------------------- ###
############################################################

classification <- read.csv('data/master/elk_classification_2025-10-17.csv')
counts <- read.csv('data/master/elk_counts_2023-07-09.csv')
preg_collar <- read.csv('data/master/elk_pregnancy_collars_2025-07-18.csv')
preg_harvest <- read.csv('data/master/elk_pregnancy_harvest_2015-03-24.csv')
age_structure <- read.csv('data/master/elk_harvest_2021.csv')

############################################################
### ---------------- DATA MANAGEMENT ------------------- ###
############################################################

classification <- classification %>%
  rename(year = Year,
         survey_date = Survey.Date,
         total_elk = Total..Elk.Classified,
         cows = X.Cows,
         calves = X.Calves,
         calf_per100cow = X.Calves..100.Cows,
         spike_per100cow = X.Yearling.bulls..100.cows,
         btb_per100cow = X.BTB..100.cows,
         bull_per100cow = X.Bulls.100.cows)

counts <- counts %>%
  rename(year = winter..Jan..,
         n_total = mean, 
         n_total_LCI = lwr.CL,
         n_total_UCI = upr.CL) %>%
  select(year, n_total, n_total_LCI, n_total_UCI)

df <- left_join(classification, counts)

df$calf_cow_ratio <- round(df$calf_per100cow / 100, 2)
df$spike_cow_ratio <- round(df$spike_per100cow / 100, 2)
df$btb_cow_ratio <- round(df$btb_per100cow/ 100, 2)
df$bull_cow_ratio <- round(df$bull_per100cow / 100, 2)

df$total_elk <- as.numeric(df$total_elk)

df$percent_cows <- df$cows / df$total_elk

df$n_cows <- round(df$n_total * df$percent_cows)
df$n_calves <- round(df$n_total * df$percent_cows * df$calf_cow_ratio)
df$n_spikes <- round(df$n_total * df$percent_cows * df$spike_cow_ratio)
df$n_btb <- round(df$n_total * df$percent_cows * df$btb_cow_ratio)
df$n_bulls <- round(df$n_total * df$percent_cows * df$bull_cow_ratio)

# incorporating park preg data from collared indivs

preg_collar_prop <- preg_collar %>%
  mutate(
    # Convert pregnant to binary 1/0
    pregnant_code = ifelse(Pregnant == "yes", 1, 0),
    
    # Extract year of capture
    year = year(as.POSIXct(strptime(Capture.Date, format = '%d-%b-%y'))),

    # Calculate age at capture
    age = year - BirthYear,
    
    # Define age classes
    AgeClass = ifelse(age <= 13, "young", "old")
  ) %>%
  group_by(year, AgeClass) %>%
  summarize(
    prop_pregnant = mean(pregnant_code, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  pivot_wider(
    names_from = AgeClass,
    values_from = prop_pregnant,
    names_glue = "{AgeClass}_prop_preg_collar"
  ) %>%
  arrange(year)%>%
  select(-NA_prop_preg_collar)

df <- left_join(df, preg_collar_prop)

# incorporating FWP preg data from harvested elk during antlerless harvest years (1997-2009)

preg_harvest_props <- preg_harvest %>%
  mutate(
    AgeClass = ifelse(ageatharvest <= 13, "young", "old")
  ) %>%
  group_by(harvestyear, AgeClass) %>%
  summarize(
    prop_pregnant = mean(pregnant_code, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  pivot_wider(
    names_from = AgeClass,
    values_from = prop_pregnant,
    names_glue = "{AgeClass}_prop_preg_harvest"
  ) %>%
  arrange(harvestyear)%>%
  rename(year = harvestyear)

preg_harvest_props

df <- left_join(df, preg_harvest_props, by='year')

### How many cows are young vs. old? Use harvest data from FWP
age_structure_prop <- age_structure %>%
  pivot_longer(
    cols = starts_with("Y"),
    names_to = "Year",
    values_to = "Count"
  ) %>%
  mutate(
    Year = str_remove(Year, "Y") %>% as.integer(),
    AgeClass = ifelse(Age <= 13, "young", "old")
  ) %>%
  group_by(Year, AgeClass) %>%
  summarize(Total = sum(Count), .groups = "drop") %>%
  group_by(Year) %>%
  mutate(Proportion = Total / sum(Total)) %>%
  select(Year, AgeClass, Proportion) %>%
  pivot_wider(
    names_from = AgeClass,
    values_from = Proportion,
    names_glue = "{AgeClass}_propN"
  ) %>%
  arrange(Year)%>%
  rename(year = Year)

df <- left_join(df, age_structure_prop, by='year')

# merges the two data sources of pregnancy (capture and harvest), 
# including harvest estimate only when the capture estimate does not exist
df <- df %>%
  rename(total_elk_classified = total_elk) %>%
  mutate(young_prop_preg = round(coalesce(young_prop_preg_collar, young_prop_preg_harvest), 2),
         old_prop_preg = round(coalesce(old_prop_preg_collar, old_prop_preg_harvest), 2),
         n_cows_young = round(n_cows * young_propN),
         n_cows_old = round(n_cows * old_propN),
         expected_calves_from_young = round(n_cows_young * young_prop_preg, 2),
         expected_calves_from_old = round(n_cows_old * old_prop_preg, 2),
         expected_calves = expected_calves_from_young + expected_calves_from_old,
         calf_surv = round(n_calves/expected_calves, 2)) %>%
  select(-survey_date, -Method, 
         -total_elk_classified, -cows, -calves,
         -calf_per100cow, -spike_per100cow, -btb_per100cow, -bull_per100cow, 
         -old_prop_preg_collar, -old_prop_preg_harvest,
         -young_prop_preg_collar, -young_prop_preg_harvest)

############################################################
### ------------------ DATA WRITING -------------------- ###
############################################################

stop('[1.3.mgmt_productivity.R] \n
All required matrices constructed. Code stopped to prevent overwriting data.\n
Continue running code beyond this line to overwrite data exports.')

write.csv(df, 'data/intermediate/productivity.csv')
