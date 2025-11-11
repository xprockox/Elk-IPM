### Calf Survival and Reproduction Data Management
### Last updated: Nov. 4, 2025
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
    pregnant_code = as.integer(Pregnant == "yes"),
    year = lubridate::year(as.POSIXct(strptime(Capture.Date, "%d-%b-%y"))),
    age = year - BirthYear,
    AgeClass = case_when(age >=2 & age <= 13 ~ "young",
                         age >  13 ~ "old",
                         age <= 1 ~ "yearling",
                         TRUE ~ NA_character_)
  ) %>%
  filter(!is.na(AgeClass)) %>%                     # avoid NA_* columns
  group_by(year, AgeClass) %>%
  summarise(
    total = n(),
    num_pregnant   = sum(pregnant_code, na.rm = TRUE),
    prop_pregnant  = mean(pregnant_code, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  pivot_wider(
    names_from  = AgeClass,
    values_from = c(total, num_pregnant, prop_pregnant),
    names_glue  = "{AgeClass}_{.value}_collar"
  ) %>%
  arrange(year)

df <- left_join(df, preg_collar_prop, by='year')

# incorporating FWP preg data from harvested elk during antlerless harvest years (1997-2009)

preg_harvest_props <- preg_harvest %>%
  mutate(
    AgeClass = case_when(
      ageatharvest >= 2 & ageatharvest <= 13 ~ "young",
      ageatharvest >= 14 ~ "old",
      ageatharvest <= 1 ~ "yearling",
      TRUE ~ NA_character_
    )
  ) %>%
  filter(!is.na(AgeClass)) %>%
  group_by(harvestyear, AgeClass) %>%
  summarise(
    total = n(),
    num_pregnant   = sum(pregnant_code, na.rm = TRUE),
    prop_pregnant  = mean(pregnant_code, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  pivot_wider(
    names_from  = AgeClass,
    values_from = c(total, num_pregnant, prop_pregnant),
    names_glue  = "{AgeClass}_{.value}_harvest"   # <-- include {.value}
  ) %>%
  arrange(harvestyear) %>%
  rename(year = harvestyear)

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
    AgeClass = case_when(
      Age >= 2 & Age <= 13 ~ "young",
      Age >= 14 ~ "old",
      Age <= 1 ~ "yearling",
      TRUE ~ NA_character_
    )
  ) %>%
  group_by(Year, AgeClass) %>%
  summarize(Total = sum(Count), .groups = "drop") %>%
  group_by(Year) %>%
  mutate(Proportion = Total / sum(Total)) %>%
  select(Year, AgeClass, Proportion) %>%
  pivot_wider(
    names_from = AgeClass,
    values_from = Proportion,
    names_glue = "percent_cows_{AgeClass}"
  ) %>%
  arrange(Year)%>%
  rename(year = Year)

df <- left_join(df, age_structure_prop, by='year')

# merges the two data sources of pregnancy (capture and harvest), 
# including harvest estimate only when the capture estimate does not exist
df <- df %>%
  rename(total_elk_classified = total_elk) %>%
  mutate(# number pregnant - use collar data when avail., harvest data when unavail. EXCEPT yearlings
         yearling_num_preg = round(coalesce(yearling_num_pregnant_harvest, yearling_num_pregnant_collar), 2), # harvest data > collar data
         young_num_preg = round(coalesce(young_num_pregnant_collar, young_num_pregnant_harvest), 2), # collar data > harvest data
         old_num_preg = round(coalesce(old_num_pregnant_collar, old_num_pregnant_harvest), 2), # collar data > harvest data
         # number captured or harvested - same logic as above
         yearling_num_capt = round(coalesce(yearling_total_harvest, yearling_total_collar), 2), 
         young_num_capt = round(coalesce(young_total_collar, young_total_harvest), 2), 
         old_num_capt = round(coalesce(old_total_collar, old_total_harvest), 2), 
         # proportion pregnant - same logic as above
         yearling_prop_pregnant = round(coalesce(yearling_prop_pregnant_harvest, yearling_prop_pregnant_collar), 2),
         young_prop_pregnant = round(coalesce(young_prop_pregnant_collar, young_prop_pregnant_harvest), 2),
         old_prop_pregnant = round(coalesce(old_prop_pregnant_collar, old_prop_pregnant_harvest), 2),
         # number of each adult class in total population (calves = yearlings in classification survey)
         n_cows_young = round(n_cows * percent_cows_young),
         n_cows_old = round(n_cows * percent_cows_old),
         # expected calves is N of a class * pregnancy rate of that class (assumes no abortion)
         expected_calves_from_young = round(n_cows_young * young_prop_pregnant, 2),
         expected_calves_from_old = round(n_cows_old * old_prop_pregnant, 2),
         expected_calves = expected_calves_from_young + expected_calves_from_old,
         # calf survival is the number of calves (what we call yearlings) counted in survey / expected calves
         calf_surv = round(n_calves/expected_calves, 2)) %>%
  select(-survey_date, -Method, 
         -total_elk_classified, -cows, -calves,
         -calf_per100cow, -spike_per100cow, -btb_per100cow, -bull_per100cow, 
         -yearling_prop_pregnant_collar, -yearling_prop_pregnant_harvest,
         -old_prop_pregnant_collar, -old_prop_pregnant_harvest,
         -young_prop_pregnant_collar, -young_prop_pregnant_harvest,
         -yearling_total_collar, -yearling_total_harvest,
         -young_total_collar, -young_total_harvest,
         -old_total_collar, -old_total_harvest,
         -yearling_num_pregnant_collar, -yearling_num_pregnant_harvest,
         -young_num_pregnant_collar, -young_num_pregnant_harvest,
         -old_num_pregnant_collar, -old_num_pregnant_harvest)

# now calculate how many of the young adults are 13 y.o.
ya_13_prop <- age_structure %>%
  filter(Age %in% 2:13) %>%                     
  pivot_longer(
    cols = starts_with("Y"),
    names_to = "year",
    values_to = "count"
  ) %>%
  mutate(year = as.numeric(substr(year, 2, 5))) %>% 
  group_by(year) %>%
  summarize(
    harvested_age13 = sum(count[Age == 13], na.rm = TRUE),
    harvested_total = sum(count, na.rm = TRUE),
    prop_age13 = round(harvested_age13 / harvested_total, 2),
    .groups = "drop"
  ) 

df <- left_join(df, ya_13_prop, by='year')

df$n_age13 <- round(df$n_total * df$prop_age13)

# there is one year where no 13 y.o. were harvested because the harvest #s were too low;
# remove this 
df$n_age13[df$n_age13 == 0] <- NA
df$prop_age13[is.na(df$n_age13)==TRUE] <- NA
df$harvested_age13[is.na(df$n_age13)==TRUE] <- NA
df$harvested_total[is.na(df$n_age13)==TRUE] <- NA

############################################################
### ------------------ DATA WRITING -------------------- ###
############################################################

stop('[1.3.mgmt_fecundity.R] \n
All required matrices constructed. Code stopped to prevent overwriting data.\n
Continue running code beyond this line to overwrite data exports.')

write.csv(df, 'data/intermediate/fecundity.csv')
