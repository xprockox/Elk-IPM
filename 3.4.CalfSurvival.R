### Calf Survival and Reproduction
### Last updated: Oct. 20, 2025
### Contact: xprockox@gmail.com

### --------------- PACKAGES ---------------- ###

library(dplyr)
library(lubridate)

### --------------- DATA IMPORT ---------------- ###

classification <- read.csv('data/master/classification_2025-10-17.csv')
counts <- read.csv('data/master/corrected_elk counts_09Jul2023.csv')
preg <- read.csv('data/master/pregnancy_collars.csv')
preg_harvest <- read.csv('data/master/pregnancy_harvest.csv')

### --------------- DATA MANAGEMENT ---------------- ###

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
         n_LCI = lwr.CL,
         n_UCI = upr.CL) %>%
  select(year, n_total, n_LCI, n_UCI)

df <- left_join(classification, counts)

df$calf_cow_ratio <- df$calf_per100cow / 100
df$spike_cow_ratio <- df$spike_per100cow / 100
df$btb_cow_ratio <- df$btb_per100cow/ 100
df$bull_cow_ratio <- df$bull_per100cow / 100

df$total_elk <- as.numeric(df$total_elk)

df$percent_cows <- df$cows / df$total_elk

df$n_cows <- round(df$n_total * df$percent_cows)
df$n_calves <- round(df$n_total * df$percent_cows * df$calf_cow_ratio)
df$n_spikes <- round(df$n_total * df$percent_cows * df$spike_cow_ratio)
df$n_btb <- round(df$n_total * df$percent_cows * df$btb_cow_ratio)
df$n_bulls <- round(df$n_total * df$percent_cows * df$bull_cow_ratio)

# incorporating park preg data from collared indivs
preg <- preg %>%
  mutate(capt_year = year(as.POSIXct(strptime(Capture.Date, format = '%d-%b-%y')))) %>%
  group_by(capt_year) %>%
  count(Pregnant)

preg <- preg[preg$Pregnant=='yes' | preg$Pregnant=='no',]

preg <- preg %>%
  group_by(capt_year) %>%
  mutate(num_capt = sum(n),
         percent_preg = n / num_capt) %>%
  rename(year = capt_year)

preg <- preg[preg$Pregnant=='yes',]

df <- left_join(df, preg)

df$expected_calves <- df$n_cows * df$percent_preg

df$calf_surv <- df$n_calves / df$expected_calves

# incorporating FWP preg data from harvested elk during antlerless harvest years (1997-2009)
preg_harvest <- preg_harvest %>%
  select(ID, harvestyear, ageatharvest, pregnant_code, winter_end) %>%
  group_by(harvestyear) %>%
  count(pregnant_code) %>%
  mutate(num_harv = sum(n),
         percent_preg = n / num_harv) %>%
  rename(year = harvestyear)

preg_harvest <- preg_harvest[preg_harvest$pregnant_code==1,]

df <- left_join(df, preg_harvest, by='year')
