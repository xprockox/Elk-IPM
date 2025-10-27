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
preg <- read.csv('data/master/elk_pregnancy_collars_2025-07-18.csv')
preg_harvest <- read.csv('data/master/elk_pregnancy_harvest_2015-03-24.csv')

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

# merges the two data sources of pregnancy (capture and harvest), 
# including harvest estimate only when the capture estimate does not exist
df <- df %>%
  rename(total_elk_classified = total_elk) %>%
  mutate(percent_preg = round(coalesce(percent_preg.x, percent_preg.y), 2),
         expected_calves = round(n_cows * percent_preg, 2),
         calf_surv = round(n_calves/expected_calves, 2)) %>%
  select(-Pregnant, -survey_date, -Method, 
         -total_elk_classified, -cows, -calves,
         -calf_per100cow, -spike_per100cow, -btb_per100cow, -bull_per100cow,
         -n.x, -num_capt, -pregnant_code, 
         -n.y, -num_harv, -percent_preg.x, -percent_preg.y)

############################################################
### ------------------ DATA WRITING -------------------- ###
############################################################

stop('[1.3.mgmt_productivity.R] \n
All required matrices constructed. Code stopped to prevent overwriting data.\n
Continue running code beyond this line to overwrite data exports.')

write.csv(df, 'data/intermediate/productivity.csv')
