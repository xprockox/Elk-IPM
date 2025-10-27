### Elk Abundance - Data management
### Last updated: Oct. 27, 2025
### Contact: xprockox@gmail.com

############################################################################################
### packages
library(tidyverse)

############################################################################################
### data loading 
abundances <- read.csv('data/master/corrected_elk counts_09Jul2023.csv')
classification <- read.csv('data/master/classification_2025-10-17.csv')

############################################################################################
### management 

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

############################################################################################
### write to .csv 
stop(
  'The following with overwrite data, are you sure you want to proceed?'
)

write.csv(abundances, 'data/intermediate/abundanceEstimates_stages.csv')
############################################################################################