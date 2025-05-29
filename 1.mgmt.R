### Elk Integrated Population Model
### Data management script
### Last updated: May 29, 2025
### Contact: xprockox@gmail.com

#######################################################################################################
### packages
library(tidyverse)

#######################################################################################################
### data loading
ad_surv <- read.csv('data/master/AdultSurvival_Long.csv')
all_elk <- read.csv('data/master/AllElk.csv')
ann_ad_surv <- read.csv('data/master/Annual_AdultSurvival.csv')
calf_haz <- read.csv('data/master/Calf_Hazards_05-07-25.csv')
calf_surv <- read.csv('data/master/CalfSurvival_SE.csv')
harvest <- read.csv('data/master/Harvest.csv')
N1995_dems <- read.csv('data/master/N1995_ElkCountsByAge.csv')
obs_elk_CIs <- read.csv('data/master/Observed_Elk_95CI.csv')
obs_elk_N <- read.csv('data/master/Observed_Elk_Abundance.csv')
obs_elk_inNout <- read.csv('data/master/Observed_Elk_InsideOutside.csv')
preg <- read.csv('data/master/Pregnancy.csv')
 