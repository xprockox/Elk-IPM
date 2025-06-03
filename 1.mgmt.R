### Elk Integrated Population Model
### Data management script
### Last updated: May 29, 2025
### Contact: xprockox@gmail.com

############################################################################################
### packages
library(tidyverse)

############################################################################################
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

############################################################################################
### visualization

# adult survival by age
ggplot(ad_surv, aes(x = Age, y = Mean, color = Parameter)) +
  geom_path(size = 1) +
  geom_point(size = 2) +
  geom_errorbar(aes(ymin = Lower, ymax = Upper), width = 0.2) +
  labs(title = "Age-Specific Adult Elk Survival by Parameter",
       y = "Survival Probability",
       x = "Age",
       color = "Parameter") +
  theme_minimal()

# adult survival by demographic class - prime (2-13 y.o.) vs. old (14+ y.o.)
ann_ad_surv <- ann_ad_surv %>%
  mutate(vjust_pos = ifelse(Cow_category == "Prime", -1.1, 1.5))

ggplot(ann_ad_surv, aes(x = Year, y = Survival, color = Cow_category)) +
  geom_line(size = 1) +
  geom_point(size = 2) +
  geom_errorbar(aes(ymin = Survival - Survival.se, ymax = Survival + Survival.se),
                width = 0.2, color = "black") +
  geom_text(aes(label = paste0("n = ", N), vjust = vjust_pos), hjust = -0.3, size = 3, show.legend = FALSE) +
  labs(title = "Annual Adult Elk Survival",
       y = "Survival Probability",
       x = "Year") +
  theme_minimal()

# calf survival and fecundity
calf_plot_df <- calf_haz %>%
  select(Year, S.mu, FemBirthPerCow) %>%
  pivot_longer(cols = c(S.mu, FemBirthPerCow),
               names_to = "Metric",
               values_to = "Value") %>%
  mutate(Metric = recode(Metric,
                         S.mu = "Calf Survival",
                         FemBirthPerCow = "Fecundity (female births per cow)"))

ggplot(calf_plot_df, aes(x = Year, y = Value, color = Metric, fill = Metric)) +
  geom_line(size = 1.2) +
  geom_point(size = 3, shape = 21, stroke = 0.4, color = "black") +
  scale_color_manual(values = c("Calf Survival" = "#1b9e77",
                                "Fecundity (female births per cow)" = "#d95f02")) +
  scale_fill_manual(values = c("Calf Survival" = "#1b9e77",
                               "Fecundity (female births per cow)" = "#d95f02")) +
  labs(title = "Annual Calf Survival and Fecundity",
       y = "Rate",
       x = "Year",
       color = NULL,
       fill = NULL) +
  theme_minimal(base_size = 14) +
  theme(legend.position = "top",
        legend.box = "horizontal",
        panel.grid.minor = element_blank())

# calf survival with CIs
ggplot(calf_surv, aes(x = Year, y = P0)) +
  geom_line(size = 1) +
  geom_point(size = 2) +
  geom_errorbar(aes(ymin = Lower, ymax = Upper),
                width = 0.2, color = "black") +
  labs(title = "Annual Calf Survival",
       y = "Survival Probability",
       x = "Year") +
  theme_minimal()

# harvest
harvest_long <- harvest %>%
  pivot_longer(cols = starts_with("Y"),
               names_to = "Year",
               names_prefix = "Y",
               values_to = "Harvest") %>%
  mutate(Year = as.integer(Year))

ggplot(harvest_long, aes(x = Year, y = Harvest, color = factor(Age))) +
  geom_line(size = 0.6) +
  labs(title = "Elk Harvest by Age Over Time",
       x = "Year",
       y = "Number Harvested",
       color = "Age (Years)") +
  theme_minimal(base_size = 14)

############################################################################################
### formatting

# Base year df
df <- data.frame(Year = 1995:2025)

# Prime adult survival
prime_ann_ad_surv <- ann_ad_surv %>%
  filter(Cow_category == "Prime") %>%
  select(Cow_category, Year,
         Wolf.mort, Wolf.se, Other.mort, Other.se,
         Survival, Survival.se) %>%
  rename(Stage = Cow_category, Year = Year)

names(prime_ann_ad_surv)[-c(1:2)] <- paste0(names(prime_ann_ad_surv)[-c(1:2)], "_prime")

# old adult survival data formatting

old_ann_ad_surv <- ann_ad_surv %>%
  filter(Cow_category == "Old") %>%
  select(Cow_category, Year,
         Wolf.mort, Wolf.se, Other.mort, Other.se,
         Survival, Survival.se) %>%
  rename(Stage = Cow_category, Year = Year)

names(old_ann_ad_surv)[-c(1:2)] <- paste0(names(old_ann_ad_surv)[-c(1:2)], "_old")

# calf survival data formatting
calf_surv_clean <- calf_surv %>%
  rename(Survival_calf = P0) %>%
  select(Year, Survival_calf)

# binding all survival dataframes together
df <- left_join(df, calf_surv_clean, by='Year')
df <- left_join(df, prime_ann_ad_surv, by='Year')
df <- left_join(df, old_ann_ad_surv, by='Year')

# clean up
df <- df %>%
  select(-c(Stage.x, Stage.y))

# combine with other data 
df <- left_join(df, obs_elk_N)
df <- left_join(df, obs_elk_inNout)

modeling_df <- df %>%
  select(Year, 
         Survival_calf, Survival_prime, Survival_old, 
         Total_Elk, Percent.female, Total_Elk_Female, Total_harvest, Inside_Elk, Outside_Elk)

# include fecundity
fec <- calf_haz %>%
  select(Year, FemBirthPerCow)

modeling_df <- left_join(modeling_df, fec)

# let's include information about stage distributions

stage_dists95 <- N1995_dems %>%
  select(Age, Age.dist, N.cows.posthunt.corrected) %>%
  mutate(Year = 1995)

stage_dists95_calves <- stage_dists95[stage_dists95$Age==1,]
stage_dists95_prime <- stage_dists95[stage_dists95$Age %in% c(2:13),]
stage_dists95_old <- stage_dists95[stage_dists95$Age >= 14,]

modeling_df$Percent.N.calves <- NA
modeling_df$Percent.N.prime <- NA
modeling_df$Percent.N.old <- NA

modeling_df[modeling_df$Year==1995, 'Percent.N.calves'] <- stage_dists95_calves$Age.dist
modeling_df[modeling_df$Year==1995, 'Percent.N.prime'] <- sum(stage_dists95_prime$Age.dist)
modeling_df[modeling_df$Year==1995, 'Percent.N.old'] <- sum(stage_dists95_old$Age.dist)

############################################################################################
# clean up workspace and write modeling dataframe to .csv
write_csv(modeling_df, 'data/intermediate/modeling_df.csv')

rm(ad_surv, all_elk, ann_ad_surv, calf_haz, calf_plot_df, calf_surv, calf_surv_clean,
   df, fec, harvest, harvest_long, N1995_dems, obs_elk_CIs, obs_elk_inNout, obs_elk_N, 
   old_ann_ad_surv, preg, prime_ann_ad_surv, stage_dists95, stage_dists95_calves,
   stage_dists95_old, stage_dists95_prime)
