

elk_temporal_phi <- read.csv('data/intermediate/temp_phi_noStages.csv')
class_df <- read.csv('data/master/Classification_1995_to_2015.csv')
counts <- read.csv('data/master/corrected_elk counts_09Jul2023.csv')

elk_temporal_phi <- elk_temporal_phi %>%
  rename(cow_phi = mean,
         cow_phi_lower95ci = `X2.5.`,
         cow_phi_upper95ci = `X97.5.`) %>%
  select(-X)

counts <- counts %>%
  rename(Year = `winter..Jan..`,
         total_n = mean,
         total_n_lower95ci = `lwr.CL`,
         total_n_upper95ci = `upr.CL`) %>%
  select(Year, total_n, total_n_lower95ci, total_n_upper95ci)

class_df <- class_df %>%
  mutate(percent_cows = Adult.Females/Total.Classified,
         percent_calves = Calves/Total.Classified,
         percent_bulls = Total.Bulls/Total.Classified,
         percent_spikes = Spikes/Total.Classified) %>%
  group_by(Year) %>%
  mutate(percent_cows = mean(percent_cows),
         percent_calves = mean(percent_calves),
         percent_bulls = mean(percent_bulls),
         percent_spikes = mean(percent_spikes),
         Total.Classified = sum(Total.Classified)) %>%
  select(Year, Total.Classified, percent_cows, percent_calves, percent_bulls, percent_spikes) %>%
  ungroup()

class_df <- unique(class_df)

df <- left_join(counts, class_df)
df <- left_join(df, elk_temporal_phi)
