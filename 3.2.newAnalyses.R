### Elk IPM
### Attempt 2
### Last updated: August 27, 2025
### Contact: xprockox@gmail.com

########## ----------------------- PACKAGES ----------------------- ############

library(nimble)
library(dplyr)
library(lubridate)

########## ---------------------- DATA IMPORT --------------------- ############

# dataframe with three columns: elk_ID, birth_year, and death_year
elk_df <- read.csv('data/elk_survival_testData.csv')


########## -------------------- DATA MANAGEMENT ------------------- ############

# create required columns
elk_df <- elk_df %>%
  mutate(death_year = year(as.Date(Last.Date.Alive))) %>%
  rename(elk_ID = ID,
         birth_year = BirthYear)

# drop rows with incomplete information
elk_df <- elk_df[which(complete.cases(elk_df)==TRUE),]

# Ensure elk_ID is a factor or character (doesn’t affect the model, but useful)
elk_df$elk_ID <- as.character(elk_df$elk_ID)

# Create a column showing whether the elk is currently alive
elk_df$Alive <- ifelse(elk_df$Last.Date.Alive == '2024-02-18', "Yes", "No")

# Step 1: Determine full time span
min_year <- min(elk_df$birth_year)
max_year <- max(elk_df$death_year)
all_years <- min_year:max_year
T <- length(all_years)

# Step 2: Create mapping of years to columns
year_index <- setNames(1:T, all_years)  # e.g., year_index[["2005"]] = column 5

# Step 3: Create obs matrix
N <- nrow(elk_df)
obs <- matrix(NA, nrow = N, ncol = T)  # rows = elk, columns = years

for (i in 1:N) {
  by <- elk_df$birth_year[i]
  dy <- elk_df$death_year[i]
  first_idx <- year_index[as.character(by)]
  death_idx <- year_index[as.character(dy)]
  
  # Elk is alive from birth year through year before death
  obs[i, first_idx:(death_idx - 1)] <- 1
  
  # Elk is dead at death year
  obs[i, death_idx] <- 0
}

# Step 4: Define other required variables
birth_year_vec <- elk_df$birth_year
first_vec <- year_index[as.character(elk_df$birth_year)]
last_vec <- year_index[as.character(elk_df$death_year)]

# Resulting variables:
# - obs: N x T matrix of 1 (alive), 0 (dead), NA (not yet born or already dead)
# - birth_year_vec: vector of birth years
# - first_vec: numeric index of first year in matrix per elk
# - last_vec: numeric index of last year in matrix per elk



survival_code <- nimbleCode({
  # Priors for survival probabilities (on logit scale)
  logit_phi[1] ~ dnorm(0, sd = 1.5)  # Age 0–2
  logit_phi[2] ~ dnorm(0, sd = 1.5)  # Age 2–14
  logit_phi[3] ~ dnorm(0, sd = 1.5)  # Age 14+
  
  for (i in 1:N) {
    alive[i, first[i]] <- 1  # alive at first year (birth)
    
    for (t in (first[i] + 1):last[i]) {
      age <- t - birth_year[i]  # age in year t
      age_class <- 1 + (age > 2) + (age > 14)
      logit(p[i, t]) <- logit_phi[age_class]
      
      alive[i, t] ~ dbern(alive[i, t - 1] * p[i, t])  # survival process
    }
    
    # Observation: 0 = died in that year, 1 = lived to next year
    for (t in first[i]:(last[i] - 1)) {
      obs[i, t] ~ dbern(alive[i, t])
    }
  }
})


# Constants
constants <- list(N = n_individuals, first = first_year_vec, last = last_year_vec)

# Data
data <- list(
  obs = obs_matrix,              # 1 = observed alive, 0 = known dead
  birth_year = birth_year_vec    # e.g., 1997, 1998...
)

# Inits
inits <- list(
  logit_phi = rep(0, 3),
  alive = matrix(1, nrow = N, ncol = T_max)  # start assuming alive
)