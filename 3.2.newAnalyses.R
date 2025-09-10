### CJS Model for Elk Data
### Last updated: Sept. 9, 2025
### Contact: xprockox@gmail.com

### --------------- PACKAGES ---------------- ###
library(dplyr)
library(lubridate)
library(tidyr)
library(MCMCvis)
library(ggplot2)

### --------------- DATA IMPORT ---------------- ###
df <- read.csv('data/elk_survival_testData.csv')
df <- df[which(complete.cases(df)==TRUE),]

### --------------- DATA MANAGEMENT ---------------- ###

# clean df
df_clean <- df %>%
  mutate(
    BirthDate = ymd(BirthDate),
    Capture.Date = ymd(Capture.Date),
    Last.Date.Alive = ymd(Last.Date.Alive),
    Last.Date.Status = trimws(Last.Date.Status),
    BirthYear = year(BirthDate)
  ) %>%
  filter(!is.na(Capture.Date) & !is.na(Last.Date.Alive))  # Remove incomplete entries

# define study period and sample size
min_year <- min(df_clean$BirthYear, na.rm = TRUE)
max_year <- max(year(df_clean$Last.Date.Alive), na.rm = TRUE)
years <- min_year:max_year
n_years <- length(years)
n_indiv <- nrow(df_clean)

# intialize matrices
y <- matrix(0, nrow = n_indiv, ncol = n_years)
z <- matrix(NA, nrow = n_indiv, ncol = n_years)
colnames(y) <- colnames(z) <- years
rownames(y) <- rownames(z) <- df_clean$ID

# fill matrices
for (i in 1:n_indiv) {
  elk <- df_clean[i, ]
  birth_year <- elk$BirthYear
  capture_year <- year(elk$Capture.Date)
  last_year <- year(elk$Last.Date.Alive)
  status <- elk$Last.Date.Status
  
  # find column indices
  birth_idx <- which(years == birth_year)
  capture_idx <- which(years == capture_year)
  last_idx <- which(years == last_year)
  
  # latent state: assume alive (1) from birth until year of last confirmed alive
  if (!is.na(birth_idx) && !is.na(last_idx)) {
    z[i, birth_idx:last_idx] <- 1
  }
  
  # if the individual was confirmed dead, mark all future years as 0
  if (status == "Dead" && last_idx < n_years) {
    z[i, (last_idx + 1):n_years] <- 0
  }
  
  # detection matrix:
  if (status == "Dead" && capture_idx < last_idx) {
    y[i, capture_idx:(last_idx - 1)] <- 1
    y[i, last_idx] <- 0  # likely died that year
  } else if (status == "Live" && capture_idx <= last_idx) {
    y[i, capture_idx:last_idx] <- 1  # censored after this
  } else if (capture_idx == last_idx) {
    y[i, capture_idx] <- 1  # one-year animal
  }
}

# visual check of matrices
image(y, main = "Detection Matrix (y)", col = c("white", "black"))
image(z, main = "Latent State Matrix (z)", col = c("grey", "black"))

# vector of when the elk were first seen
first_seen <- apply(y, 1, function(row) {
  first <- which(row == 1)[1]
  if (is.na(first)) return(NA) else return(first)
})

#########################################################################
### ------------------------ MODEL ONE ------------------------------ ###
#########################################################################

### Model one is the simplest version of a CJS model, where elk survival
### is assumed to be constant across age-classes and throughout time.

### --------------- NIMBLE CODE ---------------- ###

elk_survival_constant <- nimbleCode({
  # Priors
  phi ~ dunif(0, 1) # Survival probability
  p ~ dunif(0, 1) # Detection probability
  
  for (i in 1:N) {
    z[i, first_seen[i]] <- 1 # All animals alive at first seen
    
    for (t in (first_seen[i] + 1):n_years) {
      # Only proceed if we’re not past their observed window
      z[i, t] ~ dbern(z[i, t - 1] * phi)
    }
    
    for (t in first_seen[i]:n_years) {
      y[i, t] ~ dbern(p * z[i, t])
    }
  }
})

### --------------- CONSTANTS AND SPECS ---------------- ###

data <- list(
  y = y
)

constants <- list(
  N = nrow(y),
  n_years = ncol(y),
  first_seen = first_seen
)

# create z init matrix
z_init <- z
z_init[is.na(z_init)] <- 1  # assume alive unless evidence says otherwise

inits <- list(
  phi = runif(1, 0.8, 0.99),
  p = runif(1, 0.6, 0.95),
  z = z_init
)

params <- c("phi", "p")

### --------------- RUN MODEL ---------------- ###

elk_model_constant <- nimbleMCMC(
  code = elk_survival_constant,
  constants = constants,
  data = data,
  inits = inits,
  monitors = params,
  nchains = 3,
  niter = 10000,
  nburnin = 3000,
  thin = 1,
  summary = TRUE
)

elk_model_constant

### --------------- ASSESS MODEL CONVERGENCE ---------------- ###

# traceplots
MCMCtrace(object = elk_model_constant$samples,
          pdf = FALSE,
          ind = TRUE,
          Rhat = TRUE,
          n.eff = TRUE,
          params = 'all')

# r-hats
MCMCsummary(elk_model_constant$samples, params = 'all')

#########################################################################
### ------------------------ MODEL TWO ------------------------------ ###
#########################################################################

### Model two is also not a stage-specific model, but this model does 
### incorporate time-varying survival.

### --------------- NIMBLE CODE ---------------- ###

elk_survival_temporal <- nimbleCode({
  # Priors for detection
  p ~ dunif(0, 1)
  
  # Priors for time-varying survival
  for (t in 1:(n_years - 1)) {
    phi[t] ~ dunif(0, 1)
  }
  
  for (i in 1:N) {
    z[i, first_seen[i]] <- 1  # Alive at first observation
    
    for (t in (first_seen[i] + 1):n_years) {
      z[i, t] ~ dbern(z[i, t - 1] * phi[t - 1])  # Survival depends on year t-1
    }
    
    for (t in first_seen[i]:n_years) {
      y[i, t] ~ dbern(p * z[i, t])
    }
  }
})

### --------------- CONSTANTS AND SPECS ---------------- ###

# Data
data <- list(
  y = y
)

# Constants
constants <- list(
  N = nrow(y),
  n_years = ncol(y),
  first_seen = first_seen
)

# Initial z matrix (assume alive)
z_init <- z
z_init[is.na(z_init)] <- 1

# Initial values
inits <- list(
  phi = runif(ncol(y) - 1, 0.8, 0.99),
  p = runif(1, 0.6, 0.95),
  z = z_init
)

# Parameters to monitor
params <- c("phi", "p")

### --------------- RUN MODEL ---------------- ###

elk_model_temporal <- nimbleMCMC(
  code = elk_survival_temporal,
  constants = constants,
  data = data,
  inits = inits,
  monitors = params,
  nchains = 3,
  niter = 10000,
  nburnin = 3000,
  thin = 1,
  summary = TRUE
)

elk_model_temporal

### --------------- ASSESS MODEL CONVERGENCE ---------------- ###

# check traceplots
MCMCtrace(object = elk_model_temporal$samples,
          pdf = FALSE,
          ind = TRUE,
          Rhat = TRUE,
          n.eff = TRUE,
          params = 'all')

# summary stats
MCMCsummary(elk_model_temporal$samples, params = 'all')

# save summary stats
elk_model_temporal_results <- MCMCsummary(elk_model_temporal$samples, params = 'all')

# extract phi estimates
elk_temporal_phi <- elk_model_temporal_results %>%
  filter(grepl("phi\\[", rownames(elk_model_temporal_results))) %>%
  mutate(Year = 1983:2024)  

# plot temporal trend
ggplot(elk_temporal_phi) + 
  geom_ribbon(aes(x = Year, ymin = `2.5%`, ymax = `97.5%`), fill = "grey90") +
  geom_path(aes(x = Year, y = mean), color = "red", linewidth = 1) +
  scale_y_continuous(name = "Estimated Survival (φ)") +
  labs(x = "Year", title = "Time-varying survival estimates") +
  theme_bw()

#########################################################################
### ----------------------- MODEL THREE ----------------------------- ###
#########################################################################

### Model three expands the temporally varying survival model to include
### unique survival values for each "class" or "stage" of elk. The stages are:
### 1: 0-1 year old
### 2: 2-13 years old
### 3: 14+ years old

### First we need to build a matrix that shows these stages.

# copy z matrix and change all non-1 values to NA 
# (1 is when the animal was known to be alive)
age_class <- z
age_class[age_class != 1] <- NA

# loop over rows and assign age classes to 1s based on relative position
for (i in 1:nrow(age_class)) {
  alive_years <- which(age_class[i, ] == 1)
  
  n1 <- min(2, length(alive_years)) # First 2 years
  n2 <- min(11, length(alive_years) - n1) # Next 11 years
  n3 <- length(alive_years) - n1 - n2 # Remaining years
  
  age_class[i, alive_years[1:n1]] <- 1
  if (n2 > 0) age_class[i, alive_years[(n1 + 1):(n1 + n2)]] <- 2
  if (n3 > 0) age_class[i, alive_years[(n1 + n2 + 1):length(alive_years)]] <- 3
}

### --------------- NIMBLE CODE ---------------- ###

elk_survival_stage_specific <- nimbleCode({
  
  # Priors for survival by age class
  for (k in 1:3) {
    phi[k] ~ dunif(0, 1)
  }
  
  # Detection probability
  p ~ dunif(0, 1)
  
  for (i in 1:N) {
    # Latent alive state
    z[i] ~ dbern(phi[ age_class[i] ])
    
    # Observation model
    y[i] ~ dbern(p * z[i])
  }
})

### --------------- CONSTANTS AND SPECS ---------------- ###

# data
data <- list(
  y = y
)

# constants
constants <- list(
  N = nrow(y),
  n_years = ncol(y),
  first_seen = first_seen,
  age_class = age_class
)

# initial values
z_init <- z
z_init[is.na(z_init)] <- 1

inits <- list(
  phi = runif(3, 0.6, 0.95),
  p = runif(1, 0.6, 0.95),
  z = z_init
)

# params to monitor
params <- c("phi", "p")

### --------------- RUN MODEL ---------------- ###

elk_model_stage_specific <- nimbleMCMC(
  code = elk_survival_stage_specific,
  constants = constants,
  data = data,
  inits = inits,
  monitors = params,
  nchains = 3,
  niter = 10000,
  nburnin = 3000,
  thin = 1,
  summary = TRUE
)

elk_model_stage_specific

### --------------- ASSESS MODEL CONVERGENCE ---------------- ###

# check traceplots
MCMCtrace(object = elk_model_temporal$samples,
          pdf = FALSE,
          ind = TRUE,
          Rhat = TRUE,
          n.eff = TRUE,
          params = 'all')

# summary stats
MCMCsummary(elk_model_temporal$samples, params = 'all')

# save summary stats
elk_model_temporal_results <- MCMCsummary(elk_model_temporal$samples, params = 'all')

# extract phi estimates
elk_temporal_phi <- elk_model_temporal_results %>%
  filter(grepl("phi\\[", rownames(elk_model_temporal_results))) %>%
  mutate(Year = 1983:2024)  

# plot temporal trend
ggplot(elk_temporal_phi) + 
  geom_ribbon(aes(x = Year, ymin = `2.5%`, ymax = `97.5%`), fill = "grey90") +
  geom_path(aes(x = Year, y = mean), color = "red", linewidth = 1) +
  scale_y_continuous(name = "Estimated Survival (φ)") +
  labs(x = "Year", title = "Time-varying survival estimates") +
  theme_bw()