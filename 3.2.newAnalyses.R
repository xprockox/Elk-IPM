### CJS Model for Elk Data
### Last updated: Sept. 16, 2025
### Contact: xprockox@gmail.com

### --------------- PACKAGES ---------------- ###
library(dplyr)
library(lubridate)
library(tidyr)
library(MCMCvis)
library(ggplot2)
library(nimble)

### --------------- DATA IMPORT ---------------- ###
df <- read.csv('data/elk_survival_testData.csv')
df <- df[which(complete.cases(df)==TRUE),] # drop rows with any NAs

vhf <- read.csv('data/vhf.csv')
gps <- readRDS('data/res.rds')

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

# clean VHF and GPS data
vhf <- vhf %>%
  mutate(
    Date = as.POSIXct(Date, format = "%m/%d/%Y"),
    Signal.Heard = case_when(
      is.na(Signal.Heard) ~ "YES",
      Signal.Heard == "" ~ "YES",
      Signal.Heard %in% c("Yes", "YES") ~ "YES",
      Signal.Heard %in% c("No", "NO") ~ "NO",
      TRUE ~ as.character(Signal.Heard)
    ),
    Visual. = case_when(
      Visual. %in% c("YES", "Yes", "yes", "Y") ~ "YES",
      TRUE ~ "NO"
    ),
    Year = year(Date)
  ) %>%
  filter(Signal.Heard == "YES" | Visual. == "YES") %>%
  rename(Visual = Visual.)

gps <- gps %>%
  mutate(
    dt = as.POSIXct(dt, format = "%Y-%m-%d %H:%M:%S"),
  ) %>%
  left_join(
    df_clean %>% select(ID, Last.Date.Alive), 
    by = "ID"
  ) %>%
  mutate(
    Year = year(dt),
    LastYear = year(as.POSIXct(Last.Date.Alive, format = "%Y-%m-%d"))
  ) %>%
  filter(Year <= LastYear)  # Remove GPS points after animal's death

# define study period and sample size
min_year <- min(df_clean$BirthYear, na.rm = TRUE)
max_year <- max(year(df_clean$Last.Date.Alive), na.rm = TRUE)
years <- min_year:max_year
n_years <- length(years)
n_indiv <- nrow(df_clean)

### --------------- MATRIX CONSTRUCTION (OBS. & LATENT STATES) ---------------- ###

# Initialize matrices
y <- matrix(0, nrow = n_indiv, ncol = n_years) # observations
z <- matrix(NA, nrow = n_indiv, ncol = n_years) # latent (true) states
colnames(y) <- colnames(z) <- years
rownames(y) <- rownames(z) <- df_clean$ID

# Fill in matrices
for (i in 1:n_indiv) {
  elk <- df_clean[i, ]
  elk_id <- elk$ID
  
  birth_year <- elk$BirthYear
  capture_year <- year(elk$Capture.Date)
  last_year <- year(elk$Last.Date.Alive)
  status <- elk$Last.Date.Status
  
  # Column indices
  birth_idx <- which(years == birth_year)
  capture_idx <- which(years == capture_year)
  last_idx <- which(years == last_year)
  
  # Latent state: assume alive from birth through year of last known alive
  if (!is.na(birth_idx) && !is.na(last_idx)) {
    z[i, birth_idx:last_idx] <- 1
  }
  
  # If known dead, mark all years after death as 0
  if (status == "Dead" && last_idx < n_years) {
    z[i, (last_idx + 1):n_years] <- 0
  }
  
  # Observation matrix (y) constructed based on VHF and GPS data
  
  # Extract observed years from VHF
  vhf_years <- vhf %>%
    filter(ID == elk_id) %>%
    pull(Year) %>%
    unique()
  
  # Extract observed years from GPS
  gps_years <- gps %>%
    filter(ID == elk_id) %>%
    pull(Year) %>%
    unique()
  
  # Combine and deduplicate
  observed_years <- sort(unique(c(vhf_years, gps_years)))
  
  # Fill y matrix
  for (yr in observed_years) {
    if (yr %in% years) {
      y[i, as.character(yr)] <- 1
    }
  }
}

# The above process does not necessarily mark the year of capture as a detection, 
# e.g. if the collar failed immediately and provided no GPS or VHF data,
# so we can force detections for all individuals during their year of capture:
df_clean$CaptureYear <- lubridate::year(df_clean$Capture.Date)

for (i in 1:nrow(df_clean)) {
  cap_year <- df_clean$CaptureYear[i]
  cap_col <- which(colnames(y) == cap_year)
  if (length(cap_col) == 1) {
    y[i, cap_col] <- 1
  }
}

# Now if the latent state (true state) is dead, there can't be an observation.
# Mask impossible observations in y: set to NA if z = 0
for (i in 1:n_indiv) {
  for (t in 1:n_years) {
    if (!is.na(z[i, t]) && z[i, t] == 0 && y[i, t] == 1) {
      message(paste("Setting y[", i, ",", t, "] to NA due to z = 0 but y = 1"))
      y[i, t] <- NA  # mask inconsistent observations
    }
    
    # Optional: also mask y after last known alive (if not already done)
    if (!is.na(z[i, t]) && z[i, t] == 0) {
      y[i, t] <- NA  # cannot be observed if truly dead
    }
  }
}

# And add a final check:

# If z = 0, we should not see an observation (mask it)
y[z == 0] <- 0

# Ensure that wherever y = 1, z = 1
z[y == 1] <- 1

### --------------- VISUALIZE MATRICES ---------------- ###

### first observations:
# Flip y so individual 1 is at the top
y_flip <- y[nrow(y):1, ]

# Keep only columns for years >= 2000
cols_to_keep <- which(as.numeric(colnames(y_flip)) >= 2000)
y_flip <- y_flip[, cols_to_keep]

# Plot detection matrix
image(
  x = 1:ncol(y_flip),
  y = 1:nrow(y_flip),
  z = t(y_flip),  # transpose for image()
  col = c("white", "black"),
  axes = FALSE,
  xlab = "Year",
  ylab = "Individual",
  main = "Detection Matrix (y)"
)

# Add axis labels
axis(1,
     at = 1:ncol(y_flip),
     labels = colnames(y_flip),
     las = 2,
     cex.axis = 0.7)

axis(2,
     at = 1:nrow(y_flip),
     labels = rev(rownames(y)),  # reversed to match y_flip
     las = 1,
     cex.axis = 0.4)

### then latent states:
# Flip z so that individual 1 is at the top
z_flip <- z[nrow(z):1, ]

# Set up plot
image(
  x = 1:ncol(z_flip),
  y = 1:nrow(z_flip),
  z = t(z_flip),  # transpose for correct orientation
  col = c("white", "gray", "black"),  # NA = white, 0 = gray, 1 = black
  breaks = c(-0.1, 0.1, 0.9, 1.1),    # assign color by value (NA, 0, 1)
  axes = FALSE,
  xlab = "Year",
  ylab = "Individual",
  main = "Latent State Matrix (z)"
)

# Add axis labels
axis(1, at = 1:ncol(z_flip), labels = colnames(z_flip), las = 2, cex.axis = 0.7)
axis(2, at = 1:nrow(z_flip), labels = rev(rownames(z_flip)), las = 1, cex.axis = 0.4)

rm(z_flip, y_flip)

#########################################################################
### ------------------------ MODEL ONE ------------------------------ ###
#########################################################################

### Model one is the simplest version of a CJS model, where elk survival
### is assumed to be constant across age-classes and throughout time.

# for this, we need a vector of when each elk was first seen
first_seen <- apply(y, 1, function(row) {
  first <- which(row == 1)[1]
  if (is.na(first)) return(n_years) else return(first)
})
first_seen <- as.integer(first_seen)

### --------------- NIMBLE CODE ---------------- ###

elk_survival_constant <- nimbleCode({
  phi ~ dunif(0, 1)
  p ~ dunif(0, 1)
  
  for (i in 1:N) {
    z[i, first_seen[i]] <- 1
    
    # Only define z beyond first_seen[i] if possible
    if (first_seen[i] < n_years) {
      for (t in (first_seen[i] + 1):n_years) {
        z[i, t] ~ dbern(z[i, t - 1] * phi)
      }
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

elk_model_constant$summary$all.chains

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
  # Priors
  p ~ dunif(0, 1)
  for (t in 1:(n_years - 1)) {
    phi[t] ~ dunif(0, 1)
  }
  
  for (i in 1:N) {
    z[i, first_seen[i]] <- 1  # alive at first observation
    
    for (t in (first_seen[i] + 1):n_years) {
      z[i, t] ~ dbern(z[i, t - 1] * phi[t - 1] + 1e-10 * (1 - z[i, t - 1]))
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

# Initial z matrix (assume not alive)
z_init <- z
z_init[is.na(z_init)] <- 0

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

elk_model_temporal$summary$all.chains

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
  scale_x_continuous(name = "Year", limits = c(1999, 2024)) +
  labs(x = "Year", title = "Time-varying survival estimates") +
  theme_bw()

#########################################################################
### ----------------------- MODEL THREE ----------------------------- ###
#########################################################################

### Model three expands the temporally varying survival model to include
### unique survival values for each "class" or "stage" of elk. The stages are:
### 1: 0-13 years old
### 2: 14+ years old

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

# also need last seen to identify the stage it was in when it died
years_vector <- as.numeric(colnames(y))
last_seen <- year(as.Date(df$Last.Date.Alive))
last_seen <- match(last_seen, years_vector)

### --------------- NIMBLE CODE ---------------- ###

elk_survival_stage_specific <- nimbleCode({
  
  # Priors for survival by age class (stage)
  phi_1 ~ dunif(0, 1)  # Calf
  phi_2 ~ dunif(0, 1)  # Young Adult
  phi_3 ~ dunif(0, 1)  # Old Adult
  
  # Prior for detection probability
  p ~ dunif(0, 1)
  
  # State process
  for (i in 1:N) {
    z[i, first_seen[i]] <- 1  # Animal is alive at first detection
    
    # for all timesteps after, true state is a bernoulli dist of survival, stage-specific
    for (t in 1:n_years) {
      if (t > first_seen[i] & t <= last_seen[i]) {
        z[i, t] ~ dbern(z[i, t - 1] *
                          (age_class[i, t] == 1) * phi_1 +
                          (age_class[i, t] == 2) * phi_2 +
                          (age_class[i, t] == 3) * phi_3)
      }
    }
    
    # Observation process
    for (t in 1:n_years) {
      if (t >= first_seen[i] & t <= last_seen[i]) {
        y[i, t] ~ dbern(p * z[i, t])
      }
    }
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
  last_seen = last_seen_index
)

# initial values
z_init <- z
z_init[is.na(z_init)] <- 1

inits <- list(
  phi_1 = runif(1, 0.6, 0.95),
  phi_2 = runif(1, 0.6, 0.95),
  phi_3 = runif(1, 0.6, 0.95),
  p = runif(1, 0.6, 0.95),
  z = z_init
)

# params to monitor
params <- c("phi_1", "phi_2", "phi_3", "p")

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
MCMCtrace(object = elk_model_stage_specific$samples,
          pdf = FALSE,
          ind = TRUE,
          Rhat = TRUE,
          n.eff = TRUE,
          params = 'all')

# summary stats
MCMCsummary(elk_model_stage_specific$samples, params = 'all')

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