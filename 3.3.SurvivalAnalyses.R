### CJS Model for Elk Data
### Last updated: Sept. 16, 2025
### Contact: xprockox@gmail.com

### --------------- PACKAGES AND SET-UP ---------------- ###
library(dplyr)
library(lubridate)
library(tidyr)
library(tidyverse)
library(MCMCvis)
library(ggplot2)
library(nimble)

# the following line determines whether your latent state and observation
# matrices begin at the first year an elk was born (1982) or at the 
# first year an elk was captured (2000).

# the logic here is that there is no detection data from prior to capture,
# so estimates might be impacted by the mismatch between detection data and 
# latent state data (i.e. detection probability would be zero up until capture)

# if use_clip == TRUE, models are constructed using data from 2000:present
# if use_clip == FALSE, models are constructed using the full dataset (1982:present)
use_clip <- TRUE 

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
  filter(!is.na(Capture.Date) & !is.na(Last.Date.Alive)) %>% # Remove incomplete entries
  arrange(ID, Capture.Date)

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
  cap_idx <- which(years == cap_year)
  
  if (length(cap_idx) == 1) {
    y[i, cap_idx] <- 1
  }
  
  if (!is.na(cap_idx) && cap_idx > 1) {
    y[i, 1:(cap_idx - 1)] <- 0
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

# create a new latent state matrix that does not include information about births 
# (i.e. the individual is NA up until capture despite being alive)
# (i.e.i.e we're pretending we don't have cementum info)
z_clipped <- z  # Copy original latent state matrix

for (i in 1:n_indiv) {
  cap_year <- df_clean$CaptureYear[i]
  cap_idx <- which(years == cap_year)
  
  if (length(cap_idx) == 1) {
    z_clipped[i, 1:(cap_idx - 1)] <- NA  # Mask all years before capture
  }
}

# drop years before capture
z_clipped <- z_clipped[,19:ncol(z_clipped)]

# now do the same for the observation matrix
y_clipped <- y[,19:ncol(y)]

### --------------- VISUALIZE MATRICES ---------------- ###

### first observations:
# Flip y so individual 1 is at the top
y_flip <- y[nrow(y):1, ]

# Keep only columns for years >= 2000
# cols_to_keep <- which(as.numeric(colnames(y_flip)) >= 2000)
# y_flip <- y_flip[, cols_to_keep]

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
axis(2, at = 1:nrow(z_flip), labels = rownames(z_flip), las = 1, cex.axis = 0.4)

rm(z_flip, y_flip)

#
##
###
####
#####
######
### ------- then plot latent states CLIPPED:
######
#####
####
###
##
#

### start with y_clipped
# Flip y_clipped so that individual 1 is at the top
y_clipped_flip <- y_clipped[nrow(y_clipped):1, ]

# Set up plot
image(
  x = 1:ncol(y_clipped_flip),
  y = 1:nrow(y_clipped_flip),
  z = t(y_clipped_flip),  # transpose for correct orientation
  col = c("white", "gray", "black"),  # NA = white, 0 = gray, 1 = black
  breaks = c(-0.1, 0.1, 0.9, 1.1),    # assign color by value (NA, 0, 1)
  axes = FALSE,
  xlab = "Year",
  ylab = "Individual",
  main = "Detection Matrix (y)"
)

# Add axis labels
axis(1, at = 1:ncol(y_clipped_flip), labels = colnames(y_clipped_flip), las = 2, cex.axis = 0.7)
axis(2, at = 1:nrow(y_clipped_flip), labels = rownames(y_clipped_flip), las = 1, cex.axis = 0.4)

### then plot z_clipped
# Flip z_clipped so that individual 1 is at the top
z_clipped_flip <- z_clipped[nrow(z_clipped):1, ]

# Set up plot
image(
  x = 1:ncol(z_clipped_flip),
  y = 1:nrow(z_clipped_flip),
  z = t(z_clipped_flip),  # transpose for correct orientation
  col = c("white", "gray", "black"),  # NA = white, 0 = gray, 1 = black
  breaks = c(-0.1, 0.1, 0.9, 1.1),    # assign color by value (NA, 0, 1)
  axes = FALSE,
  xlab = "Year",
  ylab = "Individual",
  main = "Latent State Matrix (z)"
)

# Add axis labels
axis(1, at = 1:ncol(z_clipped_flip), labels = colnames(z_clipped_flip), las = 2, cex.axis = 0.7)
axis(2, at = 1:nrow(z_clipped_flip), labels = rownames(z_clipped_flip), las = 1, cex.axis = 0.4)

rm(z_clipped_flip, y_clipped_flip)

### --------------- FIRST SEEN VECTOR ---------------- ###

# finally, we need a vector of when each individual was first seen
first_seen <- apply(y, 1, function(row) {
  first <- which(row == 1)[1]
  if (is.na(first)) return(n_years) else return(first)
})
first_seen <- as.integer(first_seen)

first_seen_clipped <- apply(y_clipped, 1, function(row) {
  first <- which(row == 1)[1]
  if (is.na(first)) return(ncol(y_clipped)) else return(first)
})
first_seen_clipped <- as.integer(first_seen_clipped)

stop('All required matrices constructed (line 338).')

#########################################################################
### ------------------------ MODEL ONE ------------------------------ ###
#########################################################################

### Model one is the simplest version of a CJS model, where elk survival
### is assumed to be constant across age-classes and throughout time.

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

if (use_clip == TRUE){
  
  # data
  data <- list(
    y = y_clipped
  )
  
  # constants
  constants <- list(
    N = nrow(y_clipped),
    n_years = ncol(y_clipped),
    first_seen = first_seen_clipped
  )
  
  # create z init matrix
  z_init <- z_clipped
  z_init[is.na(z_init)] <- 1  # assume alive unless evidence says otherwise
  
} else {
  
  # data
  data <- list(
    y = y
  )
  
  # constants
  constants <- list(
    N = nrow(y),
    n_years = ncol(y),
    first_seen = first_seen
  )
  
  # create z init matrix
  z_init <- z
  z_init[is.na(z_init)] <- 1  # assume alive unless evidence says otherwise
}

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

stop('(Line 446): End of model one. Continue beyond line 446 to model two.')

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

if (use_clip == TRUE) {
  
  # Data
  data <- list(
    y = y_clipped
  )
  
  # Constants
  constants <- list(
    N = nrow(y_clipped),
    n_years = ncol(y_clipped),
    first_seen = first_seen_clipped
  )
  
  # Initial z matrix (assume not alive)
  z_init <- z_clipped
  z_init[is.na(z_init)] <- 0
  
  # Initial values
  inits <- list(
    phi = runif(ncol(y_clipped) - 1, 0.8, 0.99),
    p = runif(1, 0.6, 0.95),
    z = z_init
  )
  
} else {
  
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
}

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
  mutate(Year = 2001:2024)  

# plot temporal trend
ggplot(elk_temporal_phi) + 
  geom_ribbon(aes(x = Year, ymin = `2.5%`, ymax = `97.5%`), fill = "grey90") +
  geom_path(aes(x = Year, y = mean), color = "red", linewidth = 1) +
  scale_y_continuous(name = "Estimated Survival (φ)") +
  scale_x_continuous(name = "Year", limits = c(2001, 2024)) +
  labs(x = "Year", title = "Time-varying survival estimates") +
  theme_bw()

stop('(Line 580): End of model two. Continue beyond line 580 to model three.')

#########################################################################
### ----------------------- MODEL THREE ----------------------------- ###
#########################################################################
# incorporating stage-specific survival, but not allowing for temporally varying rates

### ------------------------ BUILD STAGE MATRIX ------------------------ ###

# Copy z and mask non-alive values
age_class <- z
age_class[age_class != 1] <- NA

# Assign class 1 (0–13 y) and class 2 (14+ y)
for (i in 1:nrow(age_class)) {
  alive_years <- which(age_class[i, ] == 1)
  if (length(alive_years) == 0) next
  split_idx <- min(length(alive_years), 14)  # 13 years of class 1, then switch
  age_class[i, alive_years[1:split_idx]] <- 1
  if (length(alive_years) > split_idx) {
    age_class[i, alive_years[(split_idx + 1):length(alive_years)]] <- 2
  }
}

# Get last seen year index
years_vector <- as.numeric(colnames(y))
last_seen_index <- match(year(as.Date(df$Last.Date.Alive)), years_vector)

# Create clipped version of the stage matrix

# copy stage matrix
age_class_clipped <- age_class

for (i in 1:n_indiv) {
  cap_year <- df_clean$CaptureYear[i]
  cap_idx <- which(years == cap_year)
  
  if (length(cap_idx) == 1) {
    age_class_clipped[i, 1:(cap_idx - 1)] <- NA  # Mask all years before capture
  }
}

# clip to years 2000:2024
age_class_clipped <- age_class_clipped[,19:ncol(age_class_clipped)]

# create dummy matrices that are 0s and 1s representing each class (helps with if-then
# logic in the model block, since NIMBLE doesn't support actual if-then statements)
is_class1_clipped <- array(0, dim = c(nrow(age_class_clipped), ncol(age_class_clipped)))
is_class2_clipped <- array(0, dim = c(nrow(age_class_clipped), ncol(age_class_clipped)))

is_class1_clipped[age_class_clipped == 1] <- 1
is_class2_clipped[age_class_clipped == 2] <- 1

### ------------------------ VISUALIZE STAGE MATRIX ------------------------ ###

# reset the visualizing pane (par settings from MCMC plots still active)
dev.off()

### non-clipped:
# Recode age_class for plotting
plot_matrix <- age_class
plot_matrix[plot_matrix == 1] <- 1  # class 1 (0–13 y)
plot_matrix[plot_matrix == 2] <- 2  # class 2 (14+ y)
plot_matrix[is.na(plot_matrix)] <- 0  # non-alive = white

# Flip
plot_matrix <- plot_matrix[nrow(plot_matrix):1, ]

# Define colors: white = NA, blue = class 1, orange = class 2
plot_colors <- c("white", "skyblue", "orange")

# Plot
image(
  x = 1:ncol(plot_matrix),
  y = 1:nrow(plot_matrix),
  z = t(plot_matrix),
  col = plot_colors,
  axes = FALSE,
  xlab = "Year",
  ylab = "Individual",
  main = "Age Class Matrix"
)

# Axis labels
axis(1, at = 1:ncol(plot_matrix), labels = colnames(plot_matrix), las = 2, cex.axis = 0.7)
axis(2, at = 1:nrow(plot_matrix), labels = rev(rownames(plot_matrix)), las = 1, cex.axis = 0.4)

# Add legend
legend("topright", legend = c("Not Alive", "Age 0–13", "Age 14+"),
       fill = plot_colors, cex = 0.8, border = NA)

### clipped:
# Recode age_class_clipped for plotting
plot_matrix <- age_class_clipped
plot_matrix[plot_matrix == 1] <- 1  # class 1 (0–13 y)
plot_matrix[plot_matrix == 2] <- 2  # class 2 (14+ y)
plot_matrix[is.na(plot_matrix)] <- 0  # non-alive = white

# Flip
plot_matrix <- plot_matrix[nrow(plot_matrix):1, ]

# Define colors: white = NA, blue = class 1, orange = class 2
plot_colors <- c("white", "skyblue", "orange")

# Plot
image(
  x = 1:ncol(plot_matrix),
  y = 1:nrow(plot_matrix),
  z = t(plot_matrix),
  col = plot_colors,
  axes = FALSE,
  xlab = "Year",
  ylab = "Individual",
  main = "Age Class Matrix"
)

# Axis labels
axis(1, at = 1:ncol(plot_matrix), labels = colnames(plot_matrix), las = 2, cex.axis = 0.7)
axis(2, at = 1:nrow(plot_matrix), labels = rev(rownames(plot_matrix)), las = 1, cex.axis = 0.4)

# Add legend
legend("topright", legend = c("Not Alive", "Age 0–13", "Age 14+"),
       fill = plot_colors, cex = 0.8, border = NA)

### ------------------------ NIMBLE CODE ------------------------ ###

elk_survival_stage_specific <- nimbleCode({
  # Priors
  p ~ dunif(0, 1)
  phi_1 ~ dunif(0, 1)
  phi_2 ~ dunif(0, 1)
  
  for (i in 1:N) {
    z[i, first_seen[i]] <- 1  # known alive at first detection
    
    for (t in (first_seen[i] + 1):n_years) {
      z[i, t] ~ dbern(
        z[i, t - 1] * (
          is_class1[i, t - 1] * phi_1 +
            is_class2[i, t - 1] * phi_2
        ) + 1e-10 * (1 - z[i, t - 1])
      )
    }
    
    for (t in first_seen[i]:n_years) {
      y[i, t] ~ dbern(p * z[i, t])
    }
  }
})

### ------------------------ MODEL SPECS ------------------------ ###

if (use_clip == TRUE) {
  
  # Data + constants
  data <- list(y = y_clipped)
  
  constants <- list(
    N = nrow(y_clipped),
    n_years = ncol(y_clipped),
    first_seen = first_seen_clipped,
    is_class1 = is_class1_clipped,
    is_class2 = is_class2_clipped
  )
  
  # Initial values
  z_init <- matrix(NA, nrow = nrow(y_clipped), ncol = ncol(y_clipped))
  for (i in 1:nrow(y_clipped)) {
    if (!is.na(first_seen_clipped[i]) && first_seen_clipped[i] < ncol(y_clipped)) {
      z_init[i, (first_seen_clipped[i] + 1):ncol(y_clipped)] <- 1
    }
  }
  
} else {
  
  # Data + constants
  data <- list(y = y)
  
  constants <- list(
    N = nrow(y),
    n_years = ncol(y),
    first_seen = first_seen,
    is_class1 = is_class1,
    is_class2 = is_class2
  )
  
  # Initial values
  z_init <- matrix(NA, nrow = nrow(y), ncol = ncol(y))
  for (i in 1:nrow(y)) {
    if (!is.na(first_seen[i]) && first_seen[i] < ncol(y)) {
      z_init[i, (first_seen[i] + 1):ncol(y)] <- 1
    }
  }
  
}

inits <- list(
  phi_1 = runif(1, 0.6, 0.95),
  phi_2 = runif(1, 0.6, 0.95),
  p = runif(1, 0.6, 0.95),
  z = z_init
)

params <- c("phi_1", "phi_2", "p")

### ------------------------ RUN MODEL ------------------------ ###

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

elk_model_stage_specific$summary$all.chains

### -------------------- CHECK CONVERGENCE ---------------------- ###

# Traceplots
MCMCtrace(elk_model_stage_specific$samples,
          pdf = FALSE, ind = TRUE,
          Rhat = TRUE, n.eff = TRUE,
          params = 'all')

# Summary
MCMCsummary(elk_model_stage_specific$samples, params = 'all')

stop('(Line 813): End of model three. Continue beyond line 813 to model four.')

#########################################################################
### ----------------------- MODEL FOUR ------------------------------ ###
#########################################################################
### ideally, this is the final survival model: temporally varying, stage-specific
### estimates of survival.

### ------------------------ NIMBLE CODE ------------------------ ###

elk_survival_final <- nimbleCode({
  
  for (t in 1:n_years){
    # Priors
    p[t] ~ dunif(0, 1)
    phi_1[t] ~ dunif(0, 1)
    phi_2[t] ~ dunif(0, 1)
  }

  
  for (i in 1:N) {
    z[i, first_seen[i]] <- 1  # known alive at first detection
    
    for (t in (first_seen[i] + 1):n_years) {
      z[i, t] ~ dbern(
        z[i, t - 1] * (
          is_class1[i, t - 1] * phi_1[t] +
            is_class2[i, t - 1] * phi_2[t]
        ) + 1e-10 * (1 - z[i, t - 1])
      )
    }
    
    for (t in first_seen[i]:n_years) {
      y[i, t] ~ dbern(p[t] * z[i, t])
    }
  }
})

### ------------------------ MODEL SPECS ------------------------ ###

if (use_clip == TRUE) {
  
  # Data + constants
  data <- list(y = y_clipped)
  
  constants <- list(
    N = nrow(y_clipped),
    n_years = ncol(y_clipped),
    first_seen = first_seen_clipped,
    is_class1 = is_class1_clipped,
    is_class2 = is_class2_clipped
  )
  
  # Initial values
  z_init <- matrix(NA, nrow = nrow(y_clipped), ncol = ncol(y_clipped))
  for (i in 1:nrow(y_clipped)) {
    if (!is.na(first_seen_clipped[i]) && first_seen_clipped[i] < ncol(y_clipped)) {
      z_init[i, (first_seen_clipped[i] + 1):ncol(y_clipped)] <- 1
    }
  }
  
  inits <- list(
    phi_1 = runif(ncol(y_clipped), 0.6, 0.95),
    phi_2 = runif(ncol(y_clipped), 0.6, 0.95),
    p = runif(ncol(y_clipped), 0.6, 0.95),
    z = z_init
  )
  
} else {
  
  # Data + constants
  data <- list(y = y)
  
  constants <- list(
    N = nrow(y),
    n_years = ncol(y),
    first_seen = first_seen,
    is_class1 = is_class1,
    is_class2 = is_class2
  )
  
  # Initial values
  z_init <- matrix(NA, nrow = nrow(y), ncol = ncol(y))
  for (i in 1:nrow(y)) {
    if (!is.na(first_seen[i]) && first_seen[i] < ncol(y)) {
      z_init[i, (first_seen[i] + 1):ncol(y)] <- 1
    }
  }
  
  inits <- list(
    phi_1 = runif(ncol(y), 0.6, 0.95),
    phi_2 = runif(ncol(y), 0.6, 0.95),
    p = runif(ncol(y), 0.6, 0.95),
    z = z_init
  )
  
}



params <- c("phi_1", "phi_2", "p")

### ------------------------ RUN MODEL ------------------------ ###

elk_model_final <- nimbleMCMC(
  code = elk_survival_final,
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

elk_model_final$summary$all.chains

### -------------------- CHECK CONVERGENCE ---------------------- ###

# Traceplots
MCMCtrace(elk_model_final$samples,
          pdf = FALSE, ind = TRUE,
          Rhat = TRUE, n.eff = TRUE,
          params = 'all')

# Summary
MCMCsummary(elk_model_final$samples, params = 'all')

# Extract phi_1
phi1_df <- MCMCsummary(elk_model_final$samples, params = "phi_1") %>%
  as.data.frame() %>%
  rownames_to_column("param") %>%
  mutate(
    time_index = as.numeric(gsub("phi_1\\[|\\]", "", param)),
    Year = 2000 + time_index,  # Adjust if needed
    Class = "Young (0-13 years)"
  )

# Extract phi_2
phi2_df <- MCMCsummary(elk_model_final$samples, params = "phi_2") %>%
  as.data.frame() %>%
  rownames_to_column("param") %>%
  mutate(
    time_index = as.numeric(gsub("phi_2\\[|\\]", "", param)),
    Year = 2000 + time_index,  # Adjust if needed
    Class = "Old (14+ years)"
  )

# Combine
phi_combined <- bind_rows(phi1_df, phi2_df)

# Plot
ggplot(phi_combined, aes(x = Year, y = mean, color = Class, fill = Class)) +
  geom_ribbon(aes(ymin = `2.5%`, ymax = `97.5%`), alpha = 0.3, color = NA) +
  geom_line(linewidth = 1.2) +
  scale_color_manual(
    values = c("Young (0-13 years)" = "#1f77b4",  # blue
               "Old (14+ years)"    = "#d62728")  # red
  ) +
  scale_fill_manual(
    values = c("Young (0-13 years)" = "#1f77b4",  # blue
               "Old (14+ years)"    = "#d62728")  # red
  ) +
  scale_y_continuous(name = "Estimated Survival (φ)", limits = c(0, 1)) +
  scale_x_continuous(name = "Year", limits = c(2001, 2024)) +
  labs(title = "Time-varying Survival Estimates by Age Class") +
  theme_bw() +
  theme(legend.title = element_blank())

stop('(Line 951): End of code.')

#######################################################################