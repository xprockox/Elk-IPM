### CJS Model for Elk Adult Survival
### Last updated: Oct. 27, 2025
### Contact: xprockox@gmail.com

#########################################################################
### -------------------------- SET-UP ------------------------------- ###
#########################################################################

### --------------- PACKAGES  ---------------- ###
library(dplyr)
library(lubridate)
library(tidyr)
library(tidyverse)
library(MCMCvis)
library(ggplot2)
library(nimble)

### --------------- DATA IMPORT ---------------- ###
load('data/intermediate/adultSurvival_cjsMatrices.rData')

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

stop('(Line 107): End of model one. Continue beyond line 107 to model two.')

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
  mutate(Year = 2001:2025)  

# plot temporal trend
ggplot(elk_temporal_phi) + 
  geom_ribbon(aes(x = Year, ymin = `2.5%`, ymax = `97.5%`), fill = "grey90") +
  geom_path(aes(x = Year, y = mean), color = "red", linewidth = 1) +
  scale_y_continuous(name = "Estimated Survival (φ)") +
  scale_x_continuous(name = "Year", limits = c(2001, 2024)) +
  labs(x = "Year", title = "Time-varying survival estimates") +
  theme_bw()

stop('(Line 213): End of model two. Continue beyond line 213 to save results, or
     skip to line 228 to begin model three.')

elk_temporal_phi <- elk_temporal_phi %>%
  select(Year, mean, `2.5%`, `97.5%`)
  
write.csv(elk_temporal_phi, 'data/intermediate/temp_phi_noStages.csv')

#########################################################################
### ----------------------- MODEL THREE ----------------------------- ###
#########################################################################
# incorporating stage-specific survival, but not allowing for temporally varying rates

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

# Data + constants
data <- list(
  y = y_clipped
)
  
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

stop('(Line 312): End of model three. Continue beyond line 312 to model four.')

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
          is_class1[i, t - 1] * phi_1[t-1] +
            is_class2[i, t - 1] * phi_2[t-1]
        ) + 1e-10 * (1 - z[i, t - 1])
      )
    }
    
    for (t in (first_seen[i] + 1):n_years) {
      y[i, t] ~ dbern(p[t-1] * z[i, t])
    }
  }
})

### ------------------------ MODEL SPECS ------------------------ ###

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
    Class = "Young (0-14 years)"
  )

# Extract phi_2
phi2_df <- MCMCsummary(elk_model_final$samples, params = "phi_2") %>%
  as.data.frame() %>%
  rownames_to_column("param") %>%
  mutate(
    time_index = as.numeric(gsub("phi_2\\[|\\]", "", param)),
    Year = 2000 + time_index,  # Adjust if needed
    Class = "Old (>14 years)"
  )

# Combine
phi_combined <- bind_rows(phi1_df, phi2_df)

# Plot
ggplot(phi_combined, aes(x = Year, y = mean, color = Class, fill = Class)) +
  geom_ribbon(aes(ymin = `2.5%`, ymax = `97.5%`), alpha = 0.3, color = NA) +
  geom_line(linewidth = 1.2) +
  scale_color_manual(
    values = c("Young (0-14 years)" = "#1f77b4",  # blue
               "Old (>14 years)"    = "#d62728")  # red
  ) +
  scale_fill_manual(
    values = c("Young (0-14 years)" = "#1f77b4",  # blue
               "Old (>14 years)"    = "#d62728")  # red
  ) +
  scale_y_continuous(name = "Estimated Survival (φ)", limits = c(0, 1)) +
  scale_x_continuous(name = "Year", limits = c(2001, 2024)) +
  labs(title = "Time-varying Survival Estimates by Age Class") +
  theme_bw() +
  theme(legend.title = element_blank())

stop('(Line 449): End of code. Continue beyond line 449 to save final results.')

adult_survival <- MCMCsummary(elk_model_final$samples, params = 'all')
write.csv(adult_survival, 'data/intermediate/temp_phi_stages.csv')

#######################################################################