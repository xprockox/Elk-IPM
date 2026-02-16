# ============================================================
# Preliminary Analysis for ESA Conference Abstract
# ============================================================

# ------------------------------------------------------------
# Packages
# ------------------------------------------------------------

library(tidyverse)
library(nimble)
library(coda)

# ------------------------------------------------------------
# Data loading
# ------------------------------------------------------------

dat <- read.csv('data/master/allSpp_Abundances.csv')
elk_n <- read.csv('data/results/elk_N_byStages_2026-02-16.csv')
elk_vrates <- read.csv('data/results/elk_vrates_2026-02-16.csv')
clim_dat <- read.csv('data/intermediate/prism_monthly_snow.csv')

# ------------------------------------------------------------
# Cougar interpolation
# Some years are missing cougar abundance.
# Fit linear model to fill internal gaps only.
# ------------------------------------------------------------

m_cougar <- lm(Cougars ~ Year, data = dat, na.action = na.exclude)

dat$Cougars_pred <- predict(m_cougar, newdata = dat)

# Replace missing cougar values with regression predictions
dat$Cougars_filled <- round(
  ifelse(is.na(dat$Cougars),
         dat$Cougars_pred,
         dat$Cougars)
)

# Keep only predicted values for visualization if needed
dat$Cougars_pred <- ifelse(is.na(dat$Cougars),
                           dat$Cougars_pred,
                           NA)

# Keep only relevant columns
dat <- dat %>%
  select(Year, Wolves, Cougars_filled, Grizzly.Bears) %>%
  rename(
    wolf_n = Wolves,
    cougar_n = Cougars_filled,
    griz_n = Grizzly.Bears
  )

# Remove incomplete rows
dat <- dat[complete.cases(dat), ]


# ------------------------------------------------------------
# Standardize predator abundances (mean = 0, SD = 1)
# Allows direct comparison of regression coefficients
# ------------------------------------------------------------

dat_scaled <- dat %>%
  mutate(
    wolf_z = scale(wolf_n)[,1],
    cougar_z = scale(cougar_n)[,1],
    griz_z = scale(griz_n)[,1]
  )

rm(m_cougar)

# ------------------------------------------------------------
# Standardize climate data
# Allows direct comparison of regression coefficients
# ------------------------------------------------------------

clim_dat_scaled <- clim_dat %>%
  mutate(
    winter_precip_z = scale(winter_precip)[,1],
    april1_winter_precip_z = scale(april1_winter_precip)[,1],
    freezing_months_z = scale(freezing_months)[,1]
  ) %>%
  rename(year = water_year)

# ------------------------------------------------------------
# Convert elk stage abundance to wide format
# ------------------------------------------------------------

elk_n_wide <- elk_n %>%
  select(year, stage, mean) %>%
  pivot_wider(names_from = stage,
              values_from = mean)


# ------------------------------------------------------------
# Clean vital rate parameter names
# IPM outputs include bracketed indices (e.g. s_c[1])
# Remove bracketed index to recover base parameter name
# ------------------------------------------------------------

elk_vrates_clean <- elk_vrates %>%
  mutate(
    param_base = str_remove(param, "\\[.*\\]"),
    param_index = as.numeric(str_extract(param, "(?<=\\[)\\d+(?=\\])"))
  )

elk_vrates_wide <- elk_vrates_clean %>%
  select(year, param_base, mean) %>%
  pivot_wider(names_from = param_base,
              values_from = mean)


# ------------------------------------------------------------
# Merge predator and elk data
# ------------------------------------------------------------

dat_all <- dat_scaled %>%
  rename(year = Year) %>%
  left_join(elk_n_wide, by = "year") %>%
  left_join(elk_vrates_wide, by = "year")

# Keep only complete rows across all species + vital rates
dat_all <- dat_all[complete.cases(dat_all), ]

head(dat_all)

# ------------------------------------------------------------
# Now merge climate data
# ------------------------------------------------------------

dat_all <- left_join(dat_all, clim_dat_scaled, by='year')
# ------------------------------------------------------------
# Logit transformations
# Necessary because survival is bounded (0,1)
# ------------------------------------------------------------

logit <- function(x) log(x / (1 - x))

dat_all <- dat_all %>%
  mutate(
    s_c_logit  = logit(s_c),   # calf survival
    s_ya_logit = logit(s_ya),  # young adult survival
    s_oa_logit = logit(s_oa)   # old adult survival
  )

# ------------------------------------------------------------
# Log transformations
# because fecundity and abundance are positive, unbounded
# ------------------------------------------------------------
dat_all <- dat_all %>%
  mutate(
    log_f_ya = log(f_ya),
    log_f_oa = log(f_oa),
    log_N_total = log(`Total Females`)
  )

# ------------------------------------------------------------
# Frequentist diagnostic of collinearity
# ------------------------------------------------------------

cor(dat_all[, c("wolf_z", "cougar_z", "griz_z")])


# ------------------------------------------------------------
# construct NIMBLE model (framework, can be adapted for any response var)
# ------------------------------------------------------------
run_model <- function(response_vector) {
  
  n_years <- length(response_vector)
  
  constants <- list(n_years = n_years)
  
  data_list <- list(
    y = response_vector,
    wolf = dat_all$wolf_z,
    cougar = dat_all$cougar_z,
    griz = dat_all$griz_z,
    winter_precip = dat_all$winter_precip_z,
    freezing = dat_all$freezing_months_z
  )
  
  code <- nimbleCode({
    
    alpha  ~ dnorm(0, sd = 2)
    
    beta_w ~ dnorm(0, sd = 1)
    beta_c ~ dnorm(0, sd = 1)
    beta_g ~ dnorm(0, sd = 1)
    
    beta_winter  ~ dnorm(0, sd = 1)
    beta_freeze  ~ dnorm(0, sd = 1)
    
    sigma ~ dunif(0, 5)
    tau <- 1 / (sigma * sigma)
    
    for (t in 1:n_years) {
      
      mu[t] <- alpha +
        beta_w * wolf[t] +
        beta_c * cougar[t] +
        beta_g * griz[t] +
        beta_winter * winter_precip[t] +
        beta_freeze * freezing[t]
      
      y[t] ~ dnorm(mu[t], tau = tau)
    }
  })
  
  inits_function <- function() {
    list(
      alpha = rnorm(1, 0, 0.5),
      beta_w = rnorm(1, 0, 0.5),
      beta_c = rnorm(1, 0, 0.5),
      beta_g = rnorm(1, 0, 0.5),
      beta_winter = rnorm(1, 0, 0.5),
      beta_freeze = rnorm(1, 0, 0.5),
      sigma = runif(1, 0.5, 1.5)
    )
  }
  
  samples <- nimbleMCMC(
    code = code,
    data = data_list,
    inits = inits_function,
    constants = constants,
    niter = 20000,
    nburnin = 5000,
    nchains = 3,
    monitors = c("alpha",
                 "beta_w", "beta_c", "beta_g",
                 "beta_winter", "beta_freeze",
                 "sigma"),
    summary = TRUE
  )
  
  return(samples$summary$all.chains)
}

# ------------------------------------------------------------
# then we can actually run the models
# ------------------------------------------------------------

# calf survival
results_s_c <- run_model(dat_all$s_c_logit)
results_s_c

# young adult survival
results_s_ya <- run_model(dat_all$s_ya_logit)
results_s_ya

# old adult survival
results_s_oa <- run_model(dat_all$s_oa_logit)
results_s_oa

# young adult fecundity
results_f_ya <- run_model(dat_all$log_f_ya)
results_f_ya

# old adult fecundity
results_f_oa <- run_model(dat_all$log_f_oa)
results_f_oa

# total abundance
results_N_total <- run_model(dat_all$log_N_total)
results_N_total
