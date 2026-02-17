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

# Create indicator for interpolated years
dat <- dat %>%
  mutate(
    interpolated = ifelse(is.na(Cougars), "Interpolated", "Observed")
  )

# Plot
ggplot(dat, aes(x = Year, y = Cougars_filled)) +
  geom_line(color = "grey40", linewidth = 0.8) +
  geom_point(aes(color = interpolated), size = 1.5) +
  scale_color_manual(
    values = c("Observed" = "black",
               "Interpolated" = "firebrick")
  ) +
  theme_bw(base_size = 12) +
  labs(
    y = "Cougar abundance",
    color = NULL,
    title = "Cougar abundance over time",
    subtitle = "Red points indicate interpolated years"
  )

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
# Pulling out residuals of grizzlies and cougars
# ------------------------------------------------------------

# Since grizzlies and cougars have some level of collinearity
# with wolves, we can instead use their residuals to find the effect
# of each of these species independent of the wolf signal.

# For grizzlies, we find the residual of griz ~ wolf
griz_resid <- resid(lm(griz_z ~ wolf_z, data = dat_all))
griz_resid <- scale(griz_resid)[,1]

# Then to separate cougars from both wolves + grizzlies, 
# we find the residual of cougar ~ wolf + griz
cougar_resid <- resid(lm(cougar_z ~ wolf_z + griz_z, data = dat_all))
cougar_resid <- scale(cougar_resid)[,1]

# ------------------------------------------------------------
# construct NIMBLE model (framework, can be adapted for any response var)
# ------------------------------------------------------------
run_model <- function(response_vector) {
  
  n_years <- length(response_vector)
  
  constants <- list(n_years = n_years)
  
  data_list <- list(
    y = response_vector,
    wolf = dat_all$wolf_z,
    cougar = cougar_resid,
    griz = griz_resid,
    winter_precip = dat_all$winter_precip_z,
    freezing = dat_all$freezing_months_z
  )
  
  code <- nimbleCode({
    
    # ---------------------------
    # Fixed effects
    # ---------------------------
    
    alpha ~ dnorm(0, sd = 5)
    
    beta_w ~ dnorm(0, sd = 1)
    beta_c ~ dnorm(0, sd = 1)
    beta_g ~ dnorm(0, sd = 1)
    
    beta_winter ~ dnorm(0, sd = 1)
    beta_freeze ~ dnorm(0, sd = 1)
    
    # ---------------------------
    # Residual variance only
    # ---------------------------
    
    sigma_resid ~ dunif(0, 2)
    tau_resid <- 1 / (sigma_resid * sigma_resid)
    
    for (t in 1:n_years) {
      
      mu[t] <- alpha +
        beta_w * wolf[t] +
        beta_c * cougar[t] +
        beta_g * griz[t] +
        beta_winter * winter_precip[t] +
        beta_freeze * freezing[t]
      
      y[t] ~ dnorm(mu[t], tau = tau_resid)
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
      sigma_resid = runif(1, 0.1, 0.5)
    )
  }
  
  samples <- nimbleMCMC(
    code = code,
    data = data_list,
    inits = inits_function,
    constants = constants,
    niter = 30000,
    nburnin = 10000,
    nchains = 3,
    monitors = c("alpha",
                 "beta_w", "beta_c", "beta_g",
                 "beta_winter", "beta_freeze",
                 "sigma_resid"),
    summary = FALSE
  )
  
  # ---------------------------------------------------
  # Convert to coda for diagnostics
  # ---------------------------------------------------
  
  mcmc_list <- mcmc.list(lapply(samples, mcmc))
  
  posterior_summary <- summary(mcmc_list)
  
  means   <- posterior_summary$statistics[, "Mean"]
  sds     <- posterior_summary$statistics[, "SD"]
  medians <- posterior_summary$quantiles[, "50%"]
  lower   <- posterior_summary$quantiles[, "2.5%"]
  upper   <- posterior_summary$quantiles[, "97.5%"]
  
  rhat <- gelman.diag(mcmc_list, multivariate = FALSE)$psrf[, "Point est."]
  
  ess <- effectiveSize(mcmc_list)
  
  results_table <- data.frame(
    Mean = means,
    Median = medians,
    SD = sds,
    CI_2.5 = lower,
    CI_97.5 = upper,
    Rhat = rhat,
    ESS = ess
  )
  
  return(results_table)
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

# ------------------------------------------------------------
# and visualize the results
# ------------------------------------------------------------

# Add vital rate labels
results_s_c$rate <- "Calf survival"
results_s_ya$rate <- "Young adult survival"
results_s_oa$rate <- "Old adult survival"
results_f_ya$rate <- "Young adult fecundity"
results_f_oa$rate <- "Old adult fecundity"
results_N_total$rate <- "Total abundance"

# Combine
all_results <- bind_rows(
  results_s_c,
  results_s_ya,
  results_s_oa,
  results_f_ya,
  results_f_oa,
  results_N_total,
)

# The rownames currently contain parameter names
all_results$parameter <- rownames(all_results)
all_results$parameter <- sub("\\.\\.\\.[0-9]+$", "", all_results$parameter)


# Create predictor column
all_results <- all_results %>%
  mutate(
    predictor = case_when(
      parameter == "beta_w"        ~ "Wolves",
      parameter == "beta_g"        ~ "Grizzlies",
      parameter == "beta_c"        ~ "Cougars",
      parameter == "beta_winter"   ~ "Winter precipitation",
      parameter == "beta_freeze"   ~ "Freezing months",
      TRUE ~ NA_character_
    )
  )


all_results <- all_results %>%
  filter(!is.na(predictor))

# Then create the coefficient plot
ggplot(all_results,
       aes(x = rate,
           y = Mean,
           ymin = CI_2.5,
           ymax = CI_97.5)) +
  geom_hline(yintercept = 0,
             linetype = "dashed",
             color = "grey40") +
  geom_pointrange(size = 0.6) +
  facet_wrap(~ predictor, scales = "free_y") +
  coord_flip() +
  theme_bw(base_size = 11) +
  labs(
    x = NULL
  )+
  theme(
    strip.text = element_blank(),
    strip.background = element_blank()
  )

ggplot(dat_all)+
  geom_line(aes(x=year, y=wolf_z))+
  geom_line(aes(x=year, y=griz_z))


# ------------------------------------------------------------
# Plot time series of all standardized predictors
# ------------------------------------------------------------

predictor_ts <- dat_all %>%
  select(year,
         wolf_z,
         griz_z,
         cougar_z,
         winter_precip_z,
         freezing_months_z,
         `Total Females`) %>%
  rename(`Total Female Elk` = `Total Females`) %>%
  pivot_longer(-year,
               names_to = "predictor",
               values_to = "value") 

ggplot(predictor_ts,
       aes(x = year,
           y = value)) +
  geom_line(size = 1, color = "black") +
  facet_wrap(~ predictor, scales = "free_y") +
  theme_bw(base_size = 12) +
  labs(
    x = "Year",
    y = ""
  )


ggplot(dat_all, aes(x=wolf_n, y=`Total Females`))+
  geom_point()+
  geom_smooth(method = "lm", se = TRUE) +
  theme_bw(base_size = 12) +
  labs(
    x = "NR Wolf Abundance",
    y = "NR Elk Abundance (female-only)",
    title = "NR Wolf and Elk Abundance Data from 2000-2023"
  )
