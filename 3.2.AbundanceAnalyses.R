### Elk Abundance - State-Space Model
### Analysis script
### Last updated: Oct. 27, 2025
### Contact: xprockox@gmail.com

############################################################################################
### packages

library(nimble)
library(MCMCvis)
library(coda)
library(dplyr)
library(tidyr)
library(stringr)
library(ggplot2)

############################################################################################
### data import

dat <- read.csv('data/intermediate/abundanceEstimates_stages.csv')

############################################################################################
### nimble code

# first specify study duration
n_years <- length(dat)

# --- NIMBLE MODEL CODE ---
elk_code <- nimbleCode({
  
  ### Priors (vital rates)
  # Survival (probabilities)
  s_c ~ dbeta(1, 1) # calf survival to next year (becomes Young)
  s_y ~ dbeta(1, 1) # young adult survival
  s_o ~ dbeta(1, 1) # old adult survival
  
  # Progression from Young -> Old (conditional on survival)
  g_y ~ dbeta(1, 1) # fraction of surviving young that age into old class
  
  # Fecundity (calves produced per female/individual per year)
  # If you have sex-specific adults, use female-only or add a sex ratio.
  f_y ~ dgamma(1, 0.1) # calves per Young
  f_o ~ dgamma(1, 0.1) # calves per Old
  
  ### Observation error (lognormal, one sd per stage)
  sigma_obs_c ~ dunif(0, 2)
  sigma_obs_y ~ dunif(0, 2)
  sigma_obs_o ~ dunif(0, 2)
  tau_obs_c <- 1 / (sigma_obs_c^2)
  tau_obs_y <- 1 / (sigma_obs_y^2)
  tau_obs_o <- 1 / (sigma_obs_o^2)
  
  ### Initial abundance priors (Poisson with diffuse gamma hyperpriors)
  lambda_init_c ~ dgamma(0.001, 0.001)
  lambda_init_y ~ dgamma(0.001, 0.001)
  lambda_init_o ~ dgamma(0.001, 0.001)
  
  N_c[1] ~ dpois(lambda_init_c)
  N_y[1] ~ dpois(lambda_init_y)
  N_o[1] ~ dpois(lambda_init_o)
  
  # Observation model at t = 1
  obs_calf[1] ~ dlnorm(log(N_c[1] + 1e-6), tau_obs_c)
  obs_young[1] ~ dlnorm(log(N_y[1] + 1e-6), tau_obs_y)
  obs_old[1] ~ dlnorm(log(N_o[1] + 1e-6), tau_obs_o)
  
  ### State process + observation for t + 1, 2, 3 ... length(n_years)
  for (t in 1:(n_years-1)) {

    # Calves at t+1 are produced by young and old at t
    mu_c[t+1] <- f_y * N_y[t] + f_o * N_o[t]
    
    # Young at t+1 come from (1) surviving calves, and (2) surviving young that DON'T progress
    mu_y[t+1] <- s_c * N_c[t] + s_y * (1 - g_y) * N_y[t]
    
    # Old at t+1 come from (1) progressing surviving young, (2) surviving old
    mu_o[t+1] <- s_y * g_y * N_y[t] + s_o * N_o[t]
    
    # Process stochasticity (Poisson stage dynamics)
    N_c[t+1] ~ dpois(mu_c[t+1] + 1e-9)
    N_y[t+1] ~ dpois(mu_y[t+1] + 1e-9)
    N_o[t+1] ~ dpois(mu_o[t+1] + 1e-9)
    
    # Observation model (lognormal around true abundance)
    obs_calf[t+1] ~ dlnorm(log(N_c[t+1] + 1e-6), tau_obs_c)
    obs_young[t+1] ~ dlnorm(log(N_y[t+1] + 1e-6), tau_obs_y)
    obs_old[t+1] ~ dlnorm(log(N_o[t+1] + 1e-6), tau_obs_o)
  }
  
  # Derived total population (females only)
  for (t in 1:n_years) {
    N_tot[t] <- N_c[t] + N_y[t] + N_o[t]
  }
})

############################################################################################
### data and constants and params, etc.

# constants
elk_constants <- list(n_years = n_years)

# data
elk_data <- list(
  obs_calf  = dat$n_calf,
  obs_young = dat$n_cow_youngadult,
  obs_old   = dat$n_cow_oldadult
)

# inits
make_inits <- function() {
  init_Nc <- ifelse(is.na(dat$n_calf),  
                    pmax(1, round(mean(dat$n_calf, na.rm = TRUE))),  
                    pmax(1, round(dat$n_calf)))
  init_Ny <- ifelse(is.na(dat$n_cow_youngadult), 
                    pmax(1, round(mean(dat$n_cow_youngadult, na.rm = TRUE))), 
                    pmax(1, round(dat$n_cow_youngadult)))
  init_No <- ifelse(is.na(dat$n_cow_oldadult),   
                    pmax(1, round(mean(dat$n_cow_oldadult, na.rm = TRUE))),   
                    pmax(1, round(dat$n_cow_oldadult)))
  
  list(
    s_c = 0.6, s_y = 0.9, s_o = 0.85,
    g_y = 0.15,
    f_y = 0.2, f_o = 0.3,
    sigma_obs_c = 0.3, sigma_obs_y = 0.3, sigma_obs_o = 0.3,
    lambda_init_c = max(1, round(init_Nc[1])),
    lambda_init_y = max(1, round(init_Ny[1])),
    lambda_init_o = max(1, round(init_No[1])),
    N_c = pmax(1, init_Nc),
    N_y = pmax(1, init_Ny),
    N_o = pmax(1, init_No)
  )
}

# params
params <- c("s_c","s_y","s_o","g_y","f_y","f_o",
            "sigma_obs_c","sigma_obs_y","sigma_obs_o",
            "N_c","N_y","N_o","N_tot")

# set seed
set.seed(17)

# set model details
nc = 3
ni = 1000000
nb = 200000

# run MCMC
elk_mod1 <- nimbleMCMC(
  code = elk_code,
  data = elk_data,
  constants = elk_constants,
  inits = make_inits,
  monitors = params,
  nchains = nc,
  niter = ni,
  nburnin = nb,
  summary=TRUE
)

# summary
round(elk_mod1$summary$all.chains,2)

############################################################################################
### for some reason, the following lines:

# MCMCtrace(elk_mod1,
#           pdf   = FALSE,
#           ind   = TRUE,
#           Rhat  = TRUE,
#           n.eff = TRUE,
#           params = 'all')
# 
# 
# MCMCsummary(elk_mod1,
#             params = 'all')

### produce the following error:

# Error in coda::mcmc.list(lapply(object, function(x) coda::mcmc(x))) : 
# Different start, end or thin values in each chain

# so let's pull the chains out one-by-one, trim them to be equal, and then stick them back together
############################################################################################

# -----------------------------
# A) Normalize chains' shapes
# -----------------------------

# 1) Each element of elk_mod1$samples is a chain; convert each of these to a plain matrix
mats <- lapply(elk_mod1$samples, function(ch) as.matrix(ch))

# 2) Trim all chains to the same number of iterations (using the minimum length of chain)
lens <- sapply(mats, nrow)
L <- min(lens)
mats_trim <- lapply(mats, function(M) tail(M, L))

# 3) Rebuild a consistent mcmc.list
mlist <- mcmc.list(lapply(mats_trim, function(M) mcmc(M, start = 1, end = L, thin = 1)))

# 4) create an alias (copy of the data) for clarity in the next step
ml <- mlist  

# -----------------------------------------
# B) Drop parameters that contain any NA/NaN
# -----------------------------------------

# Convert each chain back to a matrix for NA screening
mats <- lapply(ml, as.matrix)

# Identify columns (parameters) that have no NA/NaN in ANY chain.
# We compute, per chain, the set of columns with zero NA/NaN, then intersect across chains
keep_cols <- Reduce(intersect, lapply(mats, function(M) {
  colnames(M)[colSums(is.na(M) | is.nan(M)) == 0]
}))

# Build a "clean" mcmc.list that retains only the safe parameters (same columns in every chain)
# Also set consistent start/end/thin to keep MCMCvis happy.
ml_clean <- mcmc.list(lapply(mats, function(M) {
  mcmc(M[, keep_cols, drop = FALSE], start = 1, end = nrow(M), thin = 1)
}))

# -----------------------------------------
# C) Summarize and visualize with MCMCvis
# -----------------------------------------
round(MCMCsummary(ml_clean, params = 'all'),2)
MCMCtrace(ml_clean, pdf = FALSE, ind = TRUE, Rhat = TRUE, n.eff = TRUE, params = 'all')

############################################################################################

# Get summaries for all N parameters
N_summ <- MCMCsummary(ml_clean, params=c('N_c', 'N_y', 'N_o', 'N_tot'))   # extracts mean, sd, 2.5%, 50%, 97.5%

# Convert rownames to columns
N_summ <- N_summ %>%
  tibble::rownames_to_column("param") %>%
  # Extract stage (N_c, N_y, N_o, N_tot) and time index number
  mutate(
    stage = str_extract(param, "^N_[a-zA-Z]+"),
    t     = as.integer(str_extract(param, "(?<=\\[)\\d+(?=\\])"))
  ) %>%
  select(stage, t, mean = mean, low = `2.5%`, high = `97.5%`) %>%
  arrange(stage, t)

years <- 1995:2023   # or however many years you have
N_summ <- N_summ %>% mutate(year = years[t])

N_summ <- N_summ %>%
  mutate(stage = recode(stage,
                        "N_c"   = "Calf",
                        "N_y"   = "Young Adult",
                        "N_o"   = "Old Adult",
                        "N_tot" = "Total"))


dat_long <- dat %>%
  pivot_longer(
    cols = -c(X,year),          # everything except year becomes stage/value pairs
    names_to = "stage",
    values_to = "value"
  ) %>%
  mutate(stage = recode(stage,
                        n_calf  = "Calf",
                        n_cow_youngadult = "Young Adult",
                        n_cow_oldadult   = "Old Adult",
                        n_total = "Total"))



dat_long$stage <- factor(dat_long$stage, levels = unique(N_summ$stage))

dat_long <- dat_long[dat_long$stage %in% c('Calf', 'Old Adult', 'Total', 'Young Adult'),]

ggplot(N_summ, aes(x = year, y = mean, group = stage)) +
  geom_ribbon(aes(ymin = low, ymax = high), alpha = 0.2) +
  geom_line(size = 1) +
  geom_point(data = dat_long, aes(y = value), color = "red", size = 2) +
  geom_line(data = dat_long, aes(y = value), color = "red", linetype = 2) +
  facet_wrap(~ stage, scales = "free_y") +
  theme_bw() +
  labs(x = "Year", y = "Abundance",
       title = "Posterior Population Estimates with Validation Data",
       subtitle = "Ribbon = 95% credible interval, Line = posterior mean, Red = observed")

############################################################################################
### what about demographic rates?

dem_rates <- MCMCsummary(ml_clean, 
                         params = c("s_c", "s_y", "s_o", "g_y", "f_y", "f_o")) %>%
  as.data.frame() %>%
  tibble::rownames_to_column("param")

dem_rates_plot <- dem_rates %>%
  select(param, mean, `2.5%`, `97.5%`) %>%
  rename(low = `2.5%`, high = `97.5%`)

ggplot(dem_rates_plot, aes(x = param, y = mean)) +
  geom_pointrange(aes(ymin = low, ymax = high)) +
  coord_flip() +
  theme_bw() +
  labs(x = "Vital Rate", y = "Posterior Mean (95% CI)",
       title = "Estimated Vital Rates")
