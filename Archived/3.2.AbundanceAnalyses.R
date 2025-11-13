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
library(tibble)

############################################################################################
### data import

dat <- read.csv('data/intermediate/abundanceEstimates_stages.csv')

############################################################################################
### nimble code

# ---- study duration (use a vector, not the data.frame itself) ----
n_years <- length(dat$n_calf)   # was: length(dat)

elk_code <- nimbleCode({
  
  ## -----------------------------
  ## Priors for YEAR-VARYING rates (non-centered)
  ## -----------------------------
  # Intercepts on link scales: weakly-informative around your beliefs
  alpha_sc ~ dnorm(qlogis(0.30), 1/0.5^2)   # calf survival (logit)
  alpha_sy ~ dnorm(qlogis(0.90), 1/0.5^2)   # young survival (logit)
  alpha_so ~ dnorm(qlogis(0.80), 1/0.5^2)   # old survival (logit)
  alpha_gy ~ dnorm(qlogis(0.15), 1/0.5^2)   # young->old (logit)
  
  alpha_fy ~ dnorm(log(0.85),  1/0.5^2)     # fecundity from Young (log)
  alpha_fo ~ dnorm(log(0.5),  1/0.5^2)     # fecundity from Old (log)
  
  # SD priors: half-normal via truncated normal (regularize year-to-year wiggle)
  sigma_sc ~ T(dnorm(0, 1/0.4^2), 0, )
  sigma_sy ~ T(dnorm(0, 1/0.3^2), 0, )
  sigma_so ~ T(dnorm(0, 1/0.3^2), 0, )
  sigma_gy ~ T(dnorm(0, 1/0.3^2), 0, )
  sigma_fy ~ T(dnorm(0, 1/0.5^2), 0, )
  sigma_fo ~ T(dnorm(0, 1/0.5^2), 0, )
  
  # Non-centered year effects: eps_*_std ~ N(0,1); multiply by sigma_* on link scale
  for (t in 1:(n_years-1)) {
    eps_sc_std[t] ~ dnorm(0,1)
    eps_sy_std[t] ~ dnorm(0,1)
    eps_so_std[t] ~ dnorm(0,1)
    eps_gy_std[t] ~ dnorm(0,1)
    eps_fy_std[t] ~ dnorm(0,1)
    eps_fo_std[t] ~ dnorm(0,1)
    
    # transformed yearly rates
    logit(s_c[t]) <- alpha_sc + sigma_sc * eps_sc_std[t]
    logit(s_y[t]) <- alpha_sy + sigma_sy * eps_sy_std[t]
    logit(s_o[t]) <- alpha_so + sigma_so * eps_so_std[t]
    logit(g_y[t]) <- alpha_gy + sigma_gy * eps_gy_std[t]
    
    log(f_y[t])   <- alpha_fy + sigma_fy * eps_fy_std[t]
    log(f_o[t])   <- alpha_fo + sigma_fo * eps_fo_std[t]
  }
  
  ## -----------------------------
  ## Observation error (keep away from 0)
  ## -----------------------------
  sigma_obs_c ~ dunif(0.05, 2)
  sigma_obs_y ~ dunif(0.05, 2)
  sigma_obs_o ~ dunif(0.05, 2)
  tau_obs_c <- 1 / (sigma_obs_c^2)
  tau_obs_y <- 1 / (sigma_obs_y^2)
  tau_obs_o <- 1 / (sigma_obs_o^2)
  
  ## -----------------------------
  ## Initial abundances (t = 1)
  ## -----------------------------
  lambda_init_c ~ dgamma(0.001, 0.001)
  lambda_init_y ~ dgamma(0.001, 0.001)
  lambda_init_o ~ dgamma(0.001, 0.001)
  
  N_c[1] ~ dpois(lambda_init_c)
  N_y[1] ~ dpois(lambda_init_y)
  N_o[1] ~ dpois(lambda_init_o)
  
  # Observation at t = 1
  obs_calf[1]  ~ dlnorm(log(N_c[1] + 1e-6), tau_obs_c)
  obs_young[1] ~ dlnorm(log(N_y[1] + 1e-6), tau_obs_y)
  obs_old[1]   ~ dlnorm(log(N_o[1] + 1e-6), tau_obs_o)
  
  ## -----------------------------
  ## State process + observations
  ## -----------------------------
  for (t in 1:(n_years-1)) {
    
    # Expected next-year abundances
    mu_c[t+1] <- f_y[t] * N_y[t] + f_o[t] * N_o[t]
    mu_y[t+1] <- s_c[t] * N_c[t] + s_y[t] * (1 - g_y[t]) * N_y[t]
    mu_o[t+1] <- s_y[t] * g_y[t] * N_y[t] + s_o[t] * N_o[t]
    
    # Process stochasticity (tiny positive floors)
    N_c[t+1] ~ dpois(max(1e-9, mu_c[t+1]))
    N_y[t+1] ~ dpois(max(1e-9, mu_y[t+1]))
    N_o[t+1] ~ dpois(max(1e-9, mu_o[t+1]))
    
    # Observation model
    obs_calf[t+1]  ~ dlnorm(log(N_c[t+1] + 1e-6), tau_obs_c)
    obs_young[t+1] ~ dlnorm(log(N_y[t+1] + 1e-6), tau_obs_y)
    obs_old[t+1]   ~ dlnorm(log(N_o[t+1] + 1e-6), tau_obs_o)
  }
  
  # Derived totals
  for (t in 1:n_years) {
    N_tot[t] <- N_c[t] + N_y[t] + N_o[t]
  }
})

############################################################################################

### RUN MODEL

## -----------------------------
## constants
## -----------------------------
elk_constants <- list(n_years = n_years)

## -----------------------------
## data
## -----------------------------
elk_data <- list(
  obs_calf  = dat$n_calf,
  obs_young = dat$n_cow_youngadult,
  obs_old   = dat$n_cow_oldadult
)

## -----------------------------
## inits  (match non-centered paramization)
## -----------------------------
ilogit <- function(x) 1/(1+exp(-x))

make_inits <- function() {
  # seed latent states near observations
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
    # intercepts centered on your beliefs
    alpha_sc = qlogis(0.60),
    alpha_sy = qlogis(0.90),
    alpha_so = qlogis(0.85),
    alpha_gy = qlogis(0.15),
    alpha_fy = log(0.20),
    alpha_fo = log(0.05),
    
    # SDs (start reasonably)
    sigma_sc = 0.20, sigma_sy = 0.15, sigma_so = 0.15,
    sigma_gy = 0.15, sigma_fy = 0.30, sigma_fo = 0.30,
    
    # standard-normal year effects start at 0
    eps_sc_std = rep(0, n_years - 1),
    eps_sy_std = rep(0, n_years - 1),
    eps_so_std = rep(0, n_years - 1),
    eps_gy_std = rep(0, n_years - 1),
    eps_fy_std = rep(0, n_years - 1),
    eps_fo_std = rep(0, n_years - 1),
    
    # observation SDs
    sigma_obs_c = 0.30, sigma_obs_y = 0.30, sigma_obs_o = 0.30,
    
    # initial abundances
    lambda_init_c = max(1, round(init_Nc[1])),
    lambda_init_y = max(1, round(init_Ny[1])),
    lambda_init_o = max(1, round(init_No[1])),
    N_c = pmax(1, init_Nc),
    N_y = pmax(1, init_Ny),
    N_o = pmax(1, init_No)
  )
}

## -----------------------------
## parameters to monitor
## -----------------------------
params <- c(
  # hyperparameters
  "alpha_sc","alpha_sy","alpha_so","alpha_gy","alpha_fy","alpha_fo",
  "sigma_sc","sigma_sy","sigma_so","sigma_gy","sigma_fy","sigma_fo",
  # yearly vital rates
  "s_c","s_y","s_o","g_y","f_y","f_o",
  # latent states (and total)
  "N_c","N_y","N_o","N_tot"
)



set.seed(17)
nc <- 3
ni <- 1000000
nb <- 200000
# th <- 10  # optional thinning

elk_mod1 <- nimbleMCMC(
  code      = elk_code,
  data      = elk_data,
  constants = elk_constants,
  inits     = make_inits,
  monitors  = params,
  nchains   = nc,
  niter     = ni,
  nburnin   = nb,
  # thin      = th,
  summary   = TRUE
)

## -----------------------------
## quick summary table
## -----------------------------
round(elk_mod1$summary$all.chains, 2)

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
### what about vital rates?

vrates <- MCMCsummary(ml_clean,
                      params = c("s_c","s_y","s_o","g_y","f_y","f_o")) %>%
  as.data.frame() %>%
  rownames_to_column("param") %>%
  rename(mean = mean, low = `2.5%`, high = `97.5%`)

vrates <- vrates %>%
  mutate(
    year_index = as.integer(str_extract(param, "(?<=\\[)\\d+(?=\\])")),
    rate = str_extract(param, "^[^\\[]+")  # remove bracket indices
  )

vrates$rate <- factor(vrates$rate,
                      levels = c("s_c","s_y","s_o","g_y","f_y","f_o"),
                      labels = c("Calf survival (s_c)",
                                 "Young survival (s_y)",
                                 "Old survival (s_o)",
                                 "Young→old transition (g_y)",
                                 "Fecundity (young) (f_y)",
                                 "Fecundity (old) (f_o)"))
vrates$year <- rep(1996:2023, 6)

ggplot(vrates, aes(x = year, y = mean)) +
  geom_ribbon(aes(ymin = low, ymax = high), alpha = 0.2) +
  geom_line(size = 0.9) +
  facet_wrap(~ rate, scales = "free_y") +
  theme_bw() +
  labs(x = "Year", y = "Estimated value",
       title = "Posterior Time-Varying Vital Rates (95% Credible Intervals)")

