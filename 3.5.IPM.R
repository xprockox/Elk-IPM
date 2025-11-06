### Elk IPM Main Script
### Last updated: Nov. 4, 2025
### Contact: xprockox@gmail.com

############################################################################################
################# --------------- PACKAGES AND SET-UP ---------------- #####################
############################################################################################
library(dplyr)
library(lubridate)
library(tidyr)
library(tidyverse)
library(MCMCvis)
library(ggplot2)
library(nimble)

############################################################################################
#################### --------------- DATA IMPORT ---------------- ##########################
############################################################################################
load('data/intermediate/adultSurvival_cjsMatrices.rData')
dat_n <- read.csv('data/intermediate/abundanceEstimates_stages.csv')
dat_fec <- read.csv('data/intermediate/fecundity.csv')

## 2025-10-30: Found an issue where there are latent state histories that are incorrect:
## e.g., one individual in z could be 1, 1, 1, 1, NA, NA, NA, 1, 1
## but if we know they're alive at the 8th and 9th timesteps, those NAs should also be 1s

z_fixed <- z # make a copy to work on

# fix the copy
for (i in 1:nrow(z)) {
  first_det <- which(z[i, ] == 1)[1]
  last_det <- max(which(z[i, ] == 1))
  
  if (!is.na(first_det) && !is.na(last_det) && first_det < last_det) {
    # Fill all NAs between first and last detection with 1
    z_fixed[i, first_det:last_det] <- 1
  }
}

# Use the fixed version
z <- z_fixed
rm(z_fixed)

# we have different years in the CJS and abundance data, so for now let's
# only use years that are shared across both datasets
shared_years <- intersect(as.numeric(dat_n$year),as.numeric(colnames(z)))
dat_n <- dat_n[dat_n$year %in% shared_years,]

colnames(is_class1) <- colnames(z)
colnames(is_class2) <- colnames(z)

z <- z[,which(colnames(z) %in% shared_years)]
y <- y[,which(colnames(y) %in% shared_years)]
is_class1 <- is_class1[,which(colnames(is_class1) %in% shared_years)]
is_class2 <- is_class2[,which(colnames(is_class2) %in% shared_years)]

# the same is true for dat_fec (has more years), so we need to remove years that aren't shared
dat_fec <- dat_fec[dat_fec$year %in% shared_years,]

### Nov. 5, 2025: The total number of elk (n_total) shouldn't be used because this doesn't exclude bulls
# instead, we want to approximate it as the total number of cows + half the total number of calves
dat_n$n_female <- dat_n$n_cow + (dat_n$n_calf/2)

############################################################################################
#################### --------------- NIMBLE CODE ---------------- ##########################
############################################################################################
# define years
n_years <- length(dat_n$n_calf)   

elk_ipm <- nimbleCode({
  
  ## -----------------------------
  ## (1) STATE-SPACE IPM 
  ## -----------------------------
  # priors for intercepts (mean values of vital rates)
  alpha_sc ~ dnorm(qlogis(0.22), 1/0.5^2)
  alpha_sy ~ dnorm(qlogis(0.90), 1/0.5^2)
  alpha_so ~ dnorm(qlogis(0.80), 1/0.5^2)
  alpha_gy ~ dnorm(qlogis(0.15), 1/0.5^2) # only doing survival and growth because fecundity has its own model
  
  # priors for SD of intercepts (controls variation in vital rates away from mean values)
  sigma_sc ~ T(dnorm(0, 1/0.4^2), 0, )
  sigma_sy ~ T(dnorm(0, 1/0.3^2), 0, )
  sigma_so ~ T(dnorm(0, 1/0.3^2), 0, )
  sigma_gy ~ T(dnorm(0, 1/0.3^2), 0, )
  
  # estimation of vital rates using mean value +/- some variation 
  # (sigma (SD) and eps (stochasticity))
  for (t in 1:(n_years)) {
    eps_sc_std[t] ~ dnorm(0,1)
    eps_sy_std[t] ~ dnorm(0,1)
    eps_so_std[t] ~ dnorm(0,1)
    eps_gy_std[t] ~ dnorm(0,1)
    
    logit(s_c[t]) <- alpha_sc + sigma_sc * eps_sc_std[t]
    logit(s_y[t]) <- alpha_sy + sigma_sy * eps_sy_std[t]
    logit(s_o[t]) <- alpha_so + sigma_so * eps_so_std[t]
    logit(g_y[t]) <- alpha_gy + sigma_gy * eps_gy_std[t]
  }
  
  sigma_obs_total ~ dunif(0.05, 2)
  tau_obs_total <- 1 / (sigma_obs_total^2)
  
  # initial expected values of stage-specific abundances
  lambda_init_c ~ dgamma(11.1, 0.00454)   # mean ≈ 2451 (dat_n$n_calf[1] = 2451)
  lambda_init_y ~ dgamma(11.1, 0.00118)   # mean ≈ 9418 (dat_n$n_cow_youngadult[1] = 9418)
  lambda_init_o ~ dgamma(11.1, 0.0111)    # mean ≈ 997 (dat_n$n_cow_oldadult[1] = 997)
  
  lambda_init_total ~ dgamma(11.1, 0.000863) # mean ≈ 12866 (sum of all three prev. values = 12866)
  
  # initial latent stage-specific abundances (from expected values)
  N_c[1] ~ dpois(lambda_init_c)
  N_y[1] ~ dpois(lambda_init_y)
  N_o[1] ~ dpois(lambda_init_o)
  
  N_total[1] ~ dpois(lambda_init_total)
  
  obs_total[1] ~ dlnorm(log(N_total[1] + 1e-6), tau_obs_total)
  
  # process model
  for (t in 1:(n_years-1)) {
    # expected values
    mu_y[t+1] <- s_c[t] * N_c[t] + s_y[t] * (1 - g_y[t]) * N_y[t]
    mu_o[t+1] <- s_y[t] * g_y[t] * N_y[t] + s_o[t] * N_o[t]
    
    # latent states
    N_y[t+1] ~ dpois(max(1e-9, mu_y[t+1]))
    N_o[t+1] ~ dpois(max(1e-9, mu_o[t+1]))
    
    # derive latent total abundance from latent stage-specific abundances 
    N_total[t+1] <- N_c[t+1] + N_y[t+1] + N_o[t+1]
    
    # observations of total abundance are related to latent total abundance
    obs_total[t+1] ~ dlnorm(log(N_total[t+1] + 1e-6), tau_obs_total)
  }
  
  ## -----------------------------
  ## (2) CJS MODEL
  ## -----------------------------
  
  # detection prob
  for (t in 1:n_years){
    p[t] ~ dunif(0, 1)
  }
  

  for (i in 1:N) {
    # Time 1
    z[i, 1] ~ dbern(equals(1, first_seen[i]))
    
    # Times 2 onward
    for (t in 2:n_years) {
      # Ensure at least one class is active (add small constant)
      phi[i, t] <- is_class1[i, t-1] * s_y[t-1] + 
        is_class2[i, t-1] * s_o[t-1] +
        1e-10  # prevent exactly 0
      
      z[i, t] ~ dbern(
        equals(t, first_seen[i]) +
          step(t - first_seen[i] - 0.5) * (1 - equals(t, first_seen[i])) *
          z[i, t-1] * phi[i, t]
      )
    }
    
    # observations
    for (t in 1:n_years) {
      y[i, t] ~ dbern(p[t] * z[i, t])
    }
  }
  
  ## -----------------------------
  ## (3) CALF SURVIVAL 
  ## -----------------------------
  
  for (t in 1:(n_years-1)) {
    mu_c[t+1] <- f_y[t] * s_c[t] * N_y[t] +
      f_o[t] * s_c[t] * N_o[t]
    
    N_c[t+1] ~ dpois(max(1e-6, mu_c[t+1]))
  }
  
  ## -----------------------------
  ## (4) FECUNDITY 
  ## -----------------------------

  for (t in 1:(n_years-1)){
    # priors
    f_y[t] ~ dbeta(1, 1)
    f_o[t] ~ dbeta(1, 1)
    
    # parameter estimation
    young_num_preg[t] ~ dbin(f_y[t], young_num_capt[t])
    old_num_preg[t] ~ dbin(f_o[t], old_num_capt[t])
  }
})

############################################################################################
########### --------------- MODEL SPECS, INITS, AND DATA SOURCES ---------------- ##########
############################################################################################
# constants & data 
N <- nrow(y)

elk_constants <- list(n_years = n_years, N = N)

elk_data <- list(
  # state-space
  obs_total = dat_n$n_female,  # female-only
  # CJS
  y = y,
  is_class1 = is_class1,
  is_class2 = is_class2,
  first_seen = first_seen,
  # fecundity
  young_num_preg = dat_fec$young_num_preg,
  young_num_capt = dat_fec$young_num_capt,
  old_num_preg = dat_fec$old_num_preg,
  old_num_capt = dat_fec$old_num_capt
)


## -----------------------------
## inits  
## -----------------------------

make_inits <- function() {
  # seed latent states near observations
  init_Nc <- ifelse(is.na(dat_n$n_calf), # is n_calf is NA for a given year,
                    pmax(1, round(mean(dat_n$n_calf, na.rm = TRUE))), # use the mean n_calf as initial, with na.rm
                    pmax(1, round(dat_n$n_calf))) # otherwise, use the "observed value"
  init_Ny <- ifelse(is.na(dat_n$n_cow_youngadult),
                    pmax(1, round(mean(dat_n$n_cow_youngadult, na.rm = TRUE))),
                    pmax(1, round(dat_n$n_cow_youngadult)))
  init_No <- ifelse(is.na(dat_n$n_cow_oldadult),
                    pmax(1, round(mean(dat_n$n_cow_oldadult, na.rm = TRUE))),
                    pmax(1, round(dat_n$n_cow_oldadult)))
  init_Ntotal <- ifelse(is.na(dat_n$n_total),
                    pmax(1, round(mean(dat_n$n_total, na.rm = TRUE))),
                    pmax(1, round(dat_n$n_total)))
  
  z_init <- matrix(NA, nrow = nrow(y), ncol = ncol(y))
  for (i in 1:nrow(y)) {
    detections <- which(y[i, ] == 1)
    if (length(detections) > 0) {
      first_det <- min(detections)
      last_det <- max(detections)
      # Set z=1 from first detection through last detection
      z_init[i, first_det:last_det] <- 1L
      # Everything before first detection stays NA (or set to 0)
      if (first_det > 1) {
        z_init[i, 1:(first_det-1)] <- 0L
      }
    }
  }
  
  list(
    # intercepts centered on beliefs
    alpha_sc = qlogis(0.22),
    alpha_sy = qlogis(0.90),
    alpha_so = qlogis(0.85),
    alpha_gy = qlogis(0.15),
    
    # SDs (start reasonably)
    sigma_sc = 0.20, 
    sigma_sy = 0.15, 
    sigma_so = 0.15,
    sigma_gy = 0.15, 
    
    # standard-normal year effects start at 0
    eps_sc_std = rep(0, n_years),
    eps_sy_std = rep(0, n_years),
    eps_so_std = rep(0, n_years),
    eps_gy_std = rep(0, n_years),
    
    # observation SDs
    sigma_obs_total = 0.30,
    
    # initial abundances
    lambda_init_c = max(1, round(init_Nc[1])),
    lambda_init_y = max(1, round(init_Ny[1])),
    lambda_init_o = max(1, round(init_No[1])),
    lambda_init_total = max(1, round(init_Ntotal[1])),
    N_c = pmax(1, init_Nc),
    N_y = pmax(1, init_Ny),
    N_o = pmax(1, init_No),
    N_total = pmax(1, init_Ntotal),
    
    # detection probability for cjs model
    p = runif(n_years, 0.6, 0.95),
    
    # initial z matrix for cjs model
    z = z_init,
    
    # fecundity initials
    f_y = rep(0.76, n_years - 1),  # mean(dat_fec$young_prop_pregnant, na.rm=TRUE) = 0.76
    f_o = rep(0.64, n_years - 1)  # mean(dat_fec$old_prop_pregnant, na.rm=TRUE) = 0.64
  )
}


## -----------------------------
## parameters to monitor
## -----------------------------
params <- c(
  # hyperparameters
  "alpha_sc","alpha_sy","alpha_so","alpha_gy",
  "sigma_sc","sigma_sy","sigma_so","sigma_gy",

  # yearly vital rates (shared by IPM & CJS)
  "s_c","s_y","s_o","g_y",
  "f_y", 'f_o',

  # detection (CJS)
  "p",
  
  # latent states (and total)
  "N_c","N_y","N_o","N_total"
)

## -----------------------------
## model specs
## -----------------------------
set.seed(17)
nc <- 3
ni <- 1000000
nb <- 200000
th <- 4

############################################################################################
####################### --------------- RUN MODEL ---------------- #########################
############################################################################################

# run MCMC
elk_mod1 <- nimbleMCMC(
  code      = elk_ipm,
  data      = elk_data,
  constants = elk_constants,
  inits     = make_inits,
  monitors  = params,
  nchains   = nc,
  niter     = ni,
  nburnin   = nb,
  thin      = th,
  summary   = TRUE
)

# OR IMPORT PREVIOUSLY RUN MODEL TO WORK WITH RESULTS BEYOND HERE
# elk_mod1 <- readRDS('data/results/elk_mod1_results.rds')

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
N_summ <- MCMCsummary(ml_clean, params=c('N_c', 'N_y', 'N_o', 'N_total'))   # extracts mean, sd, 2.5%, 50%, 97.5%

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
                        "N_total" = "Total"))


dat_long <- dat_n %>%
  pivot_longer(
    cols = -c(X,year),          # everything except year becomes stage/value pairs
    names_to = "stage",
    values_to = "value"
  ) %>%
  mutate(stage = recode(stage,
                        n_calf  = "Calf",
                        n_cow_youngadult = "Young Adult",
                        n_cow_oldadult   = "Old Adult",
                        n_cow = "Total"))



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
vrates$year <- rep(shared_years, 4)

ggplot(vrates, aes(x = year, y = mean)) +
  geom_ribbon(aes(ymin = low, ymax = high), alpha = 0.2) +
  geom_line(size = 0.9) +
  facet_wrap(~ rate, scales = "free_y") +
  theme_bw() +
  labs(x = "Year", y = "Estimated value",
       title = "Posterior Time-Varying Vital Rates (95% Credible Intervals)")

