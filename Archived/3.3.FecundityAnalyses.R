### Elk Abundance - Fecundity Model
### Analysis script
### Last updated: Oct. 29, 2025
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

dat <- read.csv('data/intermediate/productivity.csv')



#
#
# ### from Chat:
#
#

# -------------------------------------------------------------------
# Standalone fecundity model (pregnancy → births → late-winter calves)
# -------------------------------------------------------------------
# Data you can provide (choose what you have; missing bits can be omitted):
#  - preg_y_succ[t], preg_y_tested[t]   # young-adult pregnancy successes / tested (Dec–Feb of year t)
#  - preg_o_succ[t], preg_o_tested[t]   # old-adult pregnancy successes / tested (Dec–Feb of year t)
#  - N_y[t], N_o[t]                     # known/estimated numbers of young/old cows in year t (late winter)
#  - C_obs[t]                           # optional: observed late-winter calf count in year t (t >= 2)
#  - r_obs[t]                           # optional: observed CCR in year t as a proportion (calves per cow)
# Constants:
#  - n_years                            # number of late-winter time points (t = 1..n_years)

fec_code <- nimbleCode({
  
  ## -----------------------------
  ## Pregnancy by class & year
  ## -----------------------------
  # Intercepts on logit scale (weakly informative)
  alpha_p_y ~ dnorm(0, 0.01)                          # young pregnancy
  alpha_p_o ~ dnorm(0, 0.01)                          # old pregnancy
  
  # Year-to-year variability (half-normal via truncation)
  sigma_p_y ~ T(dnorm(0, 1/0.5^2), 0, )
  sigma_p_o ~ T(dnorm(0, 1/0.5^2), 0, )
  
  for (t in 1:n_years) {
    eps_p_y_std[t] ~ dnorm(0, 1)
    eps_p_o_std[t] ~ dnorm(0, 1)
    
    logit(pi_y[t]) <- alpha_p_y + sigma_p_y * eps_p_y_std[t]
    logit(pi_o[t]) <- alpha_p_o + sigma_p_o * eps_p_o_std[t]
    
    # If pregnancy test data exist, include these (otherwise, omit/comment them):
    # preg_y_succ[t] ~ dbin(pi_y[t], preg_y_tested[t])
    # preg_o_succ[t] ~ dbin(pi_o[t], preg_o_tested[t])
  }
  
  ## -----------------------------
  ## Female fraction at birth
  ## -----------------------------
  rho <- 0.5    # fix at 0.5 to start (can be estimated later if needed)
  
  ## -----------------------------
  ## Birth → late-winter calf survival (cohort born in t, seen at t+1)
  ## -----------------------------
  alpha_bw ~ dnorm(qlogis(0.6), 1/0.6^2)              # prior mean ~ 0.6
  sigma_bw ~ T(dnorm(0, 1/0.4^2), 0, )
  
  # Indexing: phi_bw[t] applies to cohort born in (t-1) and observed at t
  # so define for t = 2..n_years
  for (t in 2:n_years) {
    eps_bw_std[t] ~ dnorm(0, 1)
    logit(phi_bw[t]) <- alpha_bw + sigma_bw * eps_bw_std[t]
  }
  # filler for t=1 so it's defined; not used in the process
  phi_bw[1] <- ilogit(alpha_bw)
  
  ## -----------------------------
  ## Latent calf state and process
  ## -----------------------------
  # N_c[t] is late-winter female calves present at time t
  # For t >= 2, they come from pregnancies in year (t-1)
  # Assumption: every pregnant cow gives birth (no abortion/twinning)
  
  # Optionally, put a diffuse prior on the first year calf state (or tie it to data if available)
  N_c[1] ~ dpois(1.0)   # harmless weak prior; adjust if you have C_obs[1]
  
  for (t in 1:(n_years - 1)) {
    # Female births in summer of year t:
    births_fem[t] <- rho * (pi_y[t] * N_y[t] + pi_o[t] * N_o[t])
    
    # Late-winter calves at t+1:
    mu_c[t+1] <- phi_bw[t+1] * births_fem[t]
    
    # Process model:
    N_c[t+1] ~ dpois(max(1e-9, mu_c[t+1]))
  }
  
  ## -----------------------------
  ## Observation models (choose what you have)
  ## -----------------------------
  
  # A) Calf counts (late winter): Poisson around N_c
  # If you have counts C_obs[t], uncomment:
  # for (t in 1:n_years) {
  #   C_obs[t] ~ dpois(max(1e-9, N_c[t]))
  # }
  
  # B) Cow:calf ratio (CCR) as a Beta around N_c / (N_y + N_o)
  # Provide r_obs[t] in (0,1); if "per 100 cows", divide by 100 before passing in.
  # kappa_ccr controls dispersion (larger = tighter around the mean)
  # If you have CCR, uncomment this block:
  # kappa_ccr ~ dgamma(2, 2)
  # for (t in 1:n_years) {
  #   denom[t]  <- max(1e-6, N_y[t] + N_o[t])
  #   r_true[t] <- N_c[t] / denom[t]
  #   a_ccr[t]  <- max(1e-6, r_true[t] * kappa_ccr)
  #   b_ccr[t]  <- max(1e-6, (1 - r_true[t]) * kappa_ccr)
  #   r_obs[t]  ~ dbeta(a_ccr[t], b_ccr[t])
  # }
})