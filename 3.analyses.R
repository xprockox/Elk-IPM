### Elk Integrated Population Model
### Last updated: May 29, 2025
### Contact: xprockox@gmail.com

#######################################################################################################
### packages
library(nimble)
library(MCMCvis)

#######################################################################################################
### data loading

#######################################################################################################
### model structure

elk_ipm1 <- nimbleCode({
  
  for (t in 1:(n_years - 1)) {
    logit_phi_calf[t] ~ dnorm(0, 1.5)
    logit_phi_ya[t] ~ dnorm(0, 1.5)
    logit_phi_oa[t] ~ dnorm(0, 1.5)
    
    log_f_ya[t] ~ dnorm(0, 1.5)
    log_f_oa[t] ~ dnorm(0, 1.5)
    
    phi_calf[t] <- ilogit(logit_phi_calf[t])
    phi_ya[t] <- ilogit(logit_phi_ya[t])
    phi_oa[t] <- ilogit(logit_phi_oa[t])
    f_ya[t] <- exp(log_f_ya[t])
    f_oa[t] <- exp(log_f_oa[t])
  }
  
  psi ~ dunif(0, 1)
  tau_obs ~ dgamma(0.001, 0.001)
  
  for (t in 1:(n_years - 1)) {
    expected_calf_recruits[t] <- f_ya[t] * N_ya[t] + f_oa[t] * N_oa[t]
    N_calf[t+1] ~ dpois(expected_calf_recruits[t])
    
    N_ya_from_calf[t] ~ dbin(phi_calf[t], N_calf[t])
    
    N_ya_survive[t] ~ dbin(phi_ya[t], N_ya[t])
    N_ya_to_oa[t] ~ dbin(psi, N_ya_survive[t])
    N_ya[t+1] <- N_ya_survive[t] - N_ya_to_oa[t] + N_ya_from_calf[t]
    
    N_oa_survive[t] ~ dbin(phi_oa[t], N_oa[t])
    N_oa[t+1] <- N_oa_survive[t] + N_ya_to_oa[t]
  }
  
  for (t in 1:n_years) {
    N_total[t] <- N_calf[t] + N_ya[t] + N_oa[t]
    y[t] ~ dnorm(N_total[t], tau_obs)
  }
  
  N_calf[1] ~ dpois(1000)
  N_ya[1] ~ dpois(5000)
  N_oa[1] ~ dpois(2000)
})

#######################################################################################################
### model set-up and execution

ni = 1000 # iterations
nb = 100 # burn-in
nc = 3 # chains

# data
data <- list(
  # calf_marked = c(...), 
  # calf_survived = c(...),
  # n_ya_marked = ..., 
  # n_oa_marked = ..., 
  # survived = ...
)

# constants
years <- 1995:2024 # adjust as needed
constants <- list(
  years = years,
  n_years = length(years)
)

# inits
inits <- list(
  N_calf = rep(1000, constants$n_years),
  N_ya = rep(5000, constants$n_years),
  N_oa = rep(2000, constants$n_years),
  
  logit_phi_calf = rnorm(constants$n_years - 1, 0, 1),
  logit_phi_ya = rnorm(constants$n_years - 1, 1.5, 1),
  logit_phi_oa = rnorm(constants$n_years - 1, 1.5, 1),
  
  log_f_ya = rnorm(constants$n_years - 1, log(0.3), 0.2),
  log_f_oa = rnorm(constants$n_years - 1, log(0.1), 0.2),
  
  psi = 0.08,
  tau_obs = 1
)


# parameters to estimate
params <- c(
  
)

# run the model
elk_mod1 <- nimbleMCMC(
  code = elk_ipm1,
  data = data,
  constants = constants,
  inits = inits,
  monitors = params,
  nchains = nc,
  niter = ni,
  nburnin = nb
)

# view model results
mod1_results_long <- MCMCsummary(elk_mod1, params = 'all', round = 3)
write.csv(mod1_results_long, 'data/results/mod1_results_long.csv')

#######################################################################################################