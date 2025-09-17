### Elk Integrated Population Model
### Analysis script
### Last updated: June 2, 2025
### Contact: xprockox@gmail.com

############################################################################################
### packages
library(nimble)
library(MCMCvis)

############################################################################################
### data loading
modeling_df <- read.csv('data/intermediate/modeling_df.csv')

# let's try removing the last five years because they are all NAs
modeling_df <- modeling_df[1:26,]

############################################################################################
### model structure

elk_ipm1 <- nimbleCode({
  
  # priors for initial stage-specific abundances
  N_calf[1] ~ dpois(mu_calf_1)
  N_ya[1] ~ dpois(mu_ya_1)
  N_oa[1] ~ dpois(mu_oa_1)
  
  # priors for time-varying fecundity and survival
  for (t in 1:(nYears - 1)) {
    f[t] ~ T(dnorm(fem_birth[t], 0.25), 0, ) # fecundity with truncation
    
    phi_calf[t] ~ T(dnorm(surv_calf[t], 0.5), 0, 1) # calf survival
    phi_ya[t] ~ T(dnorm(surv_ya[t], 0.5), 0, 1) # young adult survival
    phi_oa[t] ~ T(dnorm(surv_oa[t], 0.5), 0, 1) # old adult survival
    
    psi[t] ~ dbeta(1, 1) # transition prob. from ya to oa
  }
  
  # obs error
  sigma_obs ~ dunif(0, 25)
  
  # process model
  for (t in 1:(nYears - 1)) {
    N_calf[t+1] ~ dpois(f[t] * (N_ya[t] + N_oa[t]))
    
    N_ya_from_calf[t] ~ dbin(phi_calf[t], N_calf[t])
    N_ya_survive[t] ~ dbin(phi_ya[t], N_ya[t])
    N_ya_to_oa[t] ~ dbin(psi[t], N_ya_survive[t])
    N_ya[t+1] <- N_ya_survive[t] - N_ya_to_oa[t] + N_ya_from_calf[t]
    
    N_oa_survive[t] ~ dbin(phi_oa[t], N_oa[t])
    N_oa[t+1] <- N_oa_survive[t] + N_ya_to_oa[t]
  }
  
  # Derived total abundance
  for (t in 1:nYears) {
    N_tot[t] <- N_calf[t] + N_ya[t] + N_oa[t]
  }
  
  # Observation model
  for (t in 1:nYears) {
    y[t] ~ dnorm(N_tot[t], sd = sigma_obs)
  }
  
})

############################################################################################
### model set-up and execution

# model specs
ni = 30000 # iterations
nb = 3000 # burn-in
nc = 3 # chains

# constants
nYears <- nrow(modeling_df) # define here first so it can be used in inits

# informative initial values based on 1995 proportions
N_calf_1 <- modeling_df$Total_Elk_Female[1] * modeling_df$Percent.N.calves[1]
N_ya_1 <- modeling_df$Total_Elk_Female[1] * modeling_df$Percent.N.prime[1]
N_oa_1 <- modeling_df$Total_Elk_Female[1] * modeling_df$Percent.N.old[1]

# data
data <- list(
  y = modeling_df$Total_Elk_Female,
  fem_birth = modeling_df$FemBirthPerCow[1:(nYears - 1)],
  surv_calf = modeling_df$Survival_calf[1:(nYears - 1)],
  surv_ya = modeling_df$Survival_prime[1:(nYears - 1)],
  surv_oa = modeling_df$Survival_old[1:(nYears - 1)]
)

# constants
constants <- list(
  nYears = nYears,
  mu_calf_1 = N_calf_1,
  mu_ya_1 = N_ya_1,
  mu_oa_1 = N_oa_1
)

# initial values
inits <- list(
  # abundances
  N_calf = pmax(1, round(rnorm(nYears, mean = N_calf_1, sd = N_calf_1 / 5))),
  N_ya = pmax(1, round(rnorm(nYears, mean = N_ya_1, sd = N_ya_1 / 5))),
  N_oa = pmax(1, round(rnorm(nYears, mean = N_oa_1, sd = N_oa_1 / 5))),
  
  # vital rates
  f = pmax(0.01, rnorm(nYears - 1, mean = mean(data$fem_birth), sd = 0.05)),
  phi_calf = pmin(pmax(rnorm(nYears - 1, mean = mean(data$surv_calf), sd = 0.05), 0.01), 0.99),
  phi_ya = pmin(pmax(rnorm(nYears - 1, mean = mean(data$surv_ya), sd = 0.05), 0.01), 0.99),
  phi_oa = pmin(pmax(rnorm(nYears - 1, mean = mean(data$surv_oa), sd = 0.05), 0.01), 0.99),
  psi = rep(0.08, nYears - 1),
  
  # obs error
  sigma_obs = round(runif(1, 0, 10))
)

# parameters to estimate
params <- c(
  "phi_calf", "phi_ya", "phi_oa",
  "psi", "sigma_obs",
  "f",
  "N_calf", "N_ya", "N_oa", "N_tot"
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
stop()
write.csv(mod1_results_long, 'data/results/mod1_results_long.csv')

############################################################################################
### trace plots

MCMCtrace(elk_mod1,
          params = c("f", "phi_calf", "phi_ya", "phi_oa", "psi", "sigma_obs"),
          ind = TRUE, 
          pdf = FALSE)
