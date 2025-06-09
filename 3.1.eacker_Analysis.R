### Elk Integrated Population Model Adapted from Eacker et al. 2017
### Last updated: June 9, 2025
### Contact: xprockox@gmail.com

############################################################################################
### packages
library(nimble)
library(MCMCvis)

############################################################################################
### data loading and manipulation
modeling_df <- read.csv('data/intermediate/modeling_df.csv')

# Approximate adult female counts (prime + old)
modeling_df$log.AF <- log(
  modeling_df$Total_Elk_Female * 
    (modeling_df$Percent.N.prime + modeling_df$Percent.N.old) + 0.000001
)

# Total counts
modeling_df$log.Y <- log(modeling_df$Total_Elk_Female + 0.000001)

############################################################################################
### model structure
eaker_ipm <- nimbleCode({
  
  # ------------------------
  # Priors
  # ------------------------
  l.tauY ~ dgamma(0.001, 0.001)
  l.tauAF ~ dgamma(0.001, 0.001)
  
  psi ~ dunif(0, 1)  # transition from young to old adult
  Sy ~ dbeta(21.92, 2.90)  # informative prior on young adult survival
  
  for (t in 1:(nYears - 1)) {
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
  
  # ------------------------
  # Initial conditions
  # ------------------------
  N_calf[1] ~ dpois(Ninit_calf)
  N_ya[1] ~ dpois(Ninit_ya)
  N_oa[1] ~ dpois(Ninit_oa)
  
  # ------------------------
  # Process model
  # ------------------------
  for (t in 1:(nYears - 1)) {
    expected_calf_recruits[t] <- f_ya[t] * N_ya[t] + f_oa[t] * N_oa[t]
    N_calf[t + 1] ~ dpois(expected_calf_recruits[t])
    
    N_ya_from_calf[t] ~ dbin(phi_calf[t], N_calf[t])
    
    N_ya_survive[t] ~ dbin(phi_ya[t], N_ya[t])
    N_ya_to_oa[t] ~ dbin(psi, N_ya_survive[t])
    N_ya[t + 1] <- N_ya_survive[t] - N_ya_to_oa[t] + N_ya_from_calf[t]
    
    N_oa_survive[t] ~ dbin(phi_oa[t], N_oa[t])
    N_oa[t + 1] <- N_oa_survive[t] + N_ya_to_oa[t]
  }
  
  # ------------------------
  # Observation model
  # ------------------------
  for (t in 1:nYears) {
    N_total[t] <- N_calf[t] + N_ya[t] + N_oa[t]
    log.Y[t] ~ dnorm(log(N_total[t] + 0.000001), l.tauY)
    log.AF[t] ~ dnorm(log(N_ya[t] + N_oa[t] + 0.000001), l.tauAF)
  }
  
  # ------------------------
  # Derived growth rates
  # ------------------------
  for (t in 1:(nYears - 1)) {
    pop_growth[t] <- (N_total[t + 1] + 0.000001) / (N_total[t] + 0.000001)
  }
})

############################################################################################
### model setup
ni <- 30000
nb <- 3000
nc <- 3

nYears <- nrow(modeling_df)

# initial abundances based on 1995 proportions
N_calf_1 <- modeling_df$Total_Elk_Female[1] * modeling_df$Percent.N.calves[1]
N_ya_1 <- modeling_df$Total_Elk_Female[1] * modeling_df$Percent.N.prime[1]
N_oa_1 <- modeling_df$Total_Elk_Female[1] * modeling_df$Percent.N.old[1]

# data
data <- list(
  log.Y = modeling_df$log.Y,
  log.AF = modeling_df$log.AF
)

constants <- list(
  nYears = nYears,
  Ninit_calf = N_calf_1,
  Ninit_ya = N_ya_1,
  Ninit_oa = N_oa_1
)

inits <- list(
  N_calf = pmax(1, round(rnorm(nYears, mean = N_calf_1, sd = N_calf_1 / 5))),
  N_ya = pmax(1, round(rnorm(nYears, mean = N_ya_1, sd = N_ya_1 / 5))),
  N_oa = pmax(1, round(rnorm(nYears, mean = N_oa_1, sd = N_oa_1 / 5))),
  
  log_f_ya = rnorm(nYears - 1, log(0.3), 0.2),
  log_f_oa = rnorm(nYears - 1, log(0.1), 0.2),
  
  logit_phi_calf = rnorm(nYears - 1, 0, 1.5),
  logit_phi_ya = rnorm(nYears - 1, 0, 1.5),
  logit_phi_oa = rnorm(nYears - 1, 0, 1.5),
  
  psi = 0.08,
  l.tauY = 1,
  l.tauAF = 1
)

params <- c("log_f_ya", "log_f_oa", "phi_calf", "phi_ya", "phi_oa", 
            "N_calf", "N_ya", "N_oa", "pop_growth", "psi", "l.tauY", "l.tauAF")

elk_mod1 <- nimbleMCMC(
  code = eacker_ipm,
  data = data,
  constants = constants,
  inits = inits,
  monitors = params,
  nchains = nc,
  niter = ni,
  nburnin = nb
)

mod1_results_long <- MCMCsummary(elk_mod1, params = 'all', round = 3)
stop()
write.csv(mod1_results_long, 'data/results/mod1_results_long.csv')

############################################################################################
### trace plots

MCMCtrace(elk_mod1,
          params = c("f", "phi_calf", "phi_ya", "phi_oa", "psi", "sigma_obs"),
          ind = TRUE, 
          pdf = FALSE)
