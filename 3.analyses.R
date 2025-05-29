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
### modeling

elk_ipm <- nimbleCode({
  
  # -----------------------------------------
  # 1. PRIORS
  # -----------------------------------------
  for (t in 1:(T-1)) {
    # survival (logit scale)
    logit_phi_calf[t] ~ dnorm(0, 1.5)
    logit_phi_ya[t]   ~ dnorm(0, 1.5)
    logit_phi_oa[t]   ~ dnorm(0, 1.5)
    
    # fecundity (log scale)
    log_f_ya[t] ~ dnorm(0, 1.5)
    log_f_oa[t] ~ dnorm(0, 1.5)
    
    # transform to probability / rate scale
    phi_calf[t] <- ilogit(logit_phi_calf[t])
    phi_ya[t]   <- ilogit(logit_phi_ya[t])
    phi_oa[t]   <- ilogit(logit_phi_oa[t])
    f_ya[t]     <- exp(log_f_ya[t])
    f_oa[t]     <- exp(log_f_oa[t])
  }
  
  psi ~ dunif(0, 1) # prop. of YAs aging to OAs
  tau_obs ~ dgamma(0.001, 0.001) # obs. error
  
  # -----------------------------------------
  # 2. STATE PROCESS MODEL
  # -----------------------------------------
  for (t in 1:(T-1)) {
    # calf recruitment
    expected_calf_recruits[t] <- f_ya[t] * N_ya[t] + f_oa[t] * N_oa[t]
    N_calf[t+1] ~ dpois(expected_calf_recruits[t])
    
    # calf survival to YA
    N_ya_from_calf[t] ~ dbin(phi_calf[t], N_calf[t])
    
    # YA survival and recruitment
    N_ya_survive[t] ~ dbin(phi_ya[t], N_ya[t])
    N_ya_to_oa[t] ~ dbin(psi, N_ya_survive[t])
    N_ya[t+1] <- N_ya_survive[t] - N_ya_to_oa[t] + N_ya_from_calf[t]
    
    # OA survival + YA recruited
    N_oa_survive[t] ~ dbin(phi_oa[t], N_oa[t])
    N_oa[t+1] <- N_oa_survive[t] + N_ya_to_oa[t]
  }
  
  # -----------------------------------------
  # 3. OBSERVATION MODEL
  # -----------------------------------------
  for (t in 1:T) {
    N_total[t] <- N_calf[t] + N_ya[t] + N_oa[t]
    y[t] ~ dnorm(N_total[t], tau_obs)
  }
  
  # -----------------------------------------
  # 4. INITIAL STATES
  # -----------------------------------------
  N_calf[1] ~ dpois(1000) # will need to update these once I see the data
  N_ya[1]   ~ dpois(5000)
  N_oa[1]   ~ dpois(2000)
})