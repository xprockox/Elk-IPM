### CJS Model for Elk Data
### Last updated: Sept. 9, 2025
### Contact: xprockox@gmail.com

### --------------- PACKAGES ---------------- ###
library(dplyr)
library(lubridate)
library(tidyr)

### --------------- DATA IMPORT ---------------- ###
df <- read.csv('data/elk_survival_testData.csv')
df <- df[which(complete.cases(df)==TRUE),]

### --------------- DATA MANAGEMENT ---------------- ###

# clean df
df_clean <- df %>%
  mutate(
    BirthDate = ymd(BirthDate),
    Capture.Date = ymd(Capture.Date),
    Last.Date.Alive = ymd(Last.Date.Alive),
    Last.Date.Status = trimws(Last.Date.Status),
    BirthYear = year(BirthDate)
  ) %>%
  filter(!is.na(Capture.Date) & !is.na(Last.Date.Alive))  # Remove incomplete entries

# define study period and sample size
min_year <- min(df_clean$BirthYear, na.rm = TRUE)
max_year <- max(year(df_clean$Last.Date.Alive), na.rm = TRUE)
years <- min_year:max_year
n_years <- length(years)
n_indiv <- nrow(df_clean)

# intialize matrices
y <- matrix(0, nrow = n_indiv, ncol = n_years)
z <- matrix(NA, nrow = n_indiv, ncol = n_years)
colnames(y) <- colnames(z) <- years
rownames(y) <- rownames(z) <- df_clean$ID

# fill matrices
for (i in 1:n_indiv) {
  elk <- df_clean[i, ]
  birth_year <- elk$BirthYear
  capture_year <- year(elk$Capture.Date)
  last_year <- year(elk$Last.Date.Alive)
  status <- elk$Last.Date.Status
  
  # find column indices
  birth_idx <- which(years == birth_year)
  capture_idx <- which(years == capture_year)
  last_idx <- which(years == last_year)
  
  # latent state: assume alive (1) from birth until year of last confirmed alive
  if (!is.na(birth_idx) && !is.na(last_idx)) {
    z[i, birth_idx:last_idx] <- 1
  }
  
  # if the individual was confirmed dead, mark all future years as 0
  if (status == "Dead" && last_idx < n_years) {
    z[i, (last_idx + 1):n_years] <- 0
  }
  
  # detection matrix:
  if (status == "Dead" && capture_idx < last_idx) {
    y[i, capture_idx:(last_idx - 1)] <- 1
    y[i, last_idx] <- 0  # likely died that year
  } else if (status == "Live" && capture_idx <= last_idx) {
    y[i, capture_idx:last_idx] <- 1  # censored after this
  } else if (capture_idx == last_idx) {
    y[i, capture_idx] <- 1  # one-year animal
  }
}

# visual check of matrices
image(y, main = "Detection Matrix (y)", col = c("white", "black"))
image(z, main = "Latent State Matrix (z)", col = c("grey", "black"))

# vector of when the elk were first seen
first_seen <- apply(y, 1, function(row) {
  first <- which(row == 1)[1]
  if (is.na(first)) return(NA) else return(first)
})

### --------------- NIMBLE CODE ---------------- ###

elk_survival_code <- nimbleCode({
  # Priors
  phi ~ dunif(0, 1) # Survival probability
  p ~ dunif(0, 1) # Detection probability
  
  for (i in 1:N) {
    z[i, first_seen[i]] <- 1 # All animals alive at first seen
    
    for (t in (first_seen[i] + 1):n_years) {
      # Only proceed if we’re not past their observed window
      z[i, t] ~ dbern(z[i, t - 1] * phi)
    }
    
    for (t in first_seen[i]:n_years) {
      y[i, t] ~ dbern(p * z[i, t])
    }
  }
})

### --------------- CONSTANTS AND SPECS ---------------- ###

data <- list(
  y = y
)

constants <- list(
  N = nrow(y),
  n_years = ncol(y),
  first_seen = first_seen
)

# create z init matrix
z_init <- z
z_init[is.na(z_init)] <- 1  # assume alive unless evidence says otherwise

inits <- list(
  phi = runif(1, 0.8, 0.99),
  p = runif(1, 0.6, 0.95),
  z = z_init
)

params <- c("phi", "p")

### --------------- RUN MODEL ---------------- ###

elk_model <- nimbleMCMC(
  code = elk_survival_code,
  constants = constants,
  data = data,
  inits = inits,
  monitors = params,
  nchains = 3,
  niter = 10000,
  nburnin = 3000,
  thin = 1,
  summary = TRUE
)

elk_model
