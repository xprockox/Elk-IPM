### Elk IPM Main Script
### Last updated: Oct. 17, 2025
### Contact: xprockox@gmail.com

### --------------- PACKAGES AND SET-UP ---------------- ###
library(dplyr)
library(lubridate)
library(tidyr)
library(tidyverse)
library(MCMCvis)
library(ggplot2)
library(nimble)

# the following line determines whether your latent state and observation
# matrices begin at the first year an elk was born (1982) or at the 
# first year an elk was captured (2000).

# the logic here is that there is no detection data from prior to capture,
# so estimates might be impacted by the mismatch between detection data and 
# latent state data (i.e. detection probability would be zero up until capture)

# if use_clip == TRUE, models are constructed using data from 2000:present
# if use_clip == FALSE, models are constructed using the full dataset (1982:present)
use_clip <- TRUE 

### --------------- DATA IMPORT ---------------- ###
df <- read.csv('data/master/elk_survival_2025-10-14.csv')
df <- df[which(complete.cases(df)==TRUE),] # drop rows with any NAs

vhf <- read.csv('data/vhf.csv')
gps <- readRDS('data/res.rds')

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
  filter(!is.na(Capture.Date) & !is.na(Last.Date.Alive)) %>% 
  filter(!is.na(BirthYear)) %>% # Remove incomplete entries
  arrange(ID, Capture.Date)

# clean VHF and GPS data
vhf <- vhf %>%
  mutate(
    Date = as.POSIXct(Date, format = "%m/%d/%Y"),
    Signal.Heard = case_when(
      is.na(Signal.Heard) ~ "YES",
      Signal.Heard == "" ~ "YES",
      Signal.Heard %in% c("Yes", "YES") ~ "YES",
      Signal.Heard %in% c("No", "NO") ~ "NO",
      TRUE ~ as.character(Signal.Heard)
    ),
    Visual. = case_when(
      Visual. %in% c("YES", "Yes", "yes", "Y") ~ "YES",
      TRUE ~ "NO"
    ),
    Year = year(Date)
  ) %>%
  filter(Signal.Heard == "YES" | Visual. == "YES") %>%
  rename(Visual = Visual.)

gps <- gps %>%
  mutate(
    dt = as.POSIXct(dt, format = "%Y-%m-%d %H:%M:%S"),
  ) %>%
  left_join(
    df_clean %>% select(ID, Last.Date.Alive), 
    by = "ID"
  ) %>%
  mutate(
    Year = year(dt),
    LastYear = year(as.POSIXct(Last.Date.Alive, format = "%Y-%m-%d"))
  ) %>%
  filter(Year <= LastYear)  # Remove GPS points after animal's death

# define study period and sample size
min_year <- min(df_clean$BirthYear, na.rm = TRUE)
max_year <- max(year(df_clean$Last.Date.Alive), na.rm = TRUE)
years <- min_year:max_year
n_years <- length(years)
n_indiv <- nrow(df_clean)

### --------------- MATRIX CONSTRUCTION (OBS. & LATENT STATES) ---------------- ###

# Initialize matrices
y <- matrix(0, nrow = n_indiv, ncol = n_years) # observations
z <- matrix(NA, nrow = n_indiv, ncol = n_years) # latent (true) states
colnames(y) <- colnames(z) <- years
rownames(y) <- rownames(z) <- df_clean$ID

# Fill in matrices
for (i in 1:n_indiv) {
  elk <- df_clean[i, ]
  elk_id <- elk$ID
  
  birth_year <- elk$BirthYear
  capture_year <- year(elk$Capture.Date)
  last_year <- year(elk$Last.Date.Alive)
  status <- elk$Last.Date.Status
  
  # Column indices
  birth_idx <- which(years == birth_year)
  capture_idx <- which(years == capture_year)
  last_idx <- which(years == last_year)
  
  # Latent state: assume alive from birth through year of last known alive
  if (!is.na(birth_idx) && !is.na(last_idx)) {
    z[i, birth_idx:last_idx] <- 1
  }
  
  # If known dead, mark all years after death as 0
  if (status == "Dead" && last_idx < n_years) {
    z[i, (last_idx + 1):n_years] <- 0
  }
  
  # Observation matrix (y) constructed based on VHF and GPS data
  
  # Extract observed years from VHF
  vhf_years <- vhf %>%
    filter(ID == elk_id) %>%
    pull(Year) %>%
    unique()
  
  # Extract observed years from GPS
  gps_years <- gps %>%
    filter(ID == elk_id) %>%
    pull(Year) %>%
    unique()
  
  # Combine and deduplicate
  observed_years <- sort(unique(c(vhf_years, gps_years)))
  
  # Fill y matrix
  for (yr in observed_years) {
    if (yr %in% years) {
      y[i, as.character(yr)] <- 1
    }
  }
}

# The above process does not necessarily mark the year of capture as a detection, 
# e.g. if the collar failed immediately and provided no GPS or VHF data,
# so we can force detections for all individuals during their year of capture:
df_clean$CaptureYear <- lubridate::year(df_clean$Capture.Date)

for (i in 1:nrow(df_clean)) {
  cap_year <- df_clean$CaptureYear[i]
  cap_idx <- which(years == cap_year)
  
  if (length(cap_idx) == 1) {
    y[i, cap_idx] <- 1
  }
  
  if (!is.na(cap_idx) && cap_idx > 1) {
    y[i, 1:(cap_idx - 1)] <- 0
  }
}

# Now if the latent state (true state) is dead, there can't be an observation.
# Mask impossible observations in y: set to NA if z = 0
for (i in 1:n_indiv) {
  for (t in 1:n_years) {
    if (!is.na(z[i, t]) && z[i, t] == 0 && y[i, t] == 1) {
      message(paste("Setting y[", i, ",", t, "] to NA due to z = 0 but y = 1"))
      y[i, t] <- NA  # mask inconsistent observations
    }
    
    # Optional: also mask y after last known alive (if not already done)
    if (!is.na(z[i, t]) && z[i, t] == 0) {
      y[i, t] <- NA  # cannot be observed if truly dead
    }
  }
}

# And add a final check:

# If z = 0, we should not see an observation (mask it)
y[z == 0] <- 0

# Ensure that wherever y = 1, z = 1
z[y == 1] <- 1

# create a new latent state matrix that does not include information about births 
# (i.e. the individual is NA up until capture despite being alive)
# (i.e.i.e we're pretending we don't have cementum info)
z_clipped <- z  # Copy original latent state matrix

for (i in 1:n_indiv) {
  cap_year <- df_clean$CaptureYear[i]
  cap_idx <- which(years == cap_year)
  
  if (length(cap_idx) == 1) {
    z_clipped[i, 1:(cap_idx - 1)] <- NA  # Mask all years before capture
  }
}

# drop years before capture
z_clipped <- z_clipped[,19:ncol(z_clipped)]

# now do the same for the observation matrix
y_clipped <- y[,19:ncol(y)]

### --------------- VISUALIZE MATRICES ---------------- ###

### first observations:
# Flip y so individual 1 is at the top
y_flip <- y[nrow(y):1, ]

# Keep only columns for years >= 2000
# cols_to_keep <- which(as.numeric(colnames(y_flip)) >= 2000)
# y_flip <- y_flip[, cols_to_keep]

# Plot detection matrix
image(
  x = 1:ncol(y_flip),
  y = 1:nrow(y_flip),
  z = t(y_flip),  # transpose for image()
  col = c("white", "black"),
  axes = FALSE,
  xlab = "Year",
  ylab = "Individual",
  main = "Detection Matrix (y)"
)

# Add axis labels
axis(1,
     at = 1:ncol(y_flip),
     labels = colnames(y_flip),
     las = 2,
     cex.axis = 0.7)

axis(2,
     at = 1:nrow(y_flip),
     labels = rev(rownames(y)),  # reversed to match y_flip
     las = 1,
     cex.axis = 0.4)

### then latent states:
# Flip z so that individual 1 is at the top
z_flip <- z[nrow(z):1, ]

# Set up plot
image(
  x = 1:ncol(z_flip),
  y = 1:nrow(z_flip),
  z = t(z_flip),  # transpose for correct orientation
  col = c("white", "gray", "black"),  # NA = white, 0 = gray, 1 = black
  breaks = c(-0.1, 0.1, 0.9, 1.1),    # assign color by value (NA, 0, 1)
  axes = FALSE,
  xlab = "Year",
  ylab = "Individual",
  main = "Latent State Matrix (z)"
)

# Add axis labels
axis(1, at = 1:ncol(z_flip), labels = colnames(z_flip), las = 2, cex.axis = 0.7)
axis(2, at = 1:nrow(z_flip), labels = rownames(z_flip), las = 1, cex.axis = 0.4)

rm(z_flip, y_flip)

#
##
###
####
#####
######
### ------- then plot latent states CLIPPED:
######
#####
####
###
##
#

### start with y_clipped
# Flip y_clipped so that individual 1 is at the top
y_clipped_flip <- y_clipped[nrow(y_clipped):1, ]

# Set up plot
image(
  x = 1:ncol(y_clipped_flip),
  y = 1:nrow(y_clipped_flip),
  z = t(y_clipped_flip),  # transpose for correct orientation
  col = c("white", "gray", "black"),  # NA = white, 0 = gray, 1 = black
  breaks = c(-0.1, 0.1, 0.9, 1.1),    # assign color by value (NA, 0, 1)
  axes = FALSE,
  xlab = "Year",
  ylab = "Individual",
  main = "Detection Matrix (y)"
)

# Add axis labels
axis(1, at = 1:ncol(y_clipped_flip), labels = colnames(y_clipped_flip), las = 2, cex.axis = 0.7)
axis(2, at = 1:nrow(y_clipped_flip), labels = rownames(y_clipped_flip), las = 1, cex.axis = 0.4)

### then plot z_clipped
# Flip z_clipped so that individual 1 is at the top
z_clipped_flip <- z_clipped[nrow(z_clipped):1, ]

# Set up plot
image(
  x = 1:ncol(z_clipped_flip),
  y = 1:nrow(z_clipped_flip),
  z = t(z_clipped_flip),  # transpose for correct orientation
  col = c("white", "gray", "black"),  # NA = white, 0 = gray, 1 = black
  breaks = c(-0.1, 0.1, 0.9, 1.1),    # assign color by value (NA, 0, 1)
  axes = FALSE,
  xlab = "Year",
  ylab = "Individual",
  main = "Latent State Matrix (z)"
)

# Add axis labels
axis(1, at = 1:ncol(z_clipped_flip), labels = colnames(z_clipped_flip), las = 2, cex.axis = 0.7)
axis(2, at = 1:nrow(z_clipped_flip), labels = rownames(z_clipped_flip), las = 1, cex.axis = 0.4)

rm(z_clipped_flip, y_clipped_flip)

### --------------- FIRST SEEN VECTOR ---------------- ###

# finally, we need a vector of when each individual was first seen
first_seen <- apply(y, 1, function(row) {
  first <- which(row == 1)[1]
  if (is.na(first)) return(n_years) else return(first)
})
first_seen <- as.integer(first_seen)

first_seen_clipped <- apply(y_clipped, 1, function(row) {
  first <- which(row == 1)[1]
  if (is.na(first)) return(ncol(y_clipped)) else return(first)
})
first_seen_clipped <- as.integer(first_seen_clipped)

stop('All required matrices constructed (line 338).')