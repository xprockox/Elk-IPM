### Constructing Survival Matrices
### Last updated: Oct. 27, 2025
### Contact: xprockox@gmail.com

############################################################
### -------------------- PACKAGES  --------------------- ###
############################################################
library(dplyr)
library(lubridate)
library(tidyr)
library(tidyverse)
library(ggplot2)

############################################################
### ------------------- DATA IMPORT -------------------- ###
############################################################
df <- read.csv('data/master/elk_survival_2025-10-14.csv')
df <- df[which(complete.cases(df)==TRUE),] # drop rows with any NAs

vhf <- read.csv('data/master/elk_vhf_2025-07-18.csv')
gps <- read.csv('data/master/elk_GPS_2025-09-03.csv')

############################################################
### ----------------- DATA MANAGEMENT ------------------ ###
############################################################

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

############################################################
### ---------------- MATRIX CONSTRUCTION --------------- ###
### -------------- (OBS. & LATENT STATES) -------------- ###
############################################################

# Initialize matrices
y <- matrix(NA, nrow = n_indiv, ncol = n_years) # observations
z <- matrix(NA, nrow = n_indiv, ncol = n_years) # latent (true) states
colnames(y) <- colnames(z) <- years
rownames(y) <- rownames(z) <- df_clean$ID

for (i in 1:n_indiv) {
  
  elk_id <- df_clean$ID[i]
  
  birth_year   <- df_clean$BirthYear[i]
  capture_year <- year(df_clean$Capture.Date[i])
  last_year    <- year(df_clean$Last.Date.Alive[i])
  status       <- df_clean$Last.Date.Status[i]
  
  # Convert to column indices
  birth_idx   <- match(birth_year, years)
  capture_idx <- match(capture_year, years)
  last_idx    <- match(last_year, years)
  
  ########################################################
  ### 1) LATENT STATE MATRIX (z)
  ########################################################
  
  # Alive from birth through last known alive year
  if (!is.na(birth_idx) && !is.na(last_idx)) {
    z[i, birth_idx:last_idx] <- 1
  }
  
  # After last known alive
  if (!is.na(last_idx) && last_idx < n_years) {
    if (status == "Dead") {
      z[i, (last_idx+1):n_years] <- 0
    } else {  # status == "Live"
      z[i, (last_idx+1):n_years] <- NA
    }
  }
  
  ########################################################
  ### 2) DETECTION MATRIX (y)
  ########################################################
  
  # Get detection years from BOTH sources
  gps_years <- gps$Year[gps$ID == elk_id]
  vhf_years <- vhf$Year[vhf$ID == elk_id]
  
  observed_years <- sort(unique(c(gps_years, vhf_years)))
  
  # Keep only years within alive window
  observed_years <- observed_years[
    observed_years >= capture_year &
      observed_years <= last_year
  ]
  
  # Convert to indices
  obs_idx <- match(observed_years, years)
  
  # 2a) Mark detections as 1
  if (length(obs_idx) > 0) {
    y[i, obs_idx] <- 1
  }
  
  # 2b) FORCE capture year to be a detection
  if (!is.na(capture_idx)) {
    y[i, capture_idx] <- 1
  }
  
  # 2c) Fill zeros BETWEEN capture and last known alive
  if (!is.na(capture_idx) && !is.na(last_idx)) {
    idx_window <- capture_idx:last_idx
    y[i, idx_window][is.na(y[i, idx_window])] <- 0
  }
  
  # 2d) Before capture = NA
  if (!is.na(capture_idx) && capture_idx > 1) {
    y[i, 1:(capture_idx-1)] <- NA
  }
  
  # 2e) After last known alive = NA
  if (!is.na(last_idx) && last_idx < n_years) {
    y[i, (last_idx+1):n_years] <- NA
  }
}

# create a new latent state matrix that does not include information about births 
# (i.e. the individual is NA up until capture despite being alive)
# (i.e.i.e we're pretending we don't have cementum info)
z_clipped <- z  # Copy original latent state matrix

# limit to study years
z_clipped <- z_clipped[,14:ncol(z_clipped)]

# now do the same for the observation matrix
y_clipped <- y[,14:ncol(y)]

############################################################
### ---------------- VISUALIZE MATRICES ---------------- ###
############################################################

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

# then plot latent states CLIPPED:

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

############################################################
### ---------------- FIRST SEEN VECTOR ----------------- ###
############################################################

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


############################################################
### --------------- BUILD STAGE MATRIX ----------------- ###
############################################################

# Copy z and mask non-alive values
age_class <- z
age_class[age_class != 1] <- NA

# Assign class 1 (0–14 y) and class 2 (> 14 y)
for (i in 1:nrow(age_class)) {
  alive_years <- which(age_class[i, ] == 1)
  if (length(alive_years) == 0) next
  split_idx <- min(length(alive_years), 15)  # 14 years of class 1, then switch
  age_class[i, alive_years[1:split_idx]] <- 1
  if (length(alive_years) > split_idx) {
    age_class[i, alive_years[(split_idx + 1):length(alive_years)]] <- 2
  }
}

# Get last seen year index
years_vector <- as.numeric(colnames(y))
last_seen_index <- match(year(as.Date(df$Last.Date.Alive)), years_vector)

# Create clipped version of the stage matrix

# copy stage matrix
age_class_clipped <- age_class

for (i in 1:n_indiv) {
  cap_year <- df_clean$CaptureYear[i]
  cap_idx <- which(years == cap_year)
  
  if (length(cap_idx) == 1) {
    age_class_clipped[i, 1:(cap_idx - 1)] <- NA  # Mask all years before capture
  }
}

# clip to years 2000:2024
age_class_clipped <- age_class_clipped[,14:ncol(age_class_clipped)]

# create dummy matrices that are 0s and 1s representing each class (helps with if-then
# logic in the model block, since NIMBLE doesn't support actual if-then statements)
is_class1_clipped <- array(0, dim = c(nrow(age_class_clipped), ncol(age_class_clipped)))
is_class2_clipped <- array(0, dim = c(nrow(age_class_clipped), ncol(age_class_clipped)))

is_class1_clipped[age_class_clipped == 1] <- 1
is_class2_clipped[age_class_clipped == 2] <- 1

############################################################
### ------------- VISUALIZE STAGE MATRIX --------------- ###
############################################################

# reset the visualizing pane (par settings from MCMC plots still active)
dev.off()

### non-clipped:
# Recode age_class for plotting
plot_matrix <- age_class
plot_matrix[plot_matrix == 1] <- 1  # class 1 (0–14 y)
plot_matrix[plot_matrix == 2] <- 2  # class 2 (>14 y)
plot_matrix[is.na(plot_matrix)] <- 0  # non-alive = white

# Flip
plot_matrix <- plot_matrix[nrow(plot_matrix):1, ]

# Define colors: white = NA, blue = class 1, orange = class 2
plot_colors <- c("white", "skyblue", "orange")

# Plot
image(
  x = 1:ncol(plot_matrix),
  y = 1:nrow(plot_matrix),
  z = t(plot_matrix),
  col = plot_colors,
  axes = FALSE,
  xlab = "Year",
  ylab = "Individual",
  main = "Age Class Matrix"
)

# Axis labels
axis(1, at = 1:ncol(plot_matrix), labels = colnames(plot_matrix), las = 2, cex.axis = 0.7)
axis(2, at = 1:nrow(plot_matrix), labels = rev(rownames(plot_matrix)), las = 1, cex.axis = 0.4)

# Add legend
legend("topright", legend = c("Not Alive", "Age 0–14", "Age >14"),
       fill = plot_colors, cex = 0.8, border = NA)

### clipped:
# Recode age_class_clipped for plotting
plot_matrix <- age_class_clipped
plot_matrix[plot_matrix == 1] <- 1  # class 1 (0–14 y)
plot_matrix[plot_matrix == 2] <- 2  # class 2 (>14 y)
plot_matrix[is.na(plot_matrix)] <- 0  # non-alive = white

# Flip
plot_matrix <- plot_matrix[nrow(plot_matrix):1, ]

# Define colors: white = NA, blue = class 1, orange = class 2
plot_colors <- c("white", "skyblue", "orange")

# Plot
image(
  x = 1:ncol(plot_matrix),
  y = 1:nrow(plot_matrix),
  z = t(plot_matrix),
  col = plot_colors,
  axes = FALSE,
  xlab = "Year",
  ylab = "Individual",
  main = "Age Class Matrix"
)

# Axis labels
axis(1, at = 1:ncol(plot_matrix), labels = colnames(plot_matrix), las = 2, cex.axis = 0.7)
axis(2, at = 1:nrow(plot_matrix), labels = rownames(plot_matrix), las = 1, cex.axis = 0.4)

# Add legend
legend("topright", legend = c("Not Alive", "Age 0–14", "Age >14"),
       fill = plot_colors, cex = 0.8, border = NA)

############################################################
### -------------------- DATA EXPORT ------------------  ###
############################################################

# drop "_clipped" syntax
is_class1 <- is_class1_clipped
is_class2 <- is_class2_clipped
y <- y_clipped
z <- z_clipped
first_seen <- first_seen_clipped

### keep the following for the model, but remove everything else:
# y_clipped
# first_seen_clipped
# z_clipped
# is_class1_clipped
# is_class2_clipped

rm(list = setdiff(ls(), c(
  "y",
  "first_seen",
  "z",
  "is_class1", 
  "is_class2" 
)))

stop('[1.1.mgmt_survivalMatrices.R] \n
All required matrices constructed. Code stopped to prevent overwriting data.\n
Continue running code beyond this line to overwrite data exports.')

save(y,
     first_seen,
     z,
     is_class1,
     is_class2, 
     file = "data/intermediate/adultSurvival_cjsMatrices.rData")
