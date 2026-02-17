# ==========================================================
# SNOW PROXIES USING PRISM (1999–PRESENT)
# 1. April 1 cumulative precip
# 2. Winter precipitation
# 3. Freezing days (Tmean ≤ 0°C)
# ==========================================================

library(prism)
library(terra)
library(sf)
library(lubridate)
library(dplyr)

# ----------------------------------------------------------
# 1. Set PRISM Download Directory
# ----------------------------------------------------------
prism_set_dl_dir("data/master/prism_data")  # Creates folder in working directory

# ----------------------------------------------------------
# 2. Load Northern Range Polygon
# ----------------------------------------------------------
nr <- vect("data/master/gis/nr_bound/Northern_Range_Bound.shp")   

# ----------------------------------------------------------
# 3. Date Range (Start 1 year early for WY 1999)
# ----------------------------------------------------------
start_year <- 1994
end_year   <- 2023

# # ----------------------------------------------------------
# # 4. Download Monthly PRISM Data 
# 
# *** Note: this is time-intensive and may have already been done. ***
#
# # ----------------------------------------------------------
# Years to download
years <- start_year:end_year

# Months to download
months <- 1:12

# Precipitation
get_prism_monthlys(
  type = "ppt",
  years = years,
  mon = months,
  keepZip = FALSE
)

# Mean temperature
get_prism_monthlys(
  type = "tmean",
  years = years,
  mon = months,
  keepZip = FALSE
)

# ----------------------------------------------------------
# 5. Load Rasters
# ----------------------------------------------------------
# List archive folders
ppt_archives <- prism_archive_ls() |> 
  grep("ppt", x = _, value = TRUE)

tmean_archives <- prism_archive_ls() |> 
  grep("tmean", x = _, value = TRUE)

# Build raster stacks properly
ppt_stack_r <- prism_stack(ppt_archives)
tmean_stack_r <- prism_stack(tmean_archives)

# Convert to terra SpatRaster
ppt_stack <- rast(ppt_stack_r)
tmean_stack <- rast(tmean_stack_r) / 100

# ----------------------------------------------------------
# 6. Crop to Northern Range
# ----------------------------------------------------------

nr <- project(nr, crs(ppt_stack))

ppt_nr <- crop(ppt_stack, nr) |> mask(nr)
tmean_nr <- crop(tmean_stack, nr) |> mask(nr)

# ----------------------------------------------------------
# 7. Extract Monthly Area Means
# ----------------------------------------------------------
ppt_monthly <- global(ppt_nr, mean, na.rm = TRUE)[,1]
tmean_monthly <- global(tmean_nr, mean, na.rm = TRUE)[,1]

# Build date index
dates <- seq(as.Date(sprintf("%d-01-01", start_year)),
             as.Date(sprintf("%d-12-01", end_year)),
             by = "month")

water_year <- ifelse(month(dates) >= 10,
                     year(dates) + 1,
                     year(dates))

# ----------------------------------------------------------
# 8. Winter Precipitation (Dec–Mar)
# ----------------------------------------------------------
winter_idx <- month(dates) %in% c(12,1,2,3)

winter_precip <- tapply(
  ppt_monthly[winter_idx],
  water_year[winter_idx],
  sum,
  na.rm = TRUE
)

# ----------------------------------------------------------
# 9. April 1 Proxy (Oct–Mar cumulative precip)
# ----------------------------------------------------------
oct_mar_idx <- month(dates) %in% c(10,11,12,1,2,3)

april1_precip <- tapply(
  ppt_monthly[oct_mar_idx],
  water_year[oct_mar_idx],
  sum,
  na.rm = TRUE
)

# ----------------------------------------------------------
# 10. Freezing Months (mean temp ≤ 0°C)
# ----------------------------------------------------------
freezing_months <- tapply(
  tmean_monthly <= 0,
  water_year,
  sum,
  na.rm = TRUE
)

# ----------------------------------------------------------
# 11. Combine
# ----------------------------------------------------------
snow_metrics <- data.frame(
  water_year = as.numeric(names(winter_precip)),
  winter_precip = as.numeric(winter_precip),
  april1_winter_precip = as.numeric(april1_precip[names(winter_precip)]),
  freezing_months = as.numeric(freezing_months[names(winter_precip)])
) |> arrange(water_year)

print(snow_metrics)

# ----------------------------------------------------------
# 12. Clean up workspace and export data
# ----------------------------------------------------------

rm(nr, ppt_nr, ppt_stack, ppt_stack_r, tmean_nr, tmean_stack, tmean_stack_r, april1_precip,
   dates, end_year, freezing_months, oct_mar_idx, ppt_archives, ppt_monthly, start_year,
   tmean_archives, tmean_monthly, water_year, winter_idx, winter_precip)

stop('The following will overwrite data. Are you sure you would like to proceed?')

write.csv(snow_metrics,
          "data/intermediate/prism_monthly_snow.csv",
          row.names = FALSE)
