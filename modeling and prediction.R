################################################################################################################################################
##################### This code is prepared by Alemneh Mekuriaw Liyew and  Kefyalew Addiss Alene
##########################################################################################
############ 1. Exploratory data analysis and covariate data preparation ############# #
##########################################################################################
####install and call the library for relevant packages 
library(snakecase)
library(janitor)
library(readxl)
library(readxl)
library(caret)
library(gam)
library(glmnet)
library(gbm)
library(mgcv)
library(kergp)

install.packages("rnaturalearth")
install.packages("rnaturalearthdata")
install.packages("ggplot2")
# Load libraries
library(rnaturalearth)
library(rnaturalearthdata)
library(ggplot2)
library(sf)
# Extract African country boundaries
africa <- ne_countries(continent = "Africa", returnclass = "sf")
eda_sf <- st_as_sf(eda, coords = c("longitude", "latitude"), crs = 4326)
eda_by_year <- eda_sf %>%
  mutate(year = year) %>%  
  group_by(year) %>%
  arrange(year)  # Arrange data by year

# Define year intervals
year_breaks <- c(1990, 2000, 2010,2020, 2023)  # Define boundaries for intervals 

labels = c("1990-1999", "2000-2009", "2010-2019", "2020-2022")

eda_by_year <- eda_by_year %>%
  mutate(year_interval = cut(year, breaks = year_breaks,
                             labels = c("1990-1999", "2000-2009", "2010-2019", "2020-2022"), 
                             right = FALSE)  # Ensure intervals are left-close

p <- ggplot() +
  geom_sf(data = africa) +  # Plot African countries as a base layer
  geom_sf(data = eda_by_year, aes(color = prop_vac), size = 1.5, alpha = 0.5) +  # Map 'cases' to color aesthetics
  scale_color_viridis_c(name = "BCG coverage per Cluster") +  # Choose a color scale suitable for proportions
  facet_wrap(~ year_interval, ncol = 2) +  # Facet by year with 2 columns per row
  theme(strip.text = element_text(size = 12),
        axis.text = element_text(size = 10),
        plot.background = element_rect(fill = "lightblue"))  # Change background color here
p <- ggplot() +
  geom_sf(data = africa) +
  geom_sf(data = eda_by_year, aes(color = prop_vac), size = 1.5, alpha = 0.5) +
  scale_color_viridis_c(name = "BCG coverage per Cluster") +
  facet_wrap(~ year_interval, ncol = 2) +
  coord_sf(expand = FALSE) +   # Important!
  theme(
    strip.text = element_text(size = 12),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.title = element_blank(),
    panel.grid = element_blank(),
    plot.background = element_rect(fill = "lightblue")
  )

####################################################################################################
###############2. extracting and combining covariates to outcome ######################################### #
######################################################################################
# Install the packages you need
install.packages("INLA",
                 repos = "https://inla.r-inla-download.org/R/stable", dep = TRUE)
## Loading required package: Matrix
library(Matrix)
library(sp)
library(INLA)
library(readxl)
library(geodata)
library(terra)
library(ggplot2)

# define the paths
path_input <- paste0(getwd(),"/INPUT")

# study area
Africa <- ne_countries(scale = "medium", continent = "Africa", returnclass = "sf")
# Filter for the African countries in your ISO3 list from BCG0
BCG0<-eda
iso3_list <- unique(BCG0$iso_a3)
myarea <- Africa %>% 
  filter(iso_a3 %in% iso3_list)
myarea1<-as. data.frame(myarea)
iso3 <- myarea1[, c(10, 11)]
# Save as CSV
write.csv(iso3, "C:/Users/ALiyew/Desktop/data/bcg_12-23m data/iso3.csv", row.names = FALSE)

# Plot the geometric boundaries of the selected study area using ggplot2
ggplot(data = myarea) +
  geom_sf(fill = "transparent", color = "black") +
  ggtitle("Study Area") +
  theme_minimal()
study_area <- ggplot() +
  geom_sf(data = myarea, fill = "red", color = NA) +  # Set fill color to red for myarea
  geom_sf(data = africa, fill = NA, color = "black") +  # Add another spatial layer
  theme(
    panel.grid.major = element_blank(),  # Remove major gridlines
    panel.grid.minor = element_blank(),  # Remove minor gridlines
    axis.text = element_blank(),         # Remove axis text
    axis.title = element_blank(),        # Remove axis titles
    axis.ticks = element_blank()         # Remove axis ticks
  )
plot(study_area)
####geolocated outcome data

BCG0$total_c <- as.integer(BCG0$total_c)
BCG0$prop_vac <- as.numeric(BCG0$prop_vac)
BCG0$'vaccinated_c' <- as.numeric(BCG0$'vaccinated_c')
BCG0$latitude <- as.numeric(BCG0$latitude )
BCG0$longitude <- as.numeric(BCG0$longitude)
BCG <-BCG0
BCG_sorted <- BCG[order(BCG$year), ]
unique(BCG_sorted$year)
#**********combine the coordinate**********************
xy <- cbind( BCG_sorted$longitude,  BCG_sorted$latitude)
###########################################################################################
###2.1. processing temporal covariates
##############################################################################################
library(raster)
###PM2.5 1998-2021################################
path_input <- paste0(getwd(), "/all_cov")
raster_files <- list.files(path_input, pattern = "HybridPM25.*\\.tif$",full.names = TRUE)
# Extract year from filenames like: V5GL03.HybridPM25.Global.199801-199812-geotiff.tif
extract_year <- function(fname) {
  # Extract the first 4-digit number in the filename (which corresponds to the start year)
  match <- regmatches(fname, regexpr("\\d{4}", fname))
  as.numeric(match)
}
basename(raster_files)
raster_years <- sapply(basename(raster_files), extract_year)
print(raster_years)
# Unique years to process
years_to_process <- sort(unique(raster_years))
print(years_to_process)


###for memory split into four periods
years_q1 <- years_to_process[1:ceiling(length(years_to_process)/4)]

pm_processed1 <- list()
pm_resampled1 <- list()

for (yr in years_q1) {
  message("Processing year (Q1): ", yr)
  idx <- which(raster_years == yr)
  if (length(idx) == 1) {
    r_path <- raster_files[idx]
    pm_raw <- raster(r_path)
    pm_cropped <- crop(pm_raw, extent(myarea))
    pm_masked  <- mask(pm_cropped, myarea)
    pm_processed1[[as.character(yr)]] <- pm_masked
    
    # Resample
    pm_resampled1[[as.character(yr)]] <- resample(pm_masked, alt, method = "bilinear")
  } else {
    warning("No raster file found for year: ", yr)
  }
}
remove(pm_processed1)###to save memory
###next quarter of year
n <- length(years_to_process)
q2_start <- ceiling(n/4) + 1
q2_end   <- ceiling(n/2)
years_q2 <- years_to_process[q2_start:q2_end]

pm_processed2 <- list()
pm_resampled2 <- list()

for (yr in years_q2) {
  message("Processing year (Q2): ", yr)
  idx <- which(raster_years == yr)
  if (length(idx) == 1) {
    r_path <- raster_files[idx]
    pm_raw <- raster(r_path)
    pm_cropped <- crop(pm_raw, extent(myarea))
    pm_masked  <- mask(pm_cropped, myarea)
    pm_processed2[[as.character(yr)]] <- pm_masked
    
    pm_resampled2[[as.character(yr)]] <- resample(pm_masked, alt, method = "bilinear")
  } else {
    warning("No raster file found for year: ", yr)
  }
}
remove(pm_processed2)
pm_resampled_combined1 <- c(pm_resampled1, pm_resampled2)
# Create list to store normalized rasters
pm_normalized1 <- list()
# Loop through each year in the combined resampled list
for (yr in names(pm_resampled_combined1)) {
  message("Normalizing year: ", yr)
  r <- pm_resampled_combined1[[yr]]
  vals <- getValues(r)
  # Remove NAs and compute quantile-normalized values
  n_vals <- length(na.omit(vals))
  norm_vals <- qnorm(seq(0.0, 1, length.out = n_vals + 2)[2:(n_vals + 1)])
  # Assign normalized ranks
  norm_ranked <- norm_vals[rank(na.omit(vals))]
  # Replace values in raster
  full_vals <- rep(NA, length(vals))
  full_vals[!is.na(vals)] <- norm_ranked
  r_norm <- r  # copy raster structure
  values(r_norm) <- full_vals
  pm_normalized1[[yr]] <- r_norm
}
remove(pm_resampled1)
remove(pm_resampled2)
remove(pm_resampled_combined1)

#####third quarter
q3_start <- ceiling(n/2) + 1
q3_end   <- ceiling(3*n/4)
years_q3 <- years_to_process[q3_start:q3_end]

pm_processed3 <- list()
pm_resampled3 <- list()

for (yr in years_q3) {
  message("Processing year (Q3): ", yr)
  idx <- which(raster_years == yr)
  if (length(idx) == 1) {
    r_path <- raster_files[idx]
    pm_raw <- raster(r_path)
    pm_cropped <- crop(pm_raw, extent(myarea))
    pm_masked  <- mask(pm_cropped, myarea)
    pm_processed3[[as.character(yr)]] <- pm_masked
    pm_resampled3[[as.character(yr)]] <- resample(pm_masked, alt, method = "bilinear")
  } else {
    warning("No raster file found for year: ", yr)
  }
}
remove(pm_processed3)
####fourth quarter
q4_start <- ceiling(3*n/4) + 1
years_q4 <- years_to_process[q4_start:n]
pm_processed4 <- list()
pm_resampled4 <- list()
for (yr in years_q4) {
  message("Processing year (Q4): ", yr)
  idx <- which(raster_years == yr)
  if (length(idx) == 1) {
    r_path <- raster_files[idx]
    pm_raw <- raster(r_path)
    pm_cropped <- crop(pm_raw, extent(myarea))
    pm_masked  <- mask(pm_cropped, myarea)
    pm_processed4[[as.character(yr)]] <- pm_masked
    pm_resampled4[[as.character(yr)]] <- resample(pm_masked, alt, method = "bilinear")
  } else {
    warning("No raster file found for year: ", yr)
  }
}
pm_resampled_combined2 <- c(pm_resampled3, pm_resampled4)
pm_normalized2 <- list()
# Loop through each year in the second combined resampled list
for (yr in names(pm_resampled_combined2)) {
  message("Normalizing year: ", yr)
  r <- pm_resampled_combined2[[yr]]
  vals <- getValues(r)
  # Remove NAs and compute quantile-normalized values
  n_vals <- length(na.omit(vals))
  norm_vals <- qnorm(seq(0.0, 1, length.out = n_vals + 2)[2:(n_vals + 1)])
  # Assign normalized ranks
  norm_ranked <- norm_vals[rank(na.omit(vals))]
  # Replace values in raster
  full_vals <- rep(NA, length(vals))
  full_vals[!is.na(vals)] <- norm_ranked
  r_norm <- r  # copy raster structure
  values(r_norm) <- full_vals
  pm_normalized2[[yr]] <- r_norm
}
remove(pm_processed4)
pm_all<-c(pm_normalized1,pm_normalized2)


unique(BCG_sorted$year)[!unique(BCG_sorted$year) %in% names(pm_all)]
# Prepare an empty vector to store extracted values
pm_values <- numeric(nrow(xy))
# Loop over each row in BCG_sorted

for (i in seq_len(nrow(BCG_sorted))) {
  year_i <- as.character(BCG_sorted$year[i])
  year_num <- as.numeric(year_i)
  
  if (year_i %in% names(pm_all)) {
    raster_i <- pm_all[[year_i]]
  } else if (year_num < 1998 && "1998" %in% names(pm_all)) {
    raster_i <- pm_all[["1998"]]
    warning("Using 1998 raster as fallback for year: ", year_i, " at row ", i)
  } else if (year_num > 2021 && "2021" %in% names(pm_all)) {
    raster_i <- pm_all[["2021"]]
    warning("Using 2021 raster as fallback for year: ", year_i, " at row ", i)
  } else {
    pm_values[i] <- NA
    warning("No raster available for year: ", year_i, " at row ", i)
    next
  }
  
  point_i <- xy[i, , drop = FALSE]
  pm_values[i] <- extract(raster_i, point_i)
}
# Add the extracted PM values to your data frame
BCG_sorted$pm25 <- pm_values
write.csv(BCG_sorted, file = file.path(path_input, "BCG_sorted_PM25.csv"), row.names = FALSE)
all_pm_objects <- ls(pattern = "^pm_")

# Specify the ones you want to keep
keep <- c("pm_resampled3", "pm_resampled4", "pm_all","pm_values")
# Remove all except the ones to keep
rm(list = setdiff(all_pm_objects, keep))
#############################################################################
###maternal education

path_input <- paste0(getwd(), "/all_cov")

# List all maternal education files
edu_files <- list.files(path_input, pattern = "IHME_AFRICA_EDU.*\\.TIF$", full.names = TRUE)

extract_year <- function(fname) {
  match <- regmatches(fname, regexpr("MEAN_\\d{4}", fname))
  as.numeric(sub("MEAN_", "", match))
}
edu_years <- sapply(basename(edu_files), extract_year)
edu_rasters <- list()
for (i in seq_along(edu_files)) {
  yr <- as.character(edu_years[i])
  r <- raster::raster(edu_files[i])
  edu_rasters[[yr]] <- r
}
print(names(edu_rasters))

#####masking and cropping to study area

years_to_process <- sort(edu_years)

# Split into halves — first half only
n <- length(years_to_process)
half_point <- ceiling(n / 2)
years_first_half <- years_to_process[1:half_point]
edu_processed1 <- list()
edu_resampled1 <- list()
# Loop over the first half years
for (yr in years_first_half) {
  message("Processing EDU year (first half): ", yr)
  idx <- which(edu_years == yr)
  if (length(idx) == 1) {
    r_path <- edu_files[idx]
    edu_raw <- raster(r_path)
    edu_cropped <- crop(edu_raw, extent(myarea))
    edu_masked <- mask(edu_cropped, myarea)
    # Save processed raster
    edu_processed1[[as.character(yr)]] <- edu_masked
    edu_resampled1[[as.character(yr)]] <- resample(edu_masked, alt, method = "bilinear")
  } else {
    warning("No education raster found for year: ", yr)
  }
}
remove(edu_processed1)

# Second half year
years_second_half <- years_to_process[(half_point + 1):n]
edu_processed2 <- list()
edu_resampled2 <- list()
# Loop over the second half year
for (yr in years_second_half) {
  message("Processing EDU year (second half): ", yr)
  idx <- which(edu_years == yr)
  if (length(idx) == 1) {
    r_path <- edu_files[idx]
    edu_raw <- raster(r_path)
    edu_cropped <- crop(edu_raw, extent(myarea))
    edu_masked <- mask(edu_cropped, myarea)
    edu_processed2[[as.character(yr)]] <- edu_masked
    edu_resampled2[[as.character(yr)]] <- resample(edu_masked, alt, method = "bilinear")
  } else {
    warning("No education raster found for year: ", yr)
  }
}
remove(edu_processed2)
edu_resampled_all<-c(edu_resampled1,edu_resampled2)
names(edu_resampled_all)

edu_normalized_all <- list()

for (yr in names(edu_resampled_all)) {
  message("Normalizing EDU year: ", yr)
  r <- edu_resampled_all[[yr]]
  vals <- getValues(r)
  n_vals <- length(na.omit(vals))
  norm_vals <- qnorm(seq(0.0, 1, length.out = n_vals + 2)[2:(n_vals + 1)])
  norm_ranked <- norm_vals[rank(na.omit(vals))]
  full_vals <- rep(NA, length(vals))
  full_vals[!is.na(vals)] <- norm_ranked
  r_norm <- r  # copy raster structure
  values(r_norm) <- full_vals
  edu_normalized_all[[yr]] <- r_norm
}
edu_values <- numeric(nrow(xy))

for (i in seq_len(nrow(BCG_sorted))) {
  year_i <- as.character(BCG_sorted$year[i])
  year_num <- as.numeric(year_i)
  
  if (year_i %in% names(edu_normalized_all)) {
    raster_i <- edu_normalized_all[[year_i]]
  } else if (year_num < 2000 && "2000" %in% names(edu_normalized_all)) {
    raster_i <- edu_normalized_all[["2000"]]
    warning("Using 2000 raster as fallback for year: ", year_i, " at row ", i)
  } else if (year_num > 2015 && "2015" %in% names(edu_normalized_all)) {
    raster_i <- edu_normalized_all[["2015"]]
    warning("Using 2015 raster as fallback for year: ", year_i, " at row ", i)
  } else {
    edu_values[i] <- NA
    warning("No raster available for year: ", year_i, " at row ", i)
    next
  }
  
  point_i <- xy[i, , drop = FALSE]
  edu_values[i] <- extract(raster_i, point_i)
}
BCG_sorted$edu <- edu_values
summary(BCG_sorted)
all_edu_objects <- ls(pattern = "^edu_")


# Specify the ones you want to keep
keep <- c("edu_normalized_all","edu_values")
# Remove all except the ones to keep
rm(list = setdiff(all_edu_objects, keep))
#############################################################################
###exclusive breastfeeding 2000-2019

path_input <- paste0(getwd(), "/all_cov")

# List all maternal education files
ebf_files <- list.files(path_input, pattern = "GLOBAL_EBF.*\\.TIF$", full.names = TRUE)
extract_year <- function(fname) {
  match <- regmatches(fname, regexpr("MEAN_\\d{4}", fname))
  as.numeric(sub("MEAN_", "", match))
}
ebf_years <- sapply(basename(ebf_files), extract_year)

ebf_rasters <- list()
for (i in seq_along(ebf_files)) {
  yr <- as.character(ebf_years[i])
  r <- raster::raster(ebf_files[i])
  ebf_rasters[[yr]] <- r
}
print(names(ebf_rasters))
years_to_process <- sort(ebf_years)
# Split into halves — first half only
n <- length(years_to_process)
half_point <- ceiling(n / 2)
years_first_half <- years_to_process[1:half_point]
ebf_processed1 <- list()
ebf_resampled1 <- list()
# Loop over the first half years
for (yr in years_first_half) {
  message("Processing EBF year (first half): ", yr)
  idx <- which(ebf_years == yr)
  if (length(idx) == 1) {
    r_path <- ebf_files[idx]
    ebf_raw <- raster(r_path)
    ebf_cropped <- crop(ebf_raw, extent(myarea))
    ebf_masked <- mask(ebf_cropped, myarea)
    # Save processed raster
    ebf_processed1[[as.character(yr)]] <- ebf_masked
    ebf_resampled1[[as.character(yr)]] <- resample(ebf_masked, alt, method = "bilinear")
  } else {
    warning("No EBF raster found for year: ", yr)
  }
}
remove(ebf_processed1)
# Identify remaining (unprocessed) years
years_remaining <- setdiff(years_to_process, years_first_half)
ebf_processed_remaining <- list()
ebf_resampled_remaining <- list()
# Loop over remaining years
for (yr in years_remaining) {
  message("Processing EBF year (remaining): ", yr)
  idx <- which(ebf_years == yr)
  if (length(idx) == 1) {
    r_path <- ebf_files[idx]
    ebf_raw <- raster(r_path)
    ebf_cropped <- crop(ebf_raw, extent(myarea))
    ebf_masked <- mask(ebf_cropped, myarea)
    # Save processed raster
    ebf_processed_remaining[[as.character(yr)]] <- ebf_masked
    ebf_resampled_remaining[[as.character(yr)]] <- resample(ebf_masked, alt, method = "bilinear")
  } else {
    warning("No EBF raster found for year: ", yr)
  }
}
remove(ebf_processed_remaining)
ebf_resampled_all <- c(ebf_resampled1, ebf_resampled_remaining)
names(ebf_resampled_all)

ebf_normalized_all <- list()

for (yr in names(ebf_resampled_all)) {
  message("Normalizing EBF year: ", yr)
  r <- ebf_resampled_all[[yr]]
  vals <- getValues(r)
  n_vals <- length(na.omit(vals))
  norm_vals <- qnorm(seq(0.0, 1, length.out = n_vals + 2)[2:(n_vals + 1)])
  norm_ranked <- norm_vals[rank(na.omit(vals))]
  full_vals <- rep(NA, length(vals))
  full_vals[!is.na(vals)] <- norm_ranked
  r_norm <- r  # copy raster structure
  values(r_norm) <- full_vals
  ebf_normalized_all[[yr]] <- r_norm
}

ebf_values <- numeric(nrow(xy))

for (i in seq_len(nrow(BCG_sorted))) {
  year_i <- as.character(BCG_sorted$year[i])
  year_num <- as.numeric(year_i)
  
  if (year_i %in% names(ebf_normalized_all)) {
    raster_i <- ebf_normalized_all[[year_i]]
  } else if (year_num < 2000 && "2000" %in% names(ebf_normalized_all)) {
    raster_i <- ebf_normalized_all[["2000"]]
    warning("Using 2000 raster as fallback for year: ", year_i, " at row ", i)
  } else if (year_num > 2019 && "2019" %in% names(ebf_normalized_all)) {
    raster_i <- ebf_normalized_all[["2019"]]
    warning("Using 2019 raster as fallback for year: ", year_i, " at row ", i)
  } else {
    ebf_values[i] <- NA
    warning("No raster available for year: ", year_i, " at row ", i)
    next
  }
  
  point_i <- xy[i, , drop = FALSE]
  ebf_values[i] <- extract(raster_i, point_i)
}
BCG_sorted$ebf <- ebf_values
summary(BCG_sorted)

all_ebf_objects <- ls(pattern = "^ebf_")
keep <- c("ebf_normalized_all", "ebf_values")
rm(list = setdiff(all_ebf_objects, keep))
######################################################################
###for unimproved sanitation

# List all USAN files
usan_files <- list.files(path_input, pattern = "_S_UNIMP.*\\.TIF$", full.names = TRUE)

extract_year <- function(fname) {
  match <- regmatches(fname, regexpr("MEAN_\\d{4}", fname))
  as.numeric(sub("MEAN_", "", match))
}

usan_years <- sapply(basename(usan_files), extract_year)

usan_rasters <- list()
for (i in seq_along(usan_files)) {
  yr <- as.character(usan_years[i])
  r <- raster::raster(usan_files[i])
  usan_rasters[[yr]] <- r
}

print(names(usan_rasters))

years_to_process <- sort(usan_years)
### split for memory purposes
n <- length(years_to_process)
half_point <- ceiling(n / 2)
years_first_half <- years_to_process[1:half_point]

usan_processed1 <- list()
usan_resampled1 <- list()

# Loop over the first half years
for (yr in years_first_half) {
  message("Processing USAN year (first half): ", yr)
  idx <- which(usan_years == yr)
  if (length(idx) == 1) {
    r_path <- usan_files[idx]
    usan_raw <- raster(r_path)
    usan_cropped <- crop(usan_raw, extent(myarea))
    usan_masked <- mask(usan_cropped, myarea)
    # Save processed raster
    usan_processed1[[as.character(yr)]] <- usan_masked
    usan_resampled1[[as.character(yr)]] <- resample(usan_masked, alt, method = "bilinear")
  } else {
    warning("No USAN raster found for year: ", yr)
  }
}

remove(usan_processed1)

# Identify remaining (unprocessed) years
years_remaining <- setdiff(years_to_process, years_first_half)

usan_processed_remaining <- list()
usan_resampled_remaining <- list()

# Loop over remaining years
for (yr in years_remaining) {
  message("Processing USAN year (remaining): ", yr)
  idx <- which(usan_years == yr)
  if (length(idx) == 1) {
    r_path <- usan_files[idx]
    usan_raw <- raster(r_path)
    usan_cropped <- crop(usan_raw, extent(myarea))
    usan_masked <- mask(usan_cropped, myarea)
    # Save processed raster
    usan_processed_remaining[[as.character(yr)]] <- usan_masked
    usan_resampled_remaining[[as.character(yr)]] <- resample(usan_masked, alt, method = "bilinear")
  } else {
    warning("No USAN raster found for year: ", yr)
  }
}

usan_resampled_all <- c(usan_resampled1, usan_resampled_remaining)
names(usan_resampled_all)

usan_normalized_all <- list()

for (yr in names(usan_resampled_all)) {
  message("Normalizing USAN year: ", yr)
  r <- usan_resampled_all[[yr]]
  vals <- getValues(r)
  n_vals <- length(na.omit(vals))
  norm_vals <- qnorm(seq(0.0, 1, length.out = n_vals + 2)[2:(n_vals + 1)])
  norm_ranked <- norm_vals[rank(na.omit(vals))]
  full_vals <- rep(NA, length(vals))
  full_vals[!is.na(vals)] <- norm_ranked
  r_norm <- r  # copy raster structure
  values(r_norm) <- full_vals
  usan_normalized_all[[yr]] <- r_norm
}

usan_values <- numeric(nrow(xy))

for (i in seq_len(nrow(BCG_sorted))) {
  year_i <- as.character(BCG_sorted$year[i])
  year_num <- as.numeric(year_i)
  
  if (year_i %in% names(usan_normalized_all)) {
    raster_i <- usan_normalized_all[[year_i]]
  } else if (year_num < 2000 && "2000" %in% names(usan_normalized_all)) {
    raster_i <- usan_normalized_all[["2000"]]
    warning("Using 2000 raster as fallback for year: ", year_i, " at row ", i)
  } else if (year_num > 2017 && "2017" %in% names(usan_normalized_all)) {
    raster_i <- usan_normalized_all[["2017"]]
    warning("Using 2017 raster as fallback for year: ", year_i, " at row ", i)
  } else {
    usan_values[i] <- NA
    warning("No raster available for year: ", year_i, " at row ", i)
    next
  }
  
  point_i <- xy[i, , drop = FALSE]
  usan_values[i] <- extract(raster_i, point_i)
}
BCG_sorted$usan <- usan_values
summary(BCG_sorted)

all_usan_objects <- ls(pattern = "^usan_")
keep <- c("usan_normalized_all", "usan_values")
rm(list = setdiff(all_usan_objects, keep))
#####################################################################
###unimproved water supply
# List all UWAT files
uwat_files <- list.files(path_input, pattern = "_W_UNIMP.*\\.TIF$", full.names = TRUE)
extract_year <- function(fname) {
  match <- regmatches(fname, regexpr("MEAN_\\d{4}", fname))
  as.numeric(sub("MEAN_", "", match))
}
uwat_years <- sapply(basename(uwat_files), extract_year)
uwat_rasters <- list()
for (i in seq_along(uwat_files)) {
  yr <- as.character(uwat_years[i])
  r <- raster::raster(uwat_files[i])
  uwat_rasters[[yr]] <- r
}
print(names(uwat_rasters))
years_to_process <- sort(uwat_years)
### Split for memory purposes
n <- length(years_to_process)
half_point <- ceiling(n / 2)
years_first_half <- years_to_process[1:half_point]

uwat_processed1 <- list()
uwat_resampled1 <- list()

# Loop over the first half years
for (yr in years_first_half) {
  message("Processing UWAT year (first half): ", yr)
  idx <- which(uwat_years == yr)
  if (length(idx) == 1) {
    r_path <- uwat_files[idx]
    uwat_raw <- raster(r_path)
    uwat_cropped <- crop(uwat_raw, extent(myarea))
    uwat_masked <- mask(uwat_cropped, myarea)
    # Save processed raster
    uwat_processed1[[as.character(yr)]] <- uwat_masked
    uwat_resampled1[[as.character(yr)]] <- resample(uwat_masked, alt, method = "bilinear")
  } else {
    warning("No UWAT raster found for year: ", yr)
  }
}
remove(uwat_processed1)
# Identify remaining (unprocessed) years
years_remaining <- setdiff(years_to_process, years_first_half)
uwat_processed_remaining <- list()
uwat_resampled_remaining <- list()

# Loop over remaining years
for (yr in years_remaining) {
  message("Processing UWAT year (remaining): ", yr)
  idx <- which(uwat_years == yr)
  if (length(idx) == 1) {
    r_path <- uwat_files[idx]
    uwat_raw <- raster(r_path)
    uwat_cropped <- crop(uwat_raw, extent(myarea))
    uwat_masked <- mask(uwat_cropped, myarea)
    # Save processed raster
    uwat_processed_remaining[[as.character(yr)]] <- uwat_masked
    uwat_resampled_remaining[[as.character(yr)]] <- resample(uwat_masked, alt, method = "bilinear")
  } else {
    warning("No UWAT raster found for year: ", yr)
  }
}

uwat_resampled_all <- c(uwat_resampled1, uwat_resampled_remaining)
names(uwat_resampled_all)


uwat_normalized_all <- list()

for (yr in names(uwat_resampled_all)) {
  message("Normalizing UWAT year: ", yr)
  r <- uwat_resampled_all[[yr]]
  vals <- getValues(r)
  n_vals <- length(na.omit(vals))
  norm_vals <- qnorm(seq(0.0, 1, length.out = n_vals + 2)[2:(n_vals + 1)])
  norm_ranked <- norm_vals[rank(na.omit(vals))]
  full_vals <- rep(NA, length(vals))
  full_vals[!is.na(vals)] <- norm_ranked
  r_norm <- r  # copy raster structure
  values(r_norm) <- full_vals
  uwat_normalized_all[[yr]] <- r_norm
}

uwat_values <- numeric(nrow(xy))

for (i in seq_len(nrow(BCG_sorted))) {
  year_i <- as.character(BCG_sorted$year[i])
  year_num <- as.numeric(year_i)
  
  if (year_i %in% names(uwat_normalized_all)) {
    raster_i <- uwat_normalized_all[[year_i]]
  } else if (year_num < 2000 && "2000" %in% names(uwat_normalized_all)) {
    raster_i <- uwat_normalized_all[["2000"]]
    warning("Using 2000 raster as fallback for year: ", year_i, " at row ", i)
  } else if (year_num > 2017 && "2017" %in% names(uwat_normalized_all)) {
    raster_i <- uwat_normalized_all[["2017"]]
    warning("Using 2017 raster as fallback for year: ", year_i, " at row ", i)
  } else {
    uwat_values[i] <- NA
    warning("No raster available for year: ", year_i, " at row ", i)
    next
  }
  
  point_i <- xy[i, , drop = FALSE]
  uwat_values[i] <- extract(raster_i, point_i)
}
BCG_sorted$uwat <- uwat_values
summary(BCG_sorted)

all_uwat_objects <- ls(pattern = "^uwat_")
keep <- c("uwat_normalized_all", "uwat_values")
rm(list = setdiff(all_uwat_objects, keep))
#*************Temp from worldclim*******************************************************************************************
library(raster)
path_input <- paste0(getwd(), "/all_cov")
TO <- list.files(path_input, pattern = "wc2\\.1_30s_tavg_\\d{2}\\.tif$", full.names = TRUE) # List the 12 monthly tavg rasters over 1983-2000
T_stack <- stack(TO)
T0 <- calc(T_stack, fun = mean, na.rm = TRUE)##average annual temperature

T0 <- T0[[c(1)]]
#*************Prec from worldclim*****************************************************************************************
PO <- list.files(path_input, pattern = "wc2\\.1_30s_prec_\\d{2}\\.tif$", full.names = TRUE) # List the 12 monthly prece rasters over 1983-2000
P_stack <- stack(TO)
P0 <- calc(T_stack, fun = mean, na.rm = TRUE)##average annual temperature

#########################################################################################
###2.2. Static covariates
#######################################################################################

#**************altitude*******************************************************************************************************
path_input <- paste0(getwd(),"/INPUT")
alt0<-geodata::elevation_global(res=5,path=path_input)
alt0<-raster::raster(alt0)
#*************travel time to health facility from MAP*********************************************************************************************
eh0 <-raster::raster(paste0(path_input,"/Annual Global High-Resolution Extreme Heat Estimates (GEHE), 1983-2016 wbgtmax32-trend-p05.tif"))
ahf0 <- raster(paste0(path_input,"/2020_walking_only_travel_time_to_healthcare.tif"))
acc0 <- raster::raster(paste0(path_input,"/accAfrica.tif"))
#**********************sociodemographic factors**********************************************************************************************
popden0 <-raster(paste0(path_input,"/gpw-v4-population-density_2000.tif"))
rdi0 <-raster(paste0(path_input,"/Relative Deprivation Index 2010-2020-v1.tif"))
###*** Nitrogen dioxide (NO2)***
N0 <-raster(paste0(path_input,"/SURFACE_NO2_010x010_2011.tif"))
##########################################################################
path_output <- file.path(getwd(),"OUTPUT")

  plot(T0)
  T1 <- mask(crop(T0, extent(myarea)),myarea)
  plot(T1, main = "Annual Average Temperature")
  ## crop and mask
  P1 <- mask(crop(P0, extent(myarea)),myarea)
  ## Check that it worked
  plot(P1, main="Precipitation")
  alt <- mask(crop(alt0, extent(myarea)),myarea)
  plot(alt, main="Altitude")
  plot(alt, main="Altitude")
  acc <- mask(crop(acc0, extent(myarea)),myarea)
  plot(acc)
  ahf <- mask(crop(ahf0, extent(myarea)),myarea)
  plot(ahf)
  popden <- mask(crop(popden0, extent(myarea)),myarea)
  plot(popden)
  rdi <- mask(crop(rdi0, extent(myarea)),myarea)
  plot(rdi)
  no2 <- mask(crop(N0, extent(myarea)),myarea)
  plot(no2)
  eh <- mask(crop(eh0, extent(myarea)),myarea)
  #*******stack rasters using stack(raster1,raster2,...)************************************************
  # as the variables are from different sources and dimensiton*so it needs to b resampleed before stack************
  T2 <- resample(T1,alt)
  eh2 <- resample(eh1,alt)
  P2 <- resample(P1,alt, method='ngb')
  acc2 <- resample(acc,alt)
  ahf2 <- resample(ahf,alt)
  popden2 <- resample(popden,alt)
  
  rdi2 <- resample(rdi,alt)
  
  no2<- resample(no2,alt)
  rs <- stack(T2,eh2,P2,alt,acc2,ahf2,popden2,rdi2,no2)
  names(rs) <- c("TMP", "EH", "PCP", "ALT", "ACC", "ACH", "POP", "RDI" , "NO2")
  
  #*************x-mean/sd to standardize the covariates rs ****************************************
  rs2 <-scale(rs)
  #***Projected_rs2 <- crs(rs2, asText=TRUE)
  #quantile normalization (rank values and make them correspond to normally distributed data)
  norm<-list()
  new_var<-list()
  new_var_full<-list()
  st2<-list()
  x<-list()
  n<-list()
  for (i in 1:nlayers(rs)){
    # linear is my raster
    st2[[i]] <- rs[[i]]
    x[[i]] <- getValues(rs[[i]])
    n[[i]] <- length(na.omit(x[[i]]))
    norm[[i]] <- qnorm(seq(0.0, 1, length.out = n[[i]] + 2)[2:(n[[i]] + 1)])
    new_var[[i]] <- norm[[i]][base::rank(na.omit(x[[i]]))]
    new_var_full[[i]] <- rep(NA, length(x[[i]]))
    new_var_full[[i]][!is.na(x[[i]])] <- new_var[[i]]
    values(st2[[i]]) <- new_var_full[[i]]
  }
  rs2<-stack(st2)#
  plot(rs2)
  #mapping covariates
  rsmap <- rs
  names(rsmap) <- c("TMP", "EH", "PCP", "ALT" ,"ACC", "ACH", "POP", "RDI", "NO2")
  #log pop for mapping purposes
  rsmap[["POP"]] <- log(subset(rsmap, 'POP')+1)
  names(rsmap[["POP"]]) <- "logPOP"
  af.layer <- function() {plot(Africa,col = "transparent", add=TRUE)}
  plot.new()
  plot(rsmap)
  pdf(file=paste0(path_output,"/pdf/Covariates.pdf"))
 install.packages("viridisLite")
 library(viridisLite)
  library(viridis)
   plot(rsmap,col=viridis(10),nc=4,addfun=af.layer)
  if(!is.null(dev.list())) dev.off() 
  if(!is.null(dev.list())) dev.off()
   #########################################################################
  ### combine temporal and static covariates for mapping (for temporal covariates latest year was used)
  
   # Extract the latest-year raster from a temporal list
   get_latest_raster <- function(rlist) {
     years <- as.numeric(names(rlist))
     latest_year <- max(years, na.rm = TRUE)
     return(rlist[[as.character(latest_year)]])
   }
   # Stack latest-year covariates with rs
   rs_updated <- stack(
     rs2,
     get_latest_raster(ebf_normalized_all),
     get_latest_raster(pm_all),
     get_latest_raster(edu_normalized_all),
     get_latest_raster(usan_normalized_all),
     get_latest_raster(uwat_normalized_all)
   )
   
   # Rename all layers, preserving existing names in rs2
   names(rs_updated) <- c(
     names(rs2),
     "PM2.5", "EDU", "EBF", "USAN", "UWAT"
   )
   
   plot(rs_updated,col=viridis(10),nc=4,addfun=af.layer)
   if(!is.null(dev.list())) dev.off() 
   if(!is.null(dev.list())) dev.off()
   rs2_updated <-scale(rs_updated)

 
####Continue with extracting static covariate at point coordinates. Note that BCG_sorted already contains temporal covariates from the previous step.
   
   xy <- cbind( BCG_sorted$longitude,  BCG_sorted$latitude)
  covariate_all_s <-data.frame(raster::extract(rs2, xy))
  ###combine data frame to have all covariates and outcome in one
  bcgcov<-cbind(BCG_sorted,covariate_all_s)
  covariate_all <- bcgcov[, (ncol(bcgcov)-13):ncol(bcgcov)]
  na_counts <- colSums(is.na(covariate_all))
  print(na_counts)
  ###KNN (k-Nearest Neighbors) imputation is a technique used to fill in missing values in a dataset based on the values of k-nearest (5 in this case)neighboring observations. 
  library(VIM)
  covariate_all_imputed <- kNN(covariate_all, k = 5)
  na_counts <- colSums(is.na(covariate_all_imputed))
  print(na_counts)
  covariate_all<-covariate_all_imputed[, 1:14]
  na_counts <- colSums(is.na(covariate_all))
  print(na_counts)
  BCG_sorted <- BCG[order(BCG$year), ]
  bcgcov<-cbind(BCG_sorted,covariate_all)###data ready for modeling
  
  #********************************************************
  ######################################################################
  ### 3. fitting ensemble submodels
  #######################################################################################

  data<- as. data.frame(bcgcov)
  set.seed(123)  # for reproducibility
  data$group_id <- interaction(data$cluster_id, data$year) # Step 1: Create a group ID for each (cluster, year)
  unique_groups <- unique(data$group_id) # Step 2: Get unique group IDs
  shuffled_groups <- sample(unique_groups)# Step 3: Shuffle and assign groups to folds
  folds <- cut(seq_along(shuffled_groups), breaks = 5, labels = FALSE)
  group_to_fold <- data.frame(group_id = shuffled_groups, fold = folds)  # Step 4: Create a group-to-fold map
  data <- merge(data, group_to_fold, by = "group_id")# Step 5: Merge fold assignments back into original data
  data <- data[order(as.numeric(rownames(data))), ]
  all.covs <- subset(data, select = -c(cluster_id, iso_a3, vaccinated_c,unvaccinated_c, year, total_c,longitude,latitude, prop_vac))
  all.covs<-as.data.frame(all.covs)
  data<-data[complete.cases(all.covs),]
  all.covs<-all.covs[complete.cases(all.covs),]
  summary (data)
  emplogit<-function(Y,N){
    top=Y*N+0.5
    bottom=N*(1-Y)+0.5
    return(log(top/bottom))
  }
  data$total_c <- as.numeric(data$total_c)
  data$prop_vac_logit<-emplogit(data$prop_vac,data$total_c)
  Y<-data$prop_vac_logit
  X<-all.covs[, -c(1, ncol(all.covs))]
  
  ##################################################################### data
  nfolds=5
  
  covariate.names<-colnames(X)
  combined <- as. data.frame(cbind(Y, X))
  combined$folds <- data$fold
  combined <- combined[, !(names(combined) %in% "fold")]
  summary(combined)
  options(java.parameters = "-Xmxys.setenv(JAVA_H400g")
  #Sys.setenv(JAVA_HOME = "C:/Program Files/Java/jre1.8.0_421")
  Sys.setenv(JAVA_HOME = "C:/Program Files/Java/jre1.8.0_451")
  library(h2o)
  h2o.init(ip = 'localhost', port = 54321,nthreads = -1, max_mem_size = "64g")
  h2odat<-as.h2o(combined)
  # install.packages("bit64")
  #  install.packages("bit64", dependencies = TRUE)
  #library(data.table)
  #packageVersion("bit64")
  ########### GBM ########################################################### 
  sample_rate=0.5
  col_sample_rate=0.7
  stopping_metric = "deviance"
  stopping_rounds=5
  stopping_tolerance=0
  min_rows=10
  ntrees=10000
  max_depth=4
  learn_rate=0.007
  
  gbm <- h2o.gbm(y = "Y", x = covariate.names,
                 
                 distribution = "gaussian",
                 training_frame =  h2odat,keep_cross_validation_predictions=TRUE,
                 ntrees=ntrees, max_depth=max_depth, learn_rate=learn_rate,
                 sample_rate=sample_rate,col_sample_rate=col_sample_rate,fold_column="folds",score_each_iteration=TRUE,
                 stopping_rounds=stopping_rounds, stopping_tolerance=stopping_tolerance,min_rows=min_rows,stopping_metric="deviance")
  print (gbm)
  
  ntrees<-gbm@parameters$ntrees 
  tmp = h2o.predict(gbm, newdata=h2odat)
  fit.gbm=as.vector(tmp$predict) 
  print(ntrees)
  
  # Save the model to the specified path##
  h2o.saveModel(gbm, path = "/home/backup/StackedGeneraliser/models", force = TRUE)
  

  cv.gbm<-rep(NA,nrow(data))
  for(i in 1:nfolds){
    train <- h2odat[h2odat$folds!=i,]
    gbm <- h2o.gbm(y = "Y", x = covariate.names,
                   distribution ="gaussian",
                   training_frame = train,
                   ntrees=ntrees, max_depth=max_depth, learn_rate=learn_rate,
                   sample_rate=sample_rate,col_sample_rate=col_sample_rate,min_rows=min_rows)
    valid <- h2odat[h2odat$folds==i,]
    tmp <- h2o.predict(gbm, newdata=valid)
    cv.gbm[as.logical(as.vector((h2odat$folds==i)))]=as.vector(tmp$predict)
    tmp<-c()
  }	


  #### MARS  ################################################################## 
  library(earth)
  mars_covs<-covariate.names
  formula<-as.formula(paste0('Y~',paste0(mars_covs,collapse='+')))
  
  mars<-earth(formula,data=combined,penalty=-1)
  fit.mars<-mars$fitted.values
  head(fit.mars)
  mars_store = mars
  
  cv.mars<-rep(NA,nrow(data))
  for(i in 1:nfolds){
    train <- combined[combined$folds!=i,]
    mars<-earth(formula,data=train,penalty=-1)
    
    valid <- combined[combined$folds==i,]
    cv.mars[as.logical(as.vector((combined$folds==i)))] = predict(mars, newdata=valid)
  }	
  plot(cv.mars)
  print(mars_store)
  summary(cv.mars)
  
#####  GAM  ##############################################################################
  library(mgcv)
  data_a<- as. data.frame(bcgcov)
  data<-data_a
  all.covs <- subset(data, select = -c(cluster_id, iso_a3, vaccinated_c,unvaccinated_c, year, total_c,longitude,latitude, prop_vac))
  all.covs<-as.data.frame(all.covs)
  data<-data[complete.cases(all.covs),]
  all.covs<-all.covs[complete.cases(all.covs),]
  summary (data)
  samp<-sample(1:nrow(data),nrow(data),replace=F)
  data<-data[samp,]
  all. covs<-all.covs[samp,]
  emplogit<-function(Y,N){
    top=Y*N+0.5
    bottom=N*(1-Y)+0.5
    return(log(top/bottom))
  }
  data$total_c <- as.numeric(data$total_c)
  data$prop_vac_logit<-emplogit(data$prop_vac,data$total_c)
  Y<-data$prop_vac_logit
  X<-all. covs
  nfolds=5
  library(cvTools)
  folds<-cvFolds(nrow(data),nfolds)
  covariate.names<-colnames(X)
  combined<-as. data.frame(cbind(Y,X))
  combined$folds=folds$which
  gam_covs <- covariate .names
  formula <- as.formula(paste0('Y ~ ', paste0('s(', gam_covs, ')', collapse = ' + ')))
  gam <- gam(formula, data = combined)
  fit.gam <- gam$fitted.values
  gam_store <- gam
  
  cv.gam <- rep(NA, nrow(combined))
  for (i in 1:nfolds) {
    train <- combined[combined$folds != i, ]
    gm_cv <- gam(formula, data = train)
    valid <- combined[combined$folds == i, ]
    cv.gam[combined$folds == i] <- predict(gm_cv, newdata = valid)
  }

  
  plot(cv.gam, as.vector(combined$Y))
  cor(plogis(cv.gam), plogis(as.vector(combined$Y)))
  plot(plogis(cv.gam), plogis(as.vector(combined$Y)))
  plot(gam_store)
  library(ggplot2)
  
  combined$predicted <- plogis(predict(cv.gam, type = "link"))
  
  cv.all<-cbind(cv.gbm,cv.gam,cv.mars)
  
  ############################################################################################ 
  ###4. stack rasters for prediction target years
 #######################################################################################
 
  get_raster_for_year <- function(rlist, year) {
    year_char <- as.character(year)
    if (year_char %in% names(rlist)) {
      return(rlist[[year_char]])
    } else {
      stop(paste("Year", year, "not found in raster list"))
    }
  }
  ###1990
   rs_1990 <- stack(
    rs2,
    get_raster_for_year(ebf_normalized_all, 2000),
    get_raster_for_year(pm_all, 1998),
    get_raster_for_year(edu_normalized_all, 2000),
    get_raster_for_year(usan_normalized_all, 2000),
    get_raster_for_year(uwat_normalized_all, 2000)
  )
  names(rs_1990) <- c(
    "TMP", "EH", "PCP", "ALT", "ACC", "ACH", "POP", "RDI", "NO2", "pm25", "edu", "ebf", "usan", "uwat"
  )
  names(rs_1990)
  plot(rs_1990)
  head(rs_1990)
  
  ###2000
  rs_2000 <- stack(
    rs2,
    get_raster_for_year(ebf_normalized_all, 2000),
    get_raster_for_year(pm_all, 2000),
    get_raster_for_year(edu_normalized_all, 2000),
    get_raster_for_year(usan_normalized_all, 2000),
    get_raster_for_year(uwat_normalized_all, 2000)
  )
  names(rs_2000) <- c(
    "TMP", "EH", "PCP", "ALT", "ACC", "ACH", "POP", "RDI", "NO2", "pm25", "edu", "ebf", "usan", "uwat"
  )
  names(rs_2000)
  plot(rs_2000)
 
  
  ### 2010
  rs_2010 <- stack(
    rs2,
    get_raster_for_year(ebf_normalized_all, 2010),
    get_raster_for_year(pm_all, 2010),
    get_raster_for_year(edu_normalized_all, 2010),
    get_raster_for_year(usan_normalized_all, 2010),
    get_raster_for_year(uwat_normalized_all, 2010)
  )
  names(rs_2010) <- c(
    "TMP", "EH", "PCP", "ALT", "ACC", "ACH", "POP", "RDI", "NO2", "pm25", "edu", "ebf", "usan", "uwat"
  )
  names(rs_2010)
  plot(rs_2010)
  head(rs_2010)
  
  ###2020
  rs_2020 <- stack(
    rs2,
    get_raster_for_year(ebf_normalized_all, 2019),
    get_raster_for_year(pm_all, 2020),
    get_raster_for_year(edu_normalized_all, 2015),
    get_raster_for_year(usan_normalized_all, 2017),
    get_raster_for_year(uwat_normalized_all, 2017)
  )
  names(rs_2020) <- c(
    "TMP", "EH", "PCP", "ALT", "ACC", "ACH", "POP", "RDI", "NO2", "pm25", "edu", "ebf", "usan", "uwat"
  )
  names(rs_2020)
  plot(rs_2020)
  head(rs_2020)
  
  ###2022
  
  rs_2022 <- stack(
    rs2,
    get_raster_for_year(ebf_normalized_all, 2019),
    get_raster_for_year(pm_all, 2021),
    get_raster_for_year(edu_normalized_all, 2015),
    get_raster_for_year(usan_normalized_all, 2017),
    get_raster_for_year(uwat_normalized_all, 2017)
  )
  
  names(rs_2020) <- c(
    "TMP", "EH", "PCP", "ALT", "ACC", "ACH", "POP", "RDI", "NO2", "pm25", "edu", "ebf", "usan", "uwat"
  )
  names(rs_2020)
  plot(rs_2020)
  head(rs_2020)
  ### 4.1. make in sample prediction for each year##########
  library(raster)
  library(h2o)
  library(earth)
  
  years <- c(1990, 2000, 2010, 2020, 2022)
  var_names <- c("TMP", "EH", "PCP", "ALT", "ACC", "ACH", "POP", "RDI", "NO2", "pm25", "edu", "ebf", "usan", "uwat")
  inv_logit <- function(x) {
    1 / (1 + exp(-x))
  }
  path_output <- file.path(getwd(),"OUTPUT")
  for (yr in years) {
    message(paste("Processing year:", yr))
    
    # Dynamically get the raster stack for the year
    rs <- get(paste0("rs_", yr))
    names(rs) <- var_names
    rs <- stack(rs)
    #####################################################################################################################################
    ##### Make in-sample predictive surfaces at 5 x 5 km pixel using each of submodels
   ###################################################################################################################################
    ## ---- GAM Prediction ----
    gam_pred <- predict(rs, gam_store, type = "response")
    gam_prob <- calc(gam_pred, fun = inv_logit)
    writeRaster(gam_prob, filename = file.path(path_output, paste0("gam_insample_pred_", yr, ".tif")), overwrite = TRUE)
    
    ## ---- GBM Prediction (H2O) ----
    raster_df <- as. data.frame(rs, xy = TRUE, na.rm = FALSE)
    predictor_values <- raster_df[, var_names]
    valid_mask <- complete.cases(predictor_values)
    predictors_df <- predictor_values[valid_mask, ]
    predictors_h2o <- as.h2o(predictors_df)
    pred_h2o <- h2o.predict(gbm, newdata = predictors_h2o)
    pred_vals <- as.vector(as.data.frame(pred_h2o)[, 1])
    raster_df$gbm_pred <- NA
    raster_df$gbm_pred[valid_mask] <- pred_vals
    gbm_pred <- rasterFromXYZ(raster_df[, c("x", "y", "gbm_pred")], crs = crs(rs))
    gbm_prob <- calc(gbm_pred, fun = inv_logit)
    writeRaster(gbm_prob, filename = file.path(path_output, paste0("gbm_insample_pred_", yr, ".tif")), overwrite = TRUE)
    plot(gbm_prob, main = paste("GBM Prediction", yr))
    
    ## ---- MARS Prediction ----
    coords <- raster_df[, c("x", "y")]
    predictors <- raster_df[, var_names]
    valid_mask <- rowSums(!is.na(predictors)) > 0
    raster_df$mars_pred <- NA
    if (any(valid_mask)) {
      mars_preds <- predict(mars, newdata = predictors[valid_mask, ], type = "response")
      raster_df$mars_pred[valid_mask] <- mars_preds
    }
    mars_pred <- rasterFromXYZ(raster_df[, c("x", "y", "mars_pred")], crs = crs(rs))
    mars_prob <- calc(mars_pred, fun = inv_logit)
    writeRaster(mars_prob, filename = file.path(path_output, paste0("mars_insample_pred_", yr, ".tif")), overwrite = TRUE)
    plot(mars_prob, main = paste("MARS Prediction", yr))
  }
  
  ###******************************************************************************************************
  # 5.Stacked Generaliser Gaussian process using INLA
############################################################# 
###Model for the main continent
##########################################################################    
   
    data_final <- cbind(data, cv.all)
    data_1 <- data_final[!data$iso_a3 %in% c("COM", "MDG"), ]
    
    data_1$cv.gbm <- plogis(data_1$cv.gbm)   # same as inverse-logit
    data_1$cv.gam <- plogis(data_1$cv.gam)
    data_1$cv.mars <- plogis(data_1$cv.mars)
    
    xy1 <- cbind( data_1$longitude,  data_1$latitude)
    cv.all_1 <- data_1[, tail(names(data_final), 3)]
    est.cov1 <- list(
      gbm = cv.all_1[[1]],
      gam = cv.all_1[[2]],
      mars = cv.all_1[[3]]
    )
    
  #define spatial and temporal mesh 
    library(INLA)
    iso3_list2 <- setdiff(unique(BCG0$iso_a3), c("COM", "MDG"))
    myarea_2 <- Africa %>% 
      filter(iso_a3 %in% iso3_list2) 
  bdry_1 <- INLA::inla.sp2segment(myarea_2)
  bdry_1$loc <- INLA::inla.mesh.map(bdry_1$loc)
  

  mesh_1<-INLA::inla.mesh.2d(loc=xy1, boundary=bdry_1, max.edge=c(0.9,6),offset=c(0.9,6),cutoff = 0.9)
  mesh_1$n
  plot(mesh_1)
  #save mesh plot
  plot.new()
  pdf(file=paste0(path_output,"/pdf/mesh.pdf"))
  if(!is.null(dev.list())) dev.off() 
  if(!is.null(dev.list())) dev.off()
  
  spde_1 = inla.spde2.matern(mesh_1,alpha=1.5)
  library(dplyr)
  
  # Create a 1-dimensional temporal mesh
  mesh1d <- inla. mesh.1d(seq(1990, 2022, by = 10),  # start and end year in data and time Points every 10 years.Yields temporal mesh nodes at 1990,2000,2010,2020, and 2022.
                         interval = c(1990, 2022),  # Interval from 1990 to 2022
                         degree = 1,  # Linear basis functions
                         boundary = c('free'))  # Free boundary conditions
    mesh1d$m
    mesh1d$loc
    
    A1.est =
      inla.spde.make.A(mesh_1, loc=as.matrix(xy1),group=data_1[,'year'],group.mesh=mesh1d)
    
    #-- Create index matrix --#
    field.indices1 =
      inla.spde.make.index("field", n.spde=mesh_1$n,n.group=mesh1d$m)
    field.group <- field.indices1$field.group
    ### setting response variable "Y" to binary
    data_1$vaccinated_c <- as.integer(data_1$vaccinated_c)
    
    data_1$total_c <- as.integer(data_1$total_c)
    
    Y= data_1$vaccinated_c
    n = data_1$total_c
    ##data_sorted <- data[order(data$vaccinated_w), ]
    maxprev <- max(Y/data_1$total_c)
    minprev <- min(Y/data_1$total_c)
    
    stack.est1 =
      inla.stack(data=list( Y = data_1$vaccinated_c,n = data_1$total_c ),
                 A=list(A1.est,1),
                 effects=
                   list(c(field.indices1),
                        c(est.cov1)
                   ),
                 tag="est", remove.unused=TRUE,compress=TRUE)
   
    
    formula <- as.formula(paste(
      "Y ~ ",
      "f(field, model=spde_1, group=field.group, control.group=list(model='ar1')) +",
      "f(gbm, model='clinear', range=c(0,1), initial=0.2) +",
      "f(gam, model='clinear', range=c(0,1), initial=0.2) +",
      "f(mars, model='clinear', range=c(0,1), initial=0.2)"
    ))
    
    stack.est1<-stack.est1
    
    #-- Call INLA and get results --#
      
     
                     mod.pred1 =   inla(formula,
                      data=inla.stack.data(stack.est1),
                      family="binomial", Ntrials = n,
                      control.predictor=list(A=inla.stack.A(stack.est1), compute=TRUE,quantiles=NULL),
                      control.compute=list(cpo=TRUE, dic=TRUE,config=TRUE),
                      keep=FALSE, verbose=TRUE,#,
                      control.inla= list(strategy = "adaptive",  
                                         int.strategy='eb',
                                         verbose=TRUE,fast=TRUE,dz=1,
                                         step.factor=0.5,
                                        stupid.search=FALSE)
                                    
    )  
                   
############################################################################                   
#### Model for Madagascar
############################################################################                 
                     library(INLA)
                     library(dplyr)
                     
                     data_final <- cbind(data, cv.all)
                     data_2 <- data_final[data_final$iso_a3 == "MDG", ]
                     
                     data_2$cv.gbm <- plogis(data_2$cv.gbm)   # same as inverse-logit
                     data_2$cv.gam <- plogis(data_2$cv.gam)
                     data_2$cv.mars <- plogis(data_2$cv.mars)
                     
                     xy2 <- cbind( data_2$longitude,  data_2$latitude)
                     cv.all_2 <- data_2[, tail(names(data_final), 3)]
                     
                     est.cov2 <- list(
                       gbm = cv.all_2[[1]],
                       gam = cv.all_2[[2]],
                       mars = cv.all_2[[3]]
                     )
                     
                     #define spatial and temporal mesh 
                     
                     iso3_list3 <- unique(BCG0$iso_a3[BCG0$iso_a3 == "MDG"])
                     myarea_3 <- Africa %>% 
                       filter(iso_a3 %in% iso3_list3) 
                     bdry_2 <- INLA::inla.sp2segment(myarea_3)
                     bdry_2$loc <- INLA::inla.mesh.map(bdry_2$loc)
                     
                     
                     mesh_2<-INLA::inla.mesh.2d(loc=xy2, boundary=bdry_2, max.edge=c(0.9,6),offset=c(0.9,6),cutoff = 0.9)
                     
                     mesh_2$n
                     plot(mesh_2)
                     #save mesh plot
                     plot.new()
                     pdf(file=paste0(path_output,"/pdf/mesh2.pdf"))
                     
                     if(!is.null(dev.list())) dev.off() 
                     if(!is.null(dev.list())) dev.off()
                     
                     
                     spde_2 = inla.spde2.matern(mesh_2,alpha=1.5)
                   
                     # Create a 1-dimensional temporal mesh
                     mesh1d <- inla.mesh.1d(seq(1990, 2022, by = 10), # start and end year in data and time Points every 10 years.Yields temporal mesh nodes at 1990,2000,2010,2020 and 2022.  
                                            interval = c(1990, 2022),  # Interval from 1990 to 2022
                                            degree = 1,  # Linear basis functions
                                            boundary = c('free'))  # Free boundary conditions
                     
                     mesh1d$m
                     
                     
                     A2.est =
                       inla.spde.make.A(mesh_2, loc=as.matrix(xy2),group=data_2[,'year'],group.mesh=mesh1d)
                     
                     #-- Create index matrix --#
                     field.indices2 =
                       inla.spde.make.index("field", n.spde=mesh_2$n,n.group=mesh1d$m)
                     field.group <- field.indices2$field.group
                     ### setting response variable "Y" to binary
                     data_2$vaccinated_c <- as.integer(data_2$vaccinated_c)
                     
                     data_2$total_c <- as.integer(data_2$total_c)
                     
                     Y= data_2$vaccinated_c
                     n = data_2$total_c
                     
                     maxprev <- max(Y/data_2$total_c)
                     minprev <- min(Y/data_2$total_c)
                     
                     stack.est2 =
                       inla.stack(data=list( Y = data_2$vaccinated_c,n = data_2$total_c ),
                                  A=list(A2.est,1),
                                  effects=
                                    list(c(field.indices2),c(est.cov2)
                                         
                                    ),
                                  tag="est", remove.unused=TRUE,compress=TRUE)
                     
                     
                     
                     formula <- as.formula(paste(
                       "Y ~ ",
                       "f(field, model=spde_2, group=field.group, control.group=list(model='ar1')) +",
                       "f(gbm, model='clinear', range=c(0,1), initial=0.2) +",
                       "f(gam, model='clinear', range=c(0,1), initial=0.2) +",
                       "f(mars, model='clinear', range=c(0,1), initial=0.2)"
                     ))
                     
                     
                     stack.est2<-stack.est2
                     
                     
                     #-- Call INLA and get results --#
                     
                     
                     mod.pred2 =   inla(formula,
                                        data=inla.stack.data(stack.est2),
                                        family="binomial", Ntrials = n,
                                        control.predictor=list(A=inla.stack.A(stack.est2), compute=TRUE,quantiles=NULL),
                                        control.compute=list(cpo=TRUE, dic=TRUE,config=TRUE),
                                        keep=FALSE, verbose=TRUE,#num.threads=5,
                                        control.inla= list(strategy = "adaptive",  
                                                           int.strategy='eb',
                                                           verbose=TRUE,fast=TRUE,dz=1,
                                                           step.factor=0.5,
                                                           stupid.search=FALSE)
                                        
                     )  
  
  save(popden, mod.pred1,mod.pred2,data_final,myarea,myarea_2,myarea_3,BCG0,cv.all,
       data,data_1,data_2, mesh_1,mesh_2, cv.all_1,cv.all_2,A1.est,A2.est,A1.pred,A2.pred,Africa,
       bdry, bdry_1,bdry_2, field.indices1,field.indices2,madagascar_shape,mesh_1,
       mesh_2,mesh1d,spde_1,spde_2,xy1,xy2, file = "my_model_objects.RData")                   
#############################################################################################################################################
######## ##### posterior sampling and coverage prediction######################
######## Main continent model)####################
       library(raster)
       library(dplyr)
       library(INLA)
       library(sf)
       library(rnaturalearth)
       
       ####### first call in sample submodel predictions for those prediction target years#############
    
       path_output <- file.path(getwd(), "OUTPUT")
       years <- c(1990, 2000, 2010, 2020, 2022)
       models <- c("gam", "gbm", "mars")
       
       # Define setup
       path_output <- file.path(getwd(), "OUTPUT")
       years <- c(1990, 2000, 2010, 2020, 2022)
       models <- c("gam", "gbm", "mars")
       all_rasters <- list()
       
       # Loop through models and read rasters
       for (model in models) {
         message(paste("Processing model:", toupper(model)))
         file_names <- paste0(model, "_insample_pred_", years, ".tif")
         file_paths <- file.path(path_output, file_names)
         for (i in seq_along(file_paths)) {
           if (file.exists(file_paths[i])) {
             r <- raster(file_paths[i])
             name <- paste0(model, "_", years[i])
             names(r) <- name
             all_rasters[[length(all_rasters) + 1]] <- r
           } else {
             warning(paste("File not found:", file_paths[i]))
           }
         }
       }
       
       final_stack <- stack(all_rasters)
       rs_90_22<-final_stack
      
       mask1<-rs_90_22[[2]]; NAvalue(mask1)=-9999
       mask1 <- mask(crop(mask1, myarea_2), myarea_2)
       plot(mask1)
       popt<-5
       plot(mask1)
       popmask1[popmask1 < popt] <- NA
       popmask1[!is.na(popmask1),] <- 1
       popmask1 <- mask(crop(popmask1, myarea_2), myarea_2)
       popmask1 <- resample(popmask1, mask1, method = "ngb")  
       plot(popmask1)
       #code to make maps
       pred_val<-getValues(mask1)
       w<-is.na(pred_val)
       index<-1:length(w)
       index<-index[!w]
       pred_locs<-xyFromCell(mask1,1:ncell(mask1))
       pred_locs<-pred_locs[!w,]
       colnames(pred_locs)<-c('longitude','latitude')
       locs_pred1 <- pred_locs
       
       if (is.matrix(locs_pred1)) {
         locs_pred1 <- as.data.frame(locs_pred1)
       }
       
       # Convert columns to numeric
       locs_pred1$longitude <- as.numeric(locs_pred1$longitude)
       locs_pred1$latitude <- as.numeric(locs_pred1$latitude)
       
       years<-mesh1d$loc
       ####prediction matrix####
       A1.pred =
         inla.spde.make.A(mesh=mesh_1, loc=as.matrix(locs_pred1))
       A1.pred <- as.matrix(A1.pred)
       dim(A1.pred)
       rs2 <- mask(crop( rs_90_22, myarea_2), myarea_2)
       plot(rs2)
       covariate_pred1 <-data.frame(raster::extract(rs2, locs_pred1))
       
       #########sampling the posterior########################
       #######################################################################################
       
       set.seed(999)
       library(INLA)
       library(sn)
       
       my_draws1 <- inla.posterior.sample(1000, mod.pred1)
       prednames1 <- my_draws1[[i]]$latent # to see all latent parametrs stored as row (intercept,field,gbm,gam,mars)
       prednames1<-as.data.frame(prednames1)
       
       ##############################################################################
       ###prediction using linear predictor#####################################################################
       # Define the mapping of loop indices to years
       
       sp_per_year <- mesh_1$n
       
       all_pred1 <- vector("list", 5) 
       
       for (years in 1:5) {
         # Field index for that year
         start_index <- (years - 1) * sp_per_year + 1
         end_index <- years * sp_per_year
         f.indexes <- start_index:end_index
         
         # Map loop index to actual year value
         year_sequence <- c(1990, 2000, 2010, 2020, 2022)
         year_label <- as.character(year_sequence[years])
         
         #  Subset covariates for this year only 
         year_cols <- grep(paste0("_", year_label, "$"), colnames(covariates_pred1), value = TRUE)
         X_year <- covariates_pred1[, year_cols, drop = FALSE]
         model_names <- sub(paste0("_", year_label), "", year_cols)  # get 'gam', 'gbm', etc.
         
         #Initialize matrix to hold prediction samples
         num_draws1 <- length(my_draws1)
         pred1 <- matrix(NA, nrow = nrow(A1.pred), ncol = num_draws1)
         
         for (i in 1:num_draws1) {
           samp <- my_draws1[[i]]
           
           # Spatial field extraction
           field <- samp$latent[grep("field", rownames(samp$latent)), ]
           field1 <- matrix(as.numeric(field[f.indexes]), ncol = 1)
           
           # Match beta for each covariate type
           k <- ncol(X_year)
           beta <- numeric(k)
           for (j in 1:k) {
             model_prefix <- model_names[j]
             beta_idx <- grep(paste0("Beta for ", model_prefix), names(samp$hyperpar), ignore.case = TRUE)
             beta[j] <- if (length(beta_idx) > 0) samp$hyperpar[beta_idx[1]] else 0
           }
           
           # Intercept
           intercept_idx <- grep("intercept", rownames(samp$latent), ignore.case = TRUE)
           intercept <- if (length(intercept_idx) > 0) samp$latent[intercept_idx[1]] else 0
           
           #  Linear predictor
           linpred <- X_year %*% beta
           lp <- intercept + linpred + drop(A1.pred %*% field1)
           pred1[, i] <- lp
         }
         
         # Apply inverse-logit to get predicted probabilities
         pred1 <- plogis(pred1)
         all_pred1[[years]] <- pred1
       }
       
       
       all_pred1_c <- do.call(rbind, all_pred1)
       
       
       # Create the predictions for quantiles
       pred_2.5 <- apply(all_pred1_c, 1, function(x) quantile(x, probs = c(0.025), na.rm = TRUE))
       
       pred_sd <- apply(all_pred1_c, 1, sd)
       pred_mean <- apply(all_pred1_c, 1, function(x) mean(x, na.rm = TRUE))
       pred_med <- apply(all_pred1_c, 1, function(x) quantile(x, probs = c(0.5), na.rm = TRUE))
       pred_975 <- apply(all_pred1_c, 1, function(x) quantile(x, probs = c(0.975), na.rm = TRUE))
       year <- rep(1:5, each = nrow(A1.pred))
       # Saving predictive maps as raster files
       predinput <- list(pred_2.5, pred_sd, pred_mean,pred_med, pred_975,year)
       prednames <- as.list(c("LUI", "Prob_sd", "Prob_mean","prob_med", "UUI","year"))
       
  
       
       out <- list()
       
       # Total number of years
       n_years <- 5
       year_labels <- c(1990, 2000, 2010, 2020, 2022)
       # Number of spatial units per year
       n_spatial <- nrow(A1.pred)
       
       # Assign names to predinput
       names(predinput) <- c("LUI", "Prob_sd", "Prob_mean", "prob_med", "UUI", "year")
       
       # Loop over prediction metrics (exclude "year")
       for (j in 1:(length(predinput) - 1)) {
         metric_name <- names(predinput)[j]
         metric_values <- predinput[[j]]
         
         for (year in 1:n_years) {
           # Subset values for this year
           index <- which(predinput$year == year)
           values_year <- metric_values[index]
           
           # Prepare raster with predicted values
           pred_val[!w] <- values_year
           out_raster <- setValues(mask1, pred_val)
           
           # Apply the resampled population mask
           out_raster <- mask(out_raster, popmask1)
           
           # File name with metric and year label
           file_name <- paste0(metric_name,"_all_", "_Year_", year_labels[year], ".tif")
           
           # Save raster
           writeRaster(out_raster, file.path(path_output, "prevalence", file_name), overwrite = TRUE)
           
         }
       }
       
       
       
    ############################################################################################################################################
    ##### Model 2 (Madagascar)#######
    library(raster)
    library(dplyr)
    library(INLA)
    library(sf)
    library(rnaturalearth)
    madagascar_shape <- ne_countries(scale = "medium", country = "Madagascar", returnclass = "sf")   
    mask2<-rs_90_22[[2]]; NAvalue(mask2)=-9999
    mask_madagascar <- mask(crop(mask2, madagascar_shape), madagascar_shape)
    mask_madagascar <- mask(crop(mask2, myarea_3), myarea_3)
    plot(mask_madagascar)
    popt<-5
    plot(mask_madagascar)
    
    popmask[popmask < popt] <- NA
    popmask[!is.na(popmask),] <- 1
    popmask <- crop(popmask, mask_madagascar)  # Crop popmask to shape extent
    popmask <- resample(popmask, mask_madagascar, method = "ngb")  
    plot(popmask)
    #code to make maps
    pred_val<-getValues(mask_madagascar)
    w<-is.na(pred_val)
    index<-1:length(w)
    index<-index[!w]
    pred_locs<-xyFromCell(mask_madagascar,1:ncell(mask_madagascar))
    pred_locs<-pred_locs[!w,]
    colnames(pred_locs)<-c('longitude','latitude')
    locs_pred2 <- pred_locs
   
    if (is.matrix(locs_pred2)) {
      locs_pred2 <- as. data.frame(locs_pred2)
    }
    
    # Convert columns to numeric
    locs_pred2$longitude <- as.numeric(locs_pred2$longitude)
    locs_pred2$latitude <- as.numeric(locs_pred2$latitude)
    years<-mesh1d$loc
    ####prediction matrix####
    A2.pred =
      inla.spde.make.A(mesh=mesh_2, loc=as.matrix(locs_pred2))
    A2.pred <- as.matrix(A2.pred)
    dim(A2.pred)
    rs2_madagascar <- mask(crop(rs_90_22, madagascar_shape), madagascar_shape)
    plot(rs2_madagascar)
    covariate_pred2 <-data.frame(raster::extract(rs2_madagascar, pred_locs))
    covariate_pred2 <-data.frame(raster::extract(rs2_madagascar, pred_locs))
    
    #########sampling the posterior########################
    #######################################################################################
    
    set.seed(999)
    library(INLA)
    library(sn)
    
   my_draws2 <- inla. posterior.sample(1000, mod.pred2)
   prednames2 <- my_draws2[[i]]$latent # to see all latent parameters stored as row (intercept, field,gbm, gam, mars)
   prednames2<-as. data.frame(prednames2)
   
   ##############################################################################
   ### raster for each year#####################################################################
     sp_per_year <- mesh_2$n
     all_pred2 <- vector("list", 5)
     
     for (years in 1:5) {
       # Field index for that year
       start_index <- (years - 1) * sp_per_year + 1
       end_index <- years * sp_per_year
       f.indexes <- start_index:end_index
       
       # Map loop index to actual year value
       year_sequence <- c(1990, 2000, 2010, 2020, 2022)
       year_label <- as.character(year_sequence[years])
       
       #  Subset covariates for this year only 
       year_cols <- grep(paste0("_", year_label, "$"), colnames(covariate_pred2), value = TRUE)
       X_year <- covariate_pred2[, year_cols, drop = FALSE]
       model_names <- sub(paste0("_", year_label), "", year_cols)  # get 'gam', 'gbm', etc.
       
       #Initialize matrix to hold prediction samples
       num_draws2 <- length(my_draws2)
       pred2 <- matrix(NA, nrow = nrow(A2.pred), ncol = num_draws2)
       
       for (i in 1:num_draws2) {
         samp <- my_draws2[[i]]
         
         # Spatial field extraction
         field <- samp$latent[grep("field", rownames(samp$latent)), ]
         field1 <- matrix(as.numeric(field[f.indexes]), ncol = 1)
         
         # Match beta for each covariate type
         k <- ncol(X_year)
         beta <- numeric(k)
         for (j in 1:k) {
           model_prefix <- model_names[j]
           beta_idx <- grep(paste0("Beta for ", model_prefix), names(samp$hyperpar), ignore.case = TRUE)
           beta[j] <- if (length(beta_idx) > 0) samp$hyperpar[beta_idx[1]] else 0
         }
         
         # Intercept
         intercept_idx <- grep("intercept", rownames(samp$latent), ignore.case = TRUE)
         intercept <- if (length(intercept_idx) > 0) samp$latent[intercept_idx[1]] else 0
         
         #  Linear predictor
         X_year <- as.matrix(X_year)
         linpred <- X_year %*% beta
         lp <- intercept + linpred + drop(A2.pred %*% field1)
         pred2[, i] <- lp
       }
       
       # Apply inverse-logit to get predicted probabilities
       pred2 <- plogis(pred2)
       all_pred2[[years]] <- pred2
     }
     
      
     all_pred2_c <- do.call(rbind, all_pred2)
     
      # Create the predictions for quantiles
     pred_2.5 <- apply(all_pred2_c, 1, function(x) quantile(x, probs = c(0.025), na.rm = TRUE))
     pred_sd <- apply(all_pred2_c, 1, sd)
     pred_mean <- apply(all_pred2_c, 1, function(x) mean(x, na.rm = TRUE))
     pred_med <- apply(all_pred2_c, 1, function(x) quantile(x, probs = c(0.5), na.rm = TRUE))
     pred_975 <- apply(all_pred2_c, 1, function(x) quantile(x, probs = c(0.975), na.rm = TRUE))
     year <- rep(1:5, each = nrow(A2.pred))
     # Saving predictive maps as raster files
     predinput <- list(pred_2.5, pred_sd, pred_mean,pred_med, pred_975,year)
     prednames <- as.list(c("LUI", "Prob_sd", "Prob_mean","prob_med", "UUI","year"))
    
     out <- list()
     
     # Total number of years
     n_years <- 5
     year_labels <- c(1990, 2000, 2010, 2020, 2022)
     n_spatial <- nrow(A2.pred) # Number of spatial units per year
     names(predinput) <- c("LUI", "Prob_sd", "Prob_mean", "prob_med", "UUI", "year")
     
     # Loop over prediction metrics (exclude "year")
     for (j in 1:(length(predinput) - 1)) {
       metric_name <- names(predinput)[j]
       metric_values <- predinput[[j]]
       
       for (year in 1:n_years) {
         index <- which(predinput$year == year)
         values_year <- metric_values[index]
         pred_val[!w] <- values_year
         out_raster <- setValues(mask_madagascar, pred_val)
         out_raster <- mask(out_raster, popmask)
         file_name <- paste0(metric_name,"_MAD_", "_Year_", year_labels[year], ".tif")
         writeRaster(out_raster, file.path(path_output, "prevalence", file_name), overwrite = TRUE)
         
       }
     }
     #############################################################################################
     ###Merging raster estimates from both models for subsequent analysis
    ###################################################################### #################
     
     library(raster)
     path_output <- file.path(getwd(), "OUTPUT")
     
     metrics <- c("LUI", "Prob_sd", "Prob_mean", "prob_med", "UUI")
     year_labels <- c(1990, 2000, 2010, 2020, 2022)
     input_dir <- file.path(path_output, "prevalence")
     output_dir <- file.path(path_output, "prevalence_merged")
     
     # Create output directory 
     dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
     
     for (metric in metrics) {
       for (year in year_labels) {
         
         # Construct input file paths
         r1_path <- file.path(input_dir, paste0(metric, "_all__Year_", year, ".tif"))
         r2_path <- file.path(input_dir, paste0(metric, "_MAD__Year_", year, ".tif"))
         
         # Check if both files exist
         if (file.exists(r1_path) && file.exists(r2_path)) {
           
           # Load rasters
           r1 <- raster(r1_path)
           r2 <- raster(r2_path)
           
           # Align r2 to r1 grid (nearest neighbor method; change to "bilinear" if continuous data)
           r2_aligned <- projectRaster(r2, r1, method = "ngb")
           
           # Merge rasters
           merged_raster <- merge(r1, r2_aligned)
           
           # Output file path
           output_path <- file.path(output_dir, paste0(metric, "_Year_", year, ".tif"))
           
           # Save merged raster
           writeRaster(merged_raster, output_path, overwrite = TRUE)
           
           cat("Merged and saved:", output_path, "\n")
           
         } else {
           warning(paste("Missing files for", metric, "in year", year))
         }
       }
     } 
    
     
 ########################################################################################   
###number of unvaccinated children at the raster level
  ############################################################################################## 
     ###1990###
     path_input <- file.path(getwd(), "INPUT")
     estimates <- c("LUI_Year_1990.tif", "Prob_mean_Year_1990.tif","prob_med_Year_1990.tif", "UUI_Year_1990.tif")
   raster.list <- file.path(path_output, "prevalence_merged", estimates)
   #extract raster data
   rasterls<-list()
   for (i in 1:length(raster.list)){
     rasterls[[i]]<-raster::raster(raster.list[[i]])
   }
   b <- raster::brick(rasterls)
   bcg_prev_1990 <-raster::stack(b[[1]], b[[2]],b[[3]], b[[4]])
   layer_names <- names(bcg_prev_1990)
   print(layer_names)
   pop_1990 <- raster::raster(paste0(path_input,"/ihme_corrected_worldpop_2_to_4_3_1990.tif"))
   
   pop_1990 <- resample(pop_1990 , bcg_prev_1990, method = "ngb")
   unvac_prev<- (1 - bcg_prev_1990)
   unvaccinated_count_1990 <- (pop_1990*unvac_prev) 
   unvaccinated_count_1990 <- round(unvaccinated_count_1990)
   names(unvaccinated_count_1990) <- c("LUI_2.5", "mean", "median", "UUI_97.5")
   writeRaster(unvaccinated_count_1990[[1]], paste0(path_output,'/','tif/unvac_LUI_1990.tif'), overwrite=TRUE)
   writeRaster(unvaccinated_count_1990[[2]], paste0(path_output,'/','tif/unvac_mean_1990.tif'), overwrite=TRUE)
   writeRaster(unvaccinated_count_1990[[2]], paste0(path_output,'/','tif/unvac_median_1990.tif'), overwrite=TRUE)
   writeRaster(unvaccinated_count_1990[[3]], paste0(path_output,'/','tif/unvac_UUI_1990.tif'), overwrite=TRUE)
   
   ###2000###
   estimates <- c("LUI_Year_2000.tif", "Prob_mean_Year_2000.tif","Prob_med_Year_2000.tif", "UUI_Year_2000.tif")
   raster.list <- file.path(path_output, "prevalence_merged", estimates)
   #extract raster data
   rasterls<-list()
   for (i in 1:length(raster.list)){
     rasterls[[i]]<-raster::raster(raster.list[[i]])
   }
   b <- raster::brick(rasterls)
   bcg_prev_2000 <-raster::stack(b[[1]], b[[2]],b[[3]],b[[4]])
   layer_names <- names(bcg_prev_2000)
   print(layer_names)
   pop_2000 <- raster::raster(paste0(path_input,"/ihme_corrected_worldpop_2_to_4_3_2000.tif"))
   
   pop_2000 <- resample(pop_2000 , bcg_prev_2000, method = "ngb")
   unvac_prev<- (1 - bcg_prev_2000)
   unvaccinated_count_2000 <- (pop_2000*unvac_prev) 
   
   unvaccinated_count_2000 <- round(unvaccinated_count_2000)
   names(unvaccinated_count_2000) <- c("LUI_2.5", "mean", "median" ,"UUI_97.5")
   writeRaster(unvaccinated_count_2000[[1]], paste0(path_output,'/','tif/unvac_LUI_2000.tif'), overwrite=TRUE)
   writeRaster(unvaccinated_count_2000[[2]], paste0(path_output,'/','tif/unvac_mean_2000.tif'), overwrite=TRUE)
   writeRaster(unvaccinated_count_2000[[2]], paste0(path_output,'/','tif/unvac_median_2000.tif'), overwrite=TRUE)
   writeRaster(unvaccinated_count_2000[[3]], paste0(path_output,'/','tif/unvac_UUI_2000.tif'), overwrite=TRUE)
   
   
   ###2010###
   estimates <- c("LUI_Year_2010.tif", "Prob_mean_Year_2010.tif","Prob_med_Year_2010.tif", "UUI_Year_2010.tif")
   raster.list <- file.path(path_output, "prevalence_merged", estimates)
   #extract raster data
   rasterls<-list()
   for (i in 1:length(raster.list)){
     rasterls[[i]]<-raster::raster(raster.list[[i]])
   }
   b <- raster::brick(rasterls)
   bcg_prev_2010 <-raster::stack(b[[1]], b[[2]],b[[3]],b[[4]])
   layer_names <- names(bcg_prev_2010)
   print(layer_names)
   pop_2010 <- raster::raster(paste0(path_input,"/ihme_corrected_worldpop_2_to_4_3_2010.tif"))
   
   pop_2010 <- resample(pop_2010 , bcg_prev_2010, method = "ngb")
   unvac_prev<- (1 - bcg_prev_2010)
   unvaccinated_count_2010 <- (pop_2010*unvac_prev) 
   
   unvaccinated_count_2010 <- round(unvaccinated_count_2010)
   names(unvaccinated_count_2010) <- c("LUI_2.5","mean","median","UUI_97.5")
   writeRaster(unvaccinated_count_2010[[1]], paste0(path_output,'/','tif/unvac_LUI_2010.tif'), overwrite=TRUE)
   writeRaster(unvaccinated_count_2010[[2]], paste0(path_output,'/','tif/unvac_mean_2010.tif'), overwrite=TRUE)
   writeRaster(unvaccinated_count_2010[[2]], paste0(path_output,'/','tif/unvac_median_2010.tif'), overwrite=TRUE)
   writeRaster(unvaccinated_count_2010[[3]], paste0(path_output,'/','tif/unvac_UUI_2010.tif'), overwrite=TRUE)
   
   ###2020###
   estimates <- c("LUI_Year_2020.tif", "Prob_mean_Year_2020.tif","Prob_med_Year_2020.tif", "UUI_Year_2020.tif")
   raster.list <- file.path(path_output, "prevalence_merged", estimates)
   #extract raster data
   rasterls<-list()
   for (i in 1:length(raster.list)){
     rasterls[[i]]<-raster::raster(raster.list[[i]])
   }
   b <- raster::brick(rasterls)
   bcg_prev_2020 <-raster::stack(b[[1]], b[[2]],b[[3]],b[[4]])
   layer_names <- names(bcg_prev_2020)
   print(layer_names)
   pop_2020 <- raster::raster(paste0(path_input,"/ihme_corrected_worldpop_2_to_4_3_2020.tif"))
   
   pop_2020 <- resample(pop_2020 , bcg_prev_2020, method = "ngb")
   unvac_prev<- (1 - bcg_prev_2020)
   unvaccinated_count_2020 <- (pop_2020*unvac_prev) 
   
   unvaccinated_count_2020 <- round(unvaccinated_count_2020)
   names(unvaccinated_count_2020) <- c("LUI_2.5","mean","median","UUI_97.5")
   writeRaster(unvaccinated_count_2020[[1]], paste0(path_output,'/','tif/unvac_LUI_2020.tif'), overwrite=TRUE)
   writeRaster(unvaccinated_count_2020[[2]], paste0(path_output,'/','tif/unvac_mean_2020.tif'), overwrite=TRUE)
   writeRaster(unvaccinated_count_2020[[2]], paste0(path_output,'/','tif/unvac_median_2020.tif'), overwrite=TRUE)
   writeRaster(unvaccinated_count_2020[[3]], paste0(path_output,'/','tif/unvac_UUI_2020.tif'), overwrite=TRUE)
  
   ###2022###
   estimates <- c("LUI_Year_2022.tif", "Prob_mean_Year_2022.tif","Prob_med_Year_2022.tif", "UUI_Year_2022.tif")
   raster.list <- file.path(path_output, "prevalence_merged", estimates)
   #extract raster data
   rasterls<-list()
   for (i in 1:length(raster.list)){
     rasterls[[i]]<-raster::raster(raster.list[[i]])
   }
   b <- raster::brick(rasterls)
   bcg_prev_2022 <-raster::stack(b[[1]], b[[2]],b[[3]],b[[4]])
   layer_names <- names(bcg_prev_2022)
   print(layer_names)
   pop_2022 <- raster::raster(paste0(path_input,"/ihme_corrected_worldpop_2_to_4_3_2022.tif"))
   
   pop_2022 <- resample(pop_2022 , bcg_prev_2022, method = "ngb")
 
   unvac_prev<- (1 - bcg_prev_2022)
   unvaccinated_count_2022 <- (pop_2022*unvac_prev) 
   
   unvaccinated_count_2022 <- round(unvaccinated_count_2022)
   names(unvaccinated_count_2022) <- c("LUI_2.5","mean","median","UUI_97.5")
   
   writeRaster(unvaccinated_count_2022[[1]], paste0(path_output,'/','tif/unvac_LUI_2022.tif'), overwrite=TRUE)
   writeRaster(unvaccinated_count_2022[[2]], paste0(path_output,'/','tif/unvac_mean_2022.tif'), overwrite=TRUE)
   writeRaster(unvaccinated_count_2022[[2]], paste0(path_output,'/','tif/unvac_median_2022.tif'), overwrite=TRUE)
   writeRaster(unvaccinated_count_2022[[3]], paste0(path_output,'/','tif/unvac_UUI_2022.tif'), overwrite=TRUE)
   
   ######################################################################################################################
   ### Sub-national estimates for each year
   #####################################################################################################################
   ###Admin 2 coverage Estimates
   ############################################################################################################################
   library(raster)
   library(dplyr)
   library(sf)
   years <- c(1990, 2000, 2010, 2020, 2022) # List of years
   
   mycountries <- unique(myarea$iso_a3)
   mycountries <- setdiff(unique(myarea$iso_a3), "COM")##COM is excluded in the analysis
   africa_adm2 <- st_read("C:/Users/ALiyew/Desktop/data/bcg_12-23m data/africa adm2 shapefiles")
   adm2poly <- africa_adm2 %>% 
     filter(ISO %in% mycountries)
###############################################################################################################  
######################## ADM2 population weighted BCG coverage
##############################################################################################################   
   library(raster)
   library(dplyr)
   library(exactextractr)
   library(matrixStats)
   
   years <- c(1990, 2000, 2010, 2020, 2022)
   bcgprev_year_layers2 <- list()
   
   for (year in years) {
     # Get the corresponding BCG prediction raster
     raster_layer_name <- paste0("bcg_prev_", year)
     bcg_layer <- get(raster_layer_name)  # Must be a RasterStack or RasterBrick
     names(bcg_layer) <- c("LUI_2.5", "mean", "median", "UUI_97.5")
     
     # Get the corresponding population raster
     pop_layer_name <- paste0("pop_", year)
     pop_raster <- get(pop_layer_name)
     pop_raster <- projectRaster(pop_raster, bcg_layer, method = "bilinear")
     
     # Extract population-weighted medians for each layer
     weighted_medians_list <- list()
     
     for (i in 1:nlayers(bcg_layer)) {
       layer <- bcg_layer[[i]]
       
       weighted_result <- exact_extract(layer, adm2poly,
                                        weights = pop_raster,
                                        fun = function(values, coverage_fractions, weights) {
                                          valid <- !is.na(values)
                                          
                                          if (sum(valid) == 0) {
                                            return(NA_real_)
                                          }
                                          
                                          matrixStats::weightedMedian(values[valid], w = weights[valid])
                                        },
                                        progress = FALSE
       )
       weighted_medians_list[[i]] <- weighted_result
     }
     # Combine into a data frame
     weighted_df <- as.data.frame(do.call(cbind, weighted_medians_list))
     colnames(weighted_df) <- names(bcg_layer)
     
     # Join with admin2 polygons
  
     Af2_bcg <- adm2poly
     Af2_bcg <- cbind(Af2_bcg, weighted_df) %>%
       mutate(Year = year) %>%
       dplyr::select(ISO, Name_0, Name_1, Name_2, Year, LUI_2.5, mean, median, UUI_97.5, geometry)
     
     # Store for this year
     bcgprev_year_layers2[[as.character(year)]] <- Af2_bcg
   }
   
   
   # Loop through each year's raster object and rename the target layers
   bcgprev_year_layers2 <- lapply(bcgprev_year_layers2, function(r) {
     layer_names <- names(r)
     
     # Replace names
     layer_names[layer_names == "LUI_2.5"] <- "LUI"
     layer_names[layer_names == "UUI_97.5"] <- "UUI"
     
     # Assign back the new names
     names(r) <- layer_names
     
     return(r)
   })
   lapply(bcgprev_year_layers2, names)
   
   ### export shapefile for each year
   for (y in 1:length(bcgprev_year_layers2)) {
     year_data <- bcgprev_year_layers2[[y]]
     year <- year_data$Year[1]  # Access 'Year' column directly
     year_data_sf <- st_as_sf(year_data) # Ensure it's an sf object
     gpkg_file <- file.path(path_output, 'shp', paste0('bcgprev_adm2_', year, '.gpkg')) # Set output filename with .gpkg extension
     st_write(year_data_sf, dsn = gpkg_file, layer = paste0("bcgprev_adm2_", year), delete_layer = TRUE) # Write to GeoPackage (supports NA properly)
   }
  ####### #### Admin2 coverage CSV data ##########################################################
   
   years_csv2 <- list()  # Initialize the list to store results
  for (y in 1:length(bcgprev_year_layers2)) {
     bcgprev_data <- bcgprev_year_layers2[[y]]
     year <- bcgprev_data$Year[1]  # Access 'Year' column directly
     # Extract the data frame for prevalence data (remove geometry using st_drop_geometry)
     af2_bcg_data <- st_drop_geometry(bcgprev_data)  # This removes geometry, leaving attribute data
     prevadm2 <- data.frame(af2_bcg_data)
     
     # Add the processed data to the list
     years_csv2[[y]] <- prevadm2
   }
   
   final_data <- do.call(cbind, years_csv2) # Column-bind all data frames in the list
   final_data_r <- do.call(rbind, years_csv2)
   write.csv(final_data, file.path(path_output, "csv", "prevadm2_all.csv"), row.names = FALSE) # Save the combined data as a CSV
   
   ###############################################################################################################
   ############ estimating BCG coverage change from 1990 to 2022 ################################# #
   ###############################################################################################################
   
  #####changes at Raster level
   # Load raster files)
   bcg_1990 <- raster(file.path(path_output, "prevalence_merged", "Prob_med_Year_1990.tif"))
   bcg_2022 <- raster(file.path(path_output, "prevalence_merged", "Prob_med_Year_2022.tif"))
   
   
   bcg_change <- bcg_2022 - bcg_1990 # Calculate the change in coverage
   writeRaster(bcg_change, file.path(path_output, "prevalence_merged", "median_bcg_change_1990_2022.tif"), format = "GTiff", overwrite = TRUE)
   summary(bcg_change) # Summary statistics of the change
   plot(bcg_change, main = "BCG Coverage Change (1990 to 2020") # Plot the change
   
   ###### change at adm2 level
   ### For coverage change at admin 2, use population-weighted admin 2 coverage estimates for 1990 and 2022
   
   # Extract the 1990 and 2022 layers
   bcg_1990 <- bcgprev_year_layers2[["1990"]]
   bcg_2022 <- bcgprev_year_layers2[["2022"]]
   
   ### Select relevant columns (excluding geometry to simplify the join)
   bcg_1990_df <- bcg_1990 %>% 
     st_drop_geometry() %>% 
     dplyr::select(ISO, Name_0, Name_1, Name_2, LUI_2.5_1990 = LUI_2.5, mean_1990 = mean,
                   median_1990 = median, UUI_97.5_1990 = UUI_97.5)
   
   bcg_2022_df <- bcg_2022 %>%
     dplyr::select(ISO, Name_0, Name_1, Name_2, LUI_2.5_2022 = LUI_2.5, mean_2022 = mean,
                   median_2022 = median, UUI_97.5_2022 = UUI_97.5, geometry)
   
   ### Join 1990 and 2022 by polygon identifiers
   bcg_change <- bcg_2022_df %>%
     left_join(bcg_1990_df, by = c("ISO", "Name_0", "Name_1", "Name_2"))
   
  ### Calculate absolute changes
   bcg_change <- bcg_change %>%
     mutate(
       delta_LUI_2.5 = LUI_2.5_2022 - LUI_2.5_1990,
       delta_mean = mean_2022 - mean_1990,
       delta_median = median_2022 - median_1990,
       delta_UUI_97.5 = UUI_97.5_2022 - UUI_97.5_1990
     )
   ###  select only the change columns for summary or visualization
   bcg_change_summary <- bcg_change %>%
     dplyr::select(ISO, Name_0, Name_1, Name_2, starts_with("delta_"), geometry)
   ###export shapefile
   gpkg_file <- file.path(path_output, 'shp', 'bcgprev_adm2_change_1990_2022.gpkg')
   st_write(bcg_change_summary, 
            dsn = gpkg_file, 
            layer = "bcgprev_adm2_change_1990_2022", 
            delete_layer = TRUE)
  
   ###export for mapping
   st_write(bcg_cov_1990_2022,
            dsn = file.path(path_output, 'shp', 'new_bcg_coverage_change_classification_1990_2022.gpkg'),
            layer = "bcg_coverage_classification",
            delete_layer = TRUE)
   ###############################################################################################################
   ######################### 80% coverage threshold exceeding probability for districts############################################ 
   ################### 80% coverage for districts is the Global vaccine action plan set to be achieved by 2020 #######
   #####################################################################################################################################################
   
   estimates <- c("Prob_mean_Year_2020.tif","Prob_med_Year_2020.tif")
   raster.list <- file.path(path_output, "prevalence_merged", estimates)
   #extract raster data
   rasterls<-list()
   for (i in 1:length(raster.list)){
     rasterls[[i]]<-raster::raster(raster.list[[i]])
   }
   b <- raster::brick(rasterls)
   bcg_prev_2020 <-raster::stack(b[[1]], b[[2]])
   names(bcg_prev_2020) <- c("mean_2020", "median_2020")
   # Check raster alignment
   stopifnot(raster::compareRaster(pop_2020, bcg_prev_2020, extent=TRUE, rowcol=TRUE, crs=TRUE))
   # Initialize
   threshold <- 0.8
   pop_raster <- pop_2020
   exceed_prob_list <- list()
   
   # Loop through each layer
   for (layer_name in names(bcg_prev_2020)) {
     cat("Processing layer:", layer_name, "\n")
     
     metric_raster <- bcg_prev_2020[[layer_name]]
     stopifnot(inherits(metric_raster, "RasterLayer"))
     
     # Thresholding
     exceed_raster <- metric_raster > threshold
     exceed_raster <- exceed_raster * 1
     
     # Population Weighted Exceedance Probability
     exceed_prob <- exact_extract(
       exceed_raster,
       adm2poly,
       weights = pop_raster,
       fun = function(values, coverage_fractions, weights) {
         valid <- !is.na(values)
         if (sum(valid) == 0) return(NA_real_)
         weighted.mean(values[valid], w = weights[valid])
       },
       progress = FALSE
     )
     
     exceed_prob_list[[layer_name]] <- exceed_prob
   }
   exceed_prob_df <- as. data.frame(exceed_prob_list)
   
   # Bind results with admin polygons 
   adm2_exceed_2020 <- bind_cols(adm2poly, exceed_prob_df)
   
   # write to GeoPackage
   st_write(
     adm2_exceed_2020,
     dsn = file.path(path_output, "shp", "bcg_prob_exceed_2020.gpkg"),
     layer = "prob_exceed_80_2020",
     delete_layer = TRUE
   )
  ###to see countries with exceed provbability =1 for all districts
   countries_all_ones <- adm2_exceed_2020 %>%
        st_drop_geometry() %>%   # <- drops the geometry column
       group_by(ISO) %>%
       summarise(all_ones = all(median_2020 == 1, na.rm = TRUE)) %>%
        filter(all_ones == TRUE)
     print(countries_all_ones$ISO)
     ###to see the exceed probability of countries where 90% of districts have exceed probability higher than 95% lower than 100%
     countries_more_than_90pct_above_95_exclude_all_1 <- adm2_exceed_2020 %>%
       st_drop_geometry() %>%
       group_by(ISO) %>%
       summarise(
         total_units = n(),
         count_above_95 = sum(median_2020 > 0.95, na.rm = TRUE),
         count_equal_1 = sum(median_2020 == 1, na.rm = TRUE)
       ) %>%
       filter(
         count_above_95 > 0.9 * total_units,
         count_equal_1 < total_units   # excludes if all are exactly 1
       )
     
     print(countries_more_than_90pct_above_95_exclude_all_1$ISO)
     #### to get countries where half of the districts have exceeded the probability below 25%
     countries_half_below_025 <- adm2_exceed_2020 %>%
     st_drop_geometry() %>%
       group_by(ISO) %>%
       summarise(
         total_units = n(),
         count_below_025 = sum(median_2020 < 0.25, na.rm = TRUE)
       ) %>%
       filter(count_below_025 >= total_units / 2)
     
     print(countries_half_below_025$ISO)
     
     
###############################################################################################################
################### proportion adm2 with increased and decreased coverage from 1990 to 2022
###############################################################################################################
   library(dplyr)
  ###increased coverage
   
   # Logical vector: TRUE if absolute increase in coverage >= 25% percentage points
   success_logical <- bcg_change$delta_median >= 0.25 ## We used this threshold to see reasonable change for such long decades
   # Count successes (ignore NA values)
   successes <- sum(success_logical, na.rm = TRUE)
   # Total districts with non-missing coverage in both years
   total <- sum(
     !is.na(bcg_change$delta_median)
   )
   
   alpha <- 1 # Beta prior parameters (uninformative)
   beta <- 1# Beta prior parameters (uninformative)
   post_alpha <- successes + alpha # Posterior parameters
   post_beta <- total - successes + beta
   
   post_mean <- post_alpha / (post_alpha + post_beta) # Point estimate (posterior mean)
   lower_cri <- qbeta(0.025, post_alpha, post_beta) # 95% credible interval (quantiles of Beta posterior)
   upper_cri <- qbeta(0.975, post_alpha, post_beta)
   
   cat(sprintf("Proportion of adm2 districts with increased median coverage: %.2f%% (95%% CrI: %.2f%% - %.2f%%)\n",
               post_mean * 100, lower_cri * 100, upper_cri * 100))
   ###decreased coverage
   #posterior mean with Bayesian 95%CrI
   success_logical <- bcg_change$delta_median <=0.25
   successes <- sum(success_logical, na.rm = TRUE)
   alpha <- 1 # Beta prior parameters (uninformative)
   beta <- 1 # Beta prior parameters (uninformative)
   
   # Posterior parameters
   post_alpha <- successes + alpha
   post_beta <- total - successes + beta

   # Posterior mean and 95% credible interval
   post_mean <- post_alpha / (post_alpha + post_beta)
   lower_cri <- qbeta(0.025, post_alpha, post_beta)
   upper_cri <- qbeta(0.975, post_alpha, post_beta)
   # Output
   cat(sprintf("Proportion of adm2 districts with decreased median coverage: %.2f%% (95%% CrI: %.2f%% - %.2f%%)\n",
               post_mean * 100, lower_cri * 100, upper_cri * 100))
  
    ### actual coverage change between 1990 to 2022 across districts
   bcg_change_prop<-bcg_change_summary 
   bcg_change_prop <- bcg_change_prop %>%
     arrange(delta_median)
   summary(bcg_change_prop)
   bcg_decrease <- bcg_change[!is.na(bcg_change$delta_median) & bcg_change$delta_median < 0, ] ##No median change =0 in data
   
   # Delta median greater than 0
   bcg_increase <- bcg_change[bcg_change$delta_median > 0, ]
  
   
   
#########################################################################################################
############################### Population weighted national BCG coverage trend ##########################################################
########################################################################################################
   library(exactextractr)
   library(dplyr)
   extract_pop_by_year <- function(pop_raster, adm2poly, year) {
     pop_sum <- exact_extract(pop_raster, adm2poly, 'sum')
     adm2poly$population <- pop_sum
     adm2poly$Year <- year
     adm2poly %>%
       sf::st_drop_geometry() %>%
       dplyr::select(ISO, Name_2, Year, population)
   }
   pop_1990_df <- extract_pop_by_year(pop_1990, adm2poly, 1990)
   pop_2000_df <- extract_pop_by_year(pop_2000, adm2poly, 2000)
   pop_2010_df <- extract_pop_by_year(pop_2010, adm2poly, 2010)
   pop_2020_df <- extract_pop_by_year(pop_2020, adm2poly, 2020)
   pop_2022_df <- extract_pop_by_year(pop_2022, adm2poly, 2022)
   pop_all_years <- bind_rows(pop_1990_df, pop_2000_df, pop_2010_df, pop_2020_df, pop_2022_df)
   
   final_data_r <- final_data_r %>%
     left_join(pop_all_years, by = c("ISO", "Name_2", "Year"))
   
   # Define a helper function for weighted median with NA handling
   weighted_median_na <- function(x, w) {
     valid <- !is.na(x) & !is.na(w)
     if (sum(valid) == 0) {
       return(NA_real_)
     }
     weightedMedian(x[valid], w[valid])
   }
   
   adm0_data2 <- final_data_r %>%
     group_by(ISO, Year) %>%
     summarise(
       LUI_2.5 = weighted_median_na(LUI_2.5, population),
       mean = weighted_median_na(mean, population),
       median = weighted_median_na(median, population),
       UUI_97.5 = weighted_median_na(UUI_97.5, population),
       .groups = "drop"
     )
   
   national_trend<-ggplot(adm0_data2, aes(x = Year)) +
     geom_line(aes(y = median), color = "steelblue", size = 1) +
     geom_point(aes(y = median), color = "darkblue", size = 1.5) +
     facet_wrap(~ ISO) +
     theme_minimal() +
     labs(
       title = " Estimated BCG Coverage Over Time",
       x = "Year",
       y = "Proportion"
     )
   plot(national_trend)
   ggsave(
     filename = file.path(path_output, "national_trend.tif"),  # combine path and filename properly
     plot = national_trend,
     width = 10, height = 6, dpi = 300
   )
    plot(national_trend)
   ########################################################################################################
   #### population-weighted number of unvaccinated children at adm1
   ####################################################
   
   library(raster)
   library(dplyr)
   
   years <- c(1990, 2000, 2010, 2020, 2022) # List of years
   
   Africa.layer <- list("sp.polygons", Africa, col = "black")
   adm1poly <- list()
   mycountries <- unique(myarea$iso_a3)
   mycountries <- setdiff(unique(myarea$iso_a3), "COM")
   
    for (i in 1:length(mycountries)){
     
     adm1poly[[i]] <- readRDS(paste0(path_input, "/gadm36_",mycountries[[i]],"_1_sp.rds"))
   }
   #put polygons together
   adm1poly <- do.call(rbind,adm1poly)
  
   year_layers1 <- list()  # Initialize list to store results
   
   for (year in years) {
     # Get unvaccinated and population rasters
     unvaccinated_layer <- get(paste0("unvaccinated_count_", year))
     population_layer   <- get(paste0("pop_", year))  # Assumes raster names like pop_1990, etc.
     # Extract total unvaccinated counts per admin1
     unvaccinated_data <- raster::extract(unvaccinated_layer, adm1poly, fun = sum, na.rm = TRUE, sp = TRUE)
     # Extract total population per admin1
     population_data <- raster::extract(population_layer, adm1poly, fun = sum, na.rm = TRUE, sp = TRUE)
     pop_col <- grep("^ihme_corrected_worldpop", names(population_data@data), value = TRUE)
     names(population_data@data)[names(population_data@data) == pop_col] <- "Population"
     
     # Merge and calculate
     df <- unvaccinated_data@data
     df$Population <- population_data@data$Population
     df$Year <- year
     df$Rate_LUI_2.5  <- (df$LUI_2.5   / df$Population) * 100000
     df$Rate_Mean     <- (df$mean      / df$Population) * 100000
     df$Rate_Median   <- (df$median    / df$Population) * 100000
     df$Rate_UUI_97.5 <- (df$UUI_97.5  / df$Population) * 100000
     
     
     # Rename selected columns for compatibility
     df <- df %>%
       dplyr::rename(
         LUI         = LUI_2.5,
         UUI         = UUI_97.5,
         rated_mean  = Rate_Mean,
         rated_median = Rate_Median,
         rated_LUI   = Rate_LUI_2.5,
         rated_UUI   = Rate_UUI_97.5
       )
     # Select columns
     df <- df %>%
       dplyr::select(GID_0, NAME_0, GID_1, NAME_1, Year,
                     LUI, mean, median, UUI, Population,
                     rated_LUI, rated_mean, rated_median, rated_UUI)
     
     # Assign back to the spatial object
     unvaccinated_data@data <- df
     
     # Store result
     year_layers1[[as.character(year)]] <- unvaccinated_data
   }
   
   
   ### export shapefile for each year
  
   for (y in 1:length(year_layers1)) {
     year_data <- year_layers1[[y]]
     year <- year_data$Year[1]  # Access 'Year' column directly
     year_data_sf <- st_as_sf(year_data) # Convert to sf object
     gpkg_name <- file.path(path_output, 'shp', paste0('new_weighted_unvac_count_adm1_', year, '.gpkg'))
     # Save as GeoPackage
     st_write(year_data_sf, dsn = gpkg_name, layer = paste0('unvac_count_adm1_', year), delete_layer = TRUE)
   }
  ################################################################################################################################### 
   ### compute model validation metrics ################################################################################################# #   
  ###################################################################################
    ###continental model
   idx1 <- inla.stack.index(stack.est1, tag = "est")$data
   predicted_mean     <- mod.pred1$summary.fitted.values[idx1, "mean"]
   predicted_sd       <- mod.pred1$summary.fitted.values[idx1, "sd"]
   predicted_ci2.5    <- mod.pred1$summary.fitted.values[idx1, "0.025quant"]
   predicted_ci97.5   <- mod.pred1$summary.fitted.values[idx1, "0.975quant"]
   
   # Create a dataframe with the predictions
   bcg.df1 <- data.frame(
     vaccinated_c= data_1$vaccinated_c,
     total_c=data_1$total_c,
     year = data_1$year,
     x = data_1$longitude,
     y= data_1$latitude,
     predicted_mean,
     predicted_sd,
     predicted_ci2.5,
     predicted_ci97.5
   )
   bcg.df1 <- bcg.df1 %>%
     mutate(p = vaccinated_c / total_c)
   
   bcg.df1 <- data.frame(
     vaccinated_c = data_1$vaccinated_c,
     total_c      = data_1$total_c,
     year         = data_1$year,
     x            = data_1$longitude,
     y            = data_1$latitude,
     predicted_mean,
     predicted_sd,
     predicted_ci2.5,
     predicted_ci97.5
   ) 
   
   # Define RMSE function
   RMSE <- function(actual, predicted) {
     sqrt(mean((actual - predicted)^2, na.rm = TRUE))
   }
   
   # Calculate model metrics
   model_metrics1 <- data.frame(
     R_squared       = cor(bcg.df1$p, bcg.df1$predicted_mean, use = "complete.obs")^2,
     Mean_Error      = mean(bcg.df1$predicted_mean - bcg.df1$p, na.rm = TRUE),
     RMSE            = RMSE(bcg.df1$p, bcg.df1$predicted_mean),
     MAE             = mean(abs(bcg.df1$predicted_mean - bcg.df1$p), na.rm = TRUE),
     Mean_NegLogCPO  = -mean(log(mod.pred1$cpo$cpo), na.rm = TRUE)
   )
   
   write.csv(model_metrics1,paste0(path_output,"/csv/model1_metrics.csv")) 
   #################################################################################NEW#####################################################
   #####for Madagascar
   ########################################################################
   
   idx2 <- inla.stack.index(stack.est2, tag = "est")$data
   predicted_mean     <- mod.pred1$summary.fitted.values[idx2, "mean"]
   predicted_sd       <- mod.pred1$summary.fitted.values[idx2, "sd"]
   predicted_ci2.5    <- mod.pred1$summary.fitted.values[idx2, "0.025quant"]
   predicted_ci97.5   <- mod.pred1$summary.fitted.values[idx2, "0.975quant"]
   
   # Create a dataframe with the predictions
   bcg.df2 <- data.frame(
     vaccinated_c= data_2$vaccinated_c,
     total_c=data_2$total_c,
     year = data_2$year,
     x = data_2$longitude,
     y= data_2$latitude,
     predicted_mean,
     predicted_sd,
     predicted_ci2.5,
     predicted_ci97.5
   )
   bcg.df2 <- bcg.df2 %>%
     mutate(p = vaccinated_c / total_c)
   
   bcg.df2 <- data.frame(
     vaccinated_c = data_2$vaccinated_c,
     total_c      = data_2$total_c,
     year         = data_2$year,
     x            = data_2$longitude,
     y            = data_2$latitude,
     predicted_mean,
     predicted_sd,
     predicted_ci2.5,
     predicted_ci97.5
   ) 
   
   # Define RMSE function
   RMSE <- function(actual, predicted) {
     sqrt(mean((actual - predicted)^2, na.rm = TRUE))
   }
   
   # Calculate model metrics
   model_metrics2 <- data.frame(
     R_squared       = cor(bcg.df2$p, bcg.df2$predicted_mean, use = "complete.obs")^2,
     Mean_Error      = mean(bcg.df2$predicted_mean - bcg.df2$p, na.rm = TRUE),
     RMSE            = RMSE(bcg.df2$p, bcg.df2$predicted_mean),
     MAE             = mean(abs(bcg.df2$predicted_mean - bcg.df2$p), na.rm = TRUE),
     Mean_NegLogCPO  = -mean(log(mod.pred2$cpo$cpo), na.rm = TRUE)
   )
   
   write.csv(model_metrics,paste0(path_output,"/csv/model_metrics.csv")) 
   
   ## validation plots
   ########################################################################
   library(dplyr)
   
  
   bcg.df_all <- dplyr::bind_rows(bcg.df1, bcg.df2)
 
   
   ggplot(bcg.df_all, aes(x = predicted_mean, y = p)) +  # Use the original data frame
     geom_point(alpha = 0.5) +  # Use alpha for transparency to see overlapping points
     geom_abline(slope = 1, intercept = 0, color = "red") +  # Reference line for perfect prediction
     labs(x = " Predicted Probability", 
          y = "Observed Probability") +  # Adding a title
     facet_wrap(~ year, ncol = 4) +  # Arrange plots in four columns
     scale_x_continuous(breaks = seq(0, 1, by = 0.25), 
                        labels = seq(0, 1, by = 0.25),  # Format labels as percentages
                        limits = c(0, 1)) +  # Set limits for x-axis
     scale_y_continuous(breaks = seq(0, 1, by = 0.25), 
                        labels = seq(0, 1, by = 0.25),  # Format labels as percentages
                        limits = c(0, 1)) +  # Set limits for y-axis
     theme_minimal() +
     theme(
       plot.title = element_text(size = 16, hjust = 0.5),  # Center title
       axis.title.x = element_text(size = 14),  # X-axis label font size
       axis.title.y = element_text(size = 14),  # Y-axis label font size
       axis.text = element_text(size = 12),      # Axis tick labels font size
       strip.text = element_text(size = 12)       # Facet label font size
     )
   ######################  
   
