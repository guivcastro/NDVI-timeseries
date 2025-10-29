# Set working directory
setwd("C:/Users/guilh/OneDrive/Documentos/PhD/Data/Climate/Sussex")

# Load required libraries
library(terra)
library(sf)
library(ncdf4)
library(ggplot2)
library(dplyr)
library(tidyr)
library(Kendall)
library(trend)
library(stringr)

# Load Sussex shapefile
sussex_shp <- st_read("Sussex_boundaries.shp")
sussex_vect <- vect(sussex_shp)

# ------------------------------------- Temperature ----------------------------------------

# Function to read and clip .files 

read_nc_folder_clipped <- function(pattern, varname, path = ".") {
 files <- list.files(path, pattern = pattern, full.names = TRUE)
 
 if(length(files) == 0) stop(paste("No files found for pattern:", pattern))
 
 data_list <- lapply(files, function(f) {
  # Load raster
  r <- rast(f)
  
  # Select the variable
  if(!(varname %in% names(r))) stop(paste("Variable", varname, "not found in", f))
  r <- r[[varname]]
  
  # Clip to Sussex
  r_sussex <- crop(r, sussex_vect)
  r_sussex <- mask(r_sussex, sussex_vect)
  
  # Annual mean over Sussex
  ann_mean <- global(r_sussex, "mean", na.rm = TRUE)[1,1]
  
  # Extract year from filename
  year <- as.numeric(str_extract(basename(f), "\\d{4}"))
  
  data.frame(year = year, ann = ann_mean)
 })
 
 df <- bind_rows(data_list) %>% arrange(year)
 return(df)
}

mean_temp <- read_nc_folder_clipped("tas_hadukgrid_uk_1km_ann", "tas")
min_temp  <- read_nc_folder_clipped("tasmin_hadukgrid_uk_1km_ann", "tasmin")
max_temp  <- read_nc_folder_clipped("tasmax_hadukgrid_uk_1km_ann", "tasmax")

# Filter data to 1995-2024

mean_temp_filtered <- mean_temp %>% filter(year >= 1995 & year <= 2024)
max_temp_filtered  <- max_temp %>% filter(year >= 1995 & year <= 2024)
min_temp_filtered  <- min_temp %>% filter(year >= 1995 & year <= 2024)

# Mann-Kendall tests
mann_kendall_mean <- Kendall(mean_temp_filtered$year, mean_temp_filtered$ann)
mann_kendall_max  <- Kendall(max_temp_filtered$year,  max_temp_filtered$ann)
mann_kendall_min  <- Kendall(min_temp_filtered$year,  min_temp_filtered$ann)

tau_mean  <- round(mann_kendall_mean$tau, 3)
p_value_mean <- round(mann_kendall_mean$sl, 3)
tau_max   <- round(mann_kendall_max$tau, 3)
p_value_max <- round(mann_kendall_max$sl, 3)
tau_min   <- round(mann_kendall_min$tau, 3)
p_value_min <- round(mann_kendall_min$sl, 3)

temp_ribbon_df <- mean_temp_filtered %>%
 select(year, ann) %>% rename(mean_value = ann) %>%
 left_join(max_temp_filtered %>% select(year, ann) %>% rename(max_value = ann), by = "year") %>%
 left_join(min_temp_filtered %>% select(year, ann) %>% rename(min_value = ann), by = "year")

mean_temp_sd <- mean_temp_filtered %>% summarise(sd_value = sd(ann, na.rm = TRUE))

# Plot
ggplot() +
 geom_ribbon(data = temp_ribbon_df,
             aes(x = year, ymin = min_value, ymax = max_value,
                 fill = "Temperature range (minimum and maximum)"),
             alpha = 0.5, show.legend = TRUE) +
 
 geom_line(data = mean_temp_filtered,
           aes(x = year, y = ann, color = "Mean temperature (°C)"),
           size = 0.8, show.legend = TRUE) +
 
 geom_point(data = mean_temp_filtered,
            aes(x = year, y = ann, color = "Mean temperature (°C)"),
            size = 2, shape = 16, color = "#8B2323") +
 
 geom_ribbon(data = mean_temp_filtered,
             aes(x = year, ymin = ann - mean_temp_sd$sd_value,
                 ymax = ann + mean_temp_sd$sd_value,
                 fill = "Standard deviation"),
             alpha = 0.2, show.legend = TRUE) +
 
 labs(x = "Year", y = "Temperature (°C)") +
 theme_minimal() +
 theme(
  text = element_text(size = 30),
  legend.position = "bottom",
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(),
  legend.key = element_blank(),
  legend.box.spacing = unit(1, "cm"),
  legend.spacing.x = unit(1, "cm"),
  axis.line = element_line(colour = "black", linewidth = 0.5),
  axis.ticks = element_line(colour = "black", linewidth = 0.7)
 ) +
 scale_color_manual(values = c("Mean temperature (°C)" = "#8B0000"),
                    name = NULL) +
 scale_fill_manual(values = c("Temperature range (minimum and maximum)" = "#CD9B9B",
                              "Standard deviation" = "#8B475D"),
                   name = NULL) +
 guides(fill = guide_legend(override.aes = list(color = NA), keyheight = 1.5),
        color = guide_legend(override.aes = list(fill = NA), keyheight = 1.5)) +
 scale_x_continuous(
  breaks = seq(min(mean_temp_filtered$year), max(mean_temp_filtered$year), by = 5),
  limits = c(min(mean_temp_filtered$year), max(mean_temp_filtered$year) + 5)
 ) +
 scale_y_continuous(expand = expansion(add = c(0, 1))) +
 annotate("text",
          x = max(mean_temp_filtered$year) + 0.5, y = max(mean_temp_filtered$ann),
          label = paste("\nTau:", tau_mean, "\np-value:", p_value_mean),
          hjust = 0, size = 8, color = "#404040") +
 annotate("text",
          x = max(mean_temp_filtered$year) + 0.5, y = max(max_temp_filtered$ann) - 1,
          label = paste("\nTau:", tau_max, "\np-value:", p_value_max),
          hjust = 0, size = 8, color = "#404040") +
 annotate("text",
          x = max(mean_temp_filtered$year) + 0.5, y = max(min_temp_filtered$ann),
          label = paste("\nTau:", tau_min, "\np-value:", p_value_min),
          hjust = 0, size = 8, color = "#404040")

# Theil–Sen slopes
ts_mean <- sens.slope(mean_temp_filtered$ann)
ts_max  <- sens.slope(max_temp_filtered$ann)
ts_min  <- sens.slope(min_temp_filtered$ann)

ts_results <- data.frame(
 Variable = c("Mean Temperature", "Max Temperature", "Min Temperature"),
 TS_slope = c(round(ts_mean$estimates,5), round(ts_max$estimates,5), round(ts_min$estimates,5)),
 CI_lower = c(round(ts_mean$conf.int[1],5), round(ts_max$conf.int[1],5), round(ts_min$conf.int[1],5)),
 CI_upper = c(round(ts_mean$conf.int[2],5), round(ts_max$conf.int[2],5), round(ts_min$conf.int[2],5))
)
print(ts_results)

# Plot mean temperature map
# Stack all tas rasters
tas_files <- list.files(pattern = "tas_hadukgrid_uk_1km_ann", full.names = TRUE)
tas_stack <- rast(lapply(tas_files, rast))

# Clip and mask each layer
tas_stack_clipped <- lapply(1:nlyr(tas_stack), function(i) {
 r_sussex <- crop(tas_stack[[i]], sussex_vect)
 mask(r_sussex, sussex_vect)
})

# Calculate mean across all years (1995-2024)
tas_mean_sussex <- rast(tas_stack_clipped)
tas_mean_sussex <- mean(tas_mean_sussex, na.rm = TRUE)

# Plot mean map
plot(tas_mean_sussex, main = "Average Temperature (1995-2024) - Sussex",
     col = terrain.colors(20))

# ------------------------------------- Rainfall ----------------------------------------

# Function to read and clip rainfall .nc files 

read_nc_rainfall <- function(pattern, varname, path = ".") {
 files <- list.files(path, pattern = pattern, full.names = TRUE)
 if(length(files) == 0) stop(paste("No files found for pattern:", pattern))
 
 data_list <- lapply(files, function(f) {
  r <- rast(f)
  
  # Select variable
  if(!(varname %in% names(r))) stop(paste("Variable", varname, "not found in", f))
  r <- r[[varname]]
  
  # Clip to Sussex
  r_sussex <- crop(r, sussex_vect)
  r_sussex <- mask(r_sussex, sussex_vect)
  
  # Annual mean over Sussex
  ann_mean <- global(r_sussex, "mean", na.rm = TRUE)[1,1]
  
  # Extract year from filename
  year <- as.numeric(str_extract(basename(f), "\\d{4}"))
  
  data.frame(year = year, ann = ann_mean)
 })
 
 df <- bind_rows(data_list) %>% arrange(year)
 return(df)
}

# Read rainfall data
rainfall <- read_nc_rainfall("rainfall_hadukgrid_uk_1km_ann", "rainfall")

# Filter years 1995 to 2024
rainfall_filtered <- rainfall %>% filter(year >= 1995 & year <= 2024)

# Mann-Kendall test
mk_rain <- Kendall(rainfall_filtered$year, rainfall_filtered$ann)
tau_rain <- round(mk_rain$tau, 3)
p_value_rain <- round(mk_rain$sl, 3)

# Theil-Sen slope
ts_rain <- sens.slope(rainfall_filtered$ann)
ts_results_rain <- data.frame(
 Variable = "Rainfall",
 TS_slope = round(ts_rain$estimates, 5),
 CI_lower = round(ts_rain$conf.int[1], 5),
 CI_upper = round(ts_rain$conf.int[2], 5)
)
print(ts_results_rain)

# Compute standard deviation bounds for rainfall
rainfall_filtered <- rainfall_filtered %>%
 mutate(lower_bound = ann - sd(ann, na.rm = TRUE),
        upper_bound = ann + sd(ann, na.rm = TRUE)) %>%
 rename(annual_value = ann)

# Assign Mann-Kendall results
tau_value <- tau_rain
p_value <- p_value_rain

# Plot rainfall trend
ggplot() +
 # Add ribbon for standard deviation
 geom_ribbon(data = rainfall_filtered, 
             aes(x = year, ymin = lower_bound, ymax = upper_bound, fill = "Standard deviation"), 
             alpha = 0.5, show.legend = TRUE) +
 
 # Line for annual rainfall
 geom_line(data = rainfall_filtered, 
           aes(x = year, y = annual_value, color = "Mean rainfall (mm)"), 
           size = 1.2, show.legend = TRUE) +
 
 # Points for each year
 geom_point(data = rainfall_filtered, 
            aes(x = year, y = annual_value, color = "Mean rainfall (mm)"), 
            size = 2, shape = 16, color = "darkblue") +
 
 # Labels and theme
 labs(x = "Year", y = "Rainfall (mm)") +
 theme_minimal() +
 theme(text = element_text(size = 30),
       legend.position = "bottom",
       panel.grid.major = element_blank(),
       panel.grid.minor = element_blank(),
       legend.key = element_blank(),
       axis.line = element_line(colour = "black", linewidth = 0.5),
       axis.ticks = element_line(colour = "black", linewidth = 0.7)) +
 
 # Manual color scales
 scale_color_manual(values = c("Mean rainfall (mm)" = "#104E8B"), name = NULL) +
 scale_fill_manual(values = c("Standard deviation" = "lightblue"), name = NULL) +
 
 # Legends
 guides(fill = guide_legend(override.aes = list(color = NA), title = NULL), 
        color = guide_legend(override.aes = list(fill = NA))) +
 
 # X-axis breaks
 scale_x_continuous(
  breaks = seq(min(rainfall_filtered$year), max(rainfall_filtered$year), by = 5),
  limits = c(min(rainfall_filtered$year), max(rainfall_filtered$year) + 5)
 ) +
 
 # Add Mann-Kendall tau and p-value on right-hand side
 annotate("text", 
          x = max(rainfall_filtered$year) + 1, 
          y = max(rainfall_filtered$annual_value, na.rm = TRUE), 
          label = paste("\nTau:", round(tau_value, 2), "\np-value:", format(p_value, digits = 3)), 
          hjust = 0, size = 8, color = "#404040")

# Map with average rainfall
# Stack all rainfall rasters
rain_files <- list.files(pattern = "rainfall_hadukgrid_uk_1km_ann", full.names = TRUE)
rain_stack <- rast(lapply(rain_files, rast))

# Clip and mask each layer
rain_stack_clipped <- lapply(1:nlyr(rain_stack), function(i) {
 r_sussex <- crop(rain_stack[[i]], sussex_vect)
 mask(r_sussex, sussex_vect)
})

# Compute mean rainfall over 1995-2024
rain_mean_sussex <- rast(rain_stack_clipped)
rain_mean_sussex <- mean(rain_mean_sussex, na.rm = TRUE)

# Plot mean rainfall map
blue_palette <- colorRampPalette(c("#deebf7", "#3182bd"))(20)  # light to dark blue
plot(rain_mean_sussex, main = "Average Rainfall (1995-2024) - Sussex",
     col = blue_palette)
