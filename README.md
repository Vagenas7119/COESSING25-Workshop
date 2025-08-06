# Ecological & Spatial Data Analysis Lab - COESSING25  

**Author**: Georgios Vagenas (MNCN, CSIC) | Predoctoral Investigador & PhD(c) at the National Spanish Research Council (CSIC) 

**Project**: The Coastal Ocean Environment Summer School In Nigeria and Ghana (COESSING) - COESSING 2025

**Location**: University of Accra, Ghana (West Africa), August 2025 - Ecological Data Analysis Workshop


![Screenshot 2025-08-05 190504](https://github.com/user-attachments/assets/b484e54e-13d4-4e17-9b9c-51cdd53c0a0b)

## 📌 Overview  
This repository contains an R workflow for analyzing ecological data, estimating biomass using the BioTIME and FishBase database, and visualizing spatiotemporal patterns in species richness across the globe. The pipeline processes occurrence records, computes ecological metrics, and generates animated maps of biodiversity changes.

Link to download the layers: https://saco.csic.es/s/ZdE9mHqrYBTfs5F

## 0. Setup Working Environment
```r
# Set working directory
setwd("C:/Users/XXX/Desktop/COESSING25_Lab_Files/")

# Install required packages
install.packages(c("terra", "mapview", "dplyr", "ggplot2", "rfishbase", 
                  "ncdf4", "gganimate", "rnaturalearth", "sf", "tidyr", "vegan"))

# Load libraries
library(terra)
library(mapview)
library(dplyr)
library(ggplot2)
library(rfishbase)
library(ncdf4)
library(gganimate)
library(rnaturalearth)
library(sf)
library(tidyr)
library(vegan)
```

## 1. Load and prepare the data
```r
md <- read.csv("raw_data_292.csv") # Australian dataset

# Create coordinate points
md_coords <- data.frame(LATITUDE = md$LATITUDE, LONGITUDE = md$LONGITUDE)
points <- vect(md_coords, geom = c("LONGITUDE", "LATITUDE"), crs = "EPSG:4326")

# Load oceanographic data
nc_file_pp <- rast("layers/phyc_baseline_2000_2020_depthsurf_7d39_02af_cdbd_U1751549420494.nc")
```
## 2. Descriptive statistics
```r
# Calculate species abundance
species_totals <- aggregate(ABUNDANCE ~ valid_name, data = md, FUN = sum)
all_species <- species_totals[order(-species_totals$ABUNDANCE), ]

# Generate summary statistics
results <- lapply(all_species$valid_name, function(species) {
  subset_data <- md[md$valid_name == species, ]
  data.frame(
    Species = species,
    Total_Abundance = sum(subset_data$ABUNDANCE),
    Mean = mean(subset_data$ABUNDANCE),
    Median = median(subset_data$ABUNDANCE),
    SD = sd(subset_data$ABUNDANCE),
    Min = min(subset_data$ABUNDANCE),
    Max = max(subset_data$ABUNDANCE),
    N_Samples = nrow(subset_data)
  )
})
all_stats <- do.call(rbind, results)

```
## 3. Biomass estimation
```r
# Get FishBase parameters
species_list <- all_stats$Species
weight_data <- popchar(species_list)

# Calculate length-weight relationships
lw_params <- length_weight(species_list) %>%
  select(Species, a, b) %>%
  mutate(a = as.numeric(a), b = as.numeric(b)) %>%
  group_by(Species) %>%
  summarise(a_mean = mean(a, na.rm = TRUE), b_mean = mean(b, na.rm = TRUE))

# Estimate biomass
weight_estimates <- weight_data %>%
  left_join(lw_params, by = "Species") %>%
  mutate(Estimated_Wmax_kg = (a_mean * (as.numeric(Lmax)^b_mean))/1000)
```

## 4. Temporal analysis
```r
# Annual landings time series
annual_landings <- md %>%
  mutate(biomass_kg = ABUNDANCE * weight_estimates$Estimated_Wmax_kg[match(valid_name, weight_estimates$Species)]) %>%
  group_by(YEAR) %>%
  summarise(total_landings_kg = sum(biomass_kg, na.rm=TRUE))

# Plot landings
landings_plot <- ggplot(annual_landings, aes(x=YEAR, y=total_landings_kg/1e6)) +
  geom_line(linewidth=1) +
  labs(title="Total Estimated Landings by Year", y="Landings (million kg)")
  ```
![total](https://github.com/user-attachments/assets/e15144fd-3a7a-46ab-a772-d55446f717e2)

## 5. Spatial analysis
```r
# Create species richness grid
annual_richness <- md %>%
  mutate(grid_x = round(LONGITUDE/0.1)*0.1,
         grid_y = round(LATITUDE/0.1)*0.1) %>%
  group_by(grid_x, grid_y, YEAR) %>%
  summarize(n_species = n_distinct(valid_name))

# Generate animated plot
animated_plot <- ggplot() +
  geom_tile(data = annual_richness, 
           aes(x = grid_x, y = grid_y, fill = n_species)) +
  transition_time(YEAR) +
  labs(title = 'Year: {frame_time}')

# Save animation
animate(animated_plot, width = 1000, height = 800, fps = 3, res = 150)
anim_save("species_richness_animation.gif")
```

![Change](https://github.com/user-attachments/assets/c448934d-6cd1-42de-b3c3-e7bd365a285d)

## 6. Export results
```r
# Save raster stack
r_stack <- rast(lapply(split(annual_richness, annual_richness$YEAR), function(x) {
  r <- rast(x[, c("grid_x", "grid_y", "n_species")], crs = "EPSG:4326")
}))
writeRaster(r_stack, "annual_richness_raster_stack.tif")

# Save summary statistics
write.csv(all_stats, "species_summary_stats.csv")
write.csv(weight_estimates, "biomass_estimates.csv")
```
![species_richness_animation](https://github.com/user-attachments/assets/96d838bc-6af9-408b-8064-281bf38968f4)
