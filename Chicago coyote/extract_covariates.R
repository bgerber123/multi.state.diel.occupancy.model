

library(raster)
library(sf)
library(dplyr)
library(mapview)
library(cli)


source("./Chicago coyote/spatial_utilities.R")

# Need a table of site locations, which we can pull from the coyote data
coyote <- read.csv(
  "./Chicago coyote/data/day_night_detections.csv"
)

# reduce down to unique sites
site_coords <- coyote[,c("Site", "lat", "Long")]
site_coords <- site_coords[!duplicated(site_coords),]

# Create spatial points
# You must put the correct CRS for respective city
sites <- sf::st_as_sf(
  site_coords,
  coords = c("Long", "lat"),
  crs = 4326
)
colnames(sites)[1] <- c("LocationName")

# Visually inspect points to ensure the projection and location is correct
mapeview::mapview(sites)

# We will use the High-res Landcover for NE Illinois for this example
#  Data can be downloaded from:
#  browseURL("https://datahub.cmap.illinois.gov/dataset/high-resolution-land-cover-ne-illinois-and-nw-indiana-2010")
#  Load iLULC map
# REPLACE FILE PATH WITH LOCAL FILE PATH

my_raster_path <- 
  "D:/GIS/cmap/landcover_2010_chicagoregion.img"

# Check if the file exists, if not provide a warning.
check_path(my_raster_path)

# read it in
my_map <- raster::raster(my_raster_path)

#  For this example we will extract the proportion canopy cover (lulc class 1)
#    and create our own 'impervious cover' value, which is the sum of multiple
#    lulc classes.
lulc_prop <- extract_raster_prop(my_points = sites,
                                 my_buffer = 1000,
                                 my_raster_data = my_map,
                                 lulc_cats = list("tree" = 1,
                                                  "imperv" = 5:7
                                 ),
                                 point_names = sites$LocationName
)

write.csv(
  lulc_prop,
  "./Chicago coyote/data/tree_imperv.csv",
  row.names = FALSE
)

# Now we need to extract the housing density data



pop_data <- sf::st_read(
  "D:/GIS/housing_density", 
  layer = "il_blk10_Census_change_1990_2010_PLA2"
)

pop_data <- sf::st_make_valid(pop_data)

# Run function to calculate housing units, housing density, population 
#  and population density.  For this example we extract population data 
#  within a 1km radius buffer.
population_data <- extract_polygon(
  my_points = sites,
  my_buffer = 1000,
  my_shp = pop_data,
  layers = c("HU10")
)


# join all the data together
covs <- dplyr::inner_join(
  lulc_prop,
  population_data,
  by = "LocationName"
)

covs$HU10 <- as.numeric(covs$HU10)

# Now make the urbanization covariate

urb <- prcomp(
  covs[,-1],
  scale. = TRUE
)

covs$urb <- urb$x[,1]

write.csv(
  covs,
  "./Chicago coyote/data/chicago_covars.csv",
  row.names = FALSE
)
