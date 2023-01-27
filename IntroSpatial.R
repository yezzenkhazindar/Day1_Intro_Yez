packageLoad <-
  function(x) {
    for (i in 1:length(x)) {
      if (!x[i] %in% installed.packages()) {
        install.packages(x[i])
      }
      library(x[i], character.only = TRUE)
    }
  }

# create a string of package names
packages <- c('tidyverse',
              'palmerpenguins',
              'sf',
              'terra',
              'tmap',
              'rmarkdown',
              'tigris',
              'elevatr',
              'rgdal')

packageLoad(packages)
setwd("/Users/yezzenkhazindar/Downloads/Courses/ESS523_Environmental_Data_Science_Applications/Day1_Intro_Yez")
source("Setup.R")

# download county shapefile for the state of Colorado
counties <- counties(state = "CO")

roads <- roads(state = "CO", county = "Larimer")

# Set tmap mode to interactive
tmap_mode("view")

qtm(counties) +
  qtm(roads)

tm_shape(counties)+
  tm_polygons()

# look at the class of counties
class(counties)

# Filter the dataset to show Poudre Canyon Highway
poudre_hwy <- roads %>% 
  filter(FULLNAME == "Poudre Canyon Hwy")

qtm(poudre_hwy)

# point data
poudre_points <- data.frame(name = c("Mishawaka", "Rustic", "Blue Lake Trailhead"),
                            long = c(-105.35634, -105.58159, -105.85563),
                            lat = c(40.68752, 40.69687, 40.57960))
# convert to spatial
poudre_points_sf <- st_as_sf(poudre_points, coords = c("long", "lat"), crs = 4326)

qtm(poudre_hwy)+
  qtm(poudre_points_sf)



# raster data
elevation <- get_elev_raster(counties, z = 7)

qtm(elevation)

tm_shape(elevation)+
  tm_raster(style = "cont", title = "Elevation (m)")

# the terra package
## we want to convert the raster layer to a terra raster object

elevation <- rast(elevation)

names(elevation) <- "Elevation"

elevation

# we see on the map that the boundary of the raster extends outside the bounds of CO
# so we can crop it to the boundary

# see the CRS in the header metadata
counties

# check projections
st_crs(counties) # checks the projection of sf objects

crs(counties) == crs(elevation) # crs() function checks projections with terra of vectors and rasters

poudre_points_prj <- st_transform(poudre_points_sf, st_crs(counties))

#Now check that they match
st_crs(poudre_points_prj) == st_crs(counties)

# project elevation layer
elevation_prj <- project(elevation, counties) # error bc we need to use ::
elevation_prj <- terra::project(elevation, counties)

#crop elevation to counties extent
elevation_crop <- crop(elevation, ext(counties))

qtm(elevation_crop)

# final map with all the spatial data we created

tm_shape(elevation, bbox = st_bbox(poudre_hwy))+
  tm_raster(style = "cont", title = "Elevation (m)")+
  tm_shape(poudre_hwy)+
  tm_lines()+
  tm_shape(poudre_points_prj)+
  tm_dots(size = 0.2)

# read and write spatial data

# To save vector data with sf, use write_sf()
write_sf(poudre_hwy, "data/poudre_hwy.shp")

write_sf(poudre_points_prj, "data/poudre_points.shp")

# To save raster data with terra use writeRaster()

writeRaster(elevation_crop, "data/elevation_larimer.tif")

# Since the poudre_hwy and poudre_points_prj were objects you 
# created in this session, to avoid the need to recreate them 
# you can save them to an .RData file:

save(poudre_hwy, poudre_points_prj, file = "data/spatial_objects.RData")

rm(poudre_hwy, poudre_points_prj)

load("data/spatial_objects.RData")

# Note that terra objects don’t properly save to .RData files, 
# but there is a work around if you save a single terra object 
# as an .RDS file. Here is that workflow, there is just a second 
# step to ‘unpack’ the loaded .RDS object with rast().

saveRDS(elevation_crop, "data/elevation_crop.RDS")

readRDS("data/elevation_crop.RDS") %>% rast()

# To read in shapefiles, you use read_sf()

read_sf("data/poudre_hwy.shp")

rast("data/elevation_larimer.tif")

# EXERCISES:

# 1. Filter out the counties data set to only include Larimer, Denver, and Pueblo counties.

NOCO <- counties %>% 
  filter(NAME %in% c('Larimer', 'Denver', 'Pueblo'))
NOCO

# 2. Make a map of the counties data colored by county area. Make a second map of counties colored by their total area of water.
?tm_polygons
tm_shape(NOCO) + 
  tm_polygons(col = 'NAME',)

tm_shape(NOCO) + 
  tm_polygons(col = 'AWATER',)

# 3. Make a barplot comparing the elevation of your 3 points in the Poudre Canyon (note: explore the extract() function in the terra package).

#crop elevation to counties extent
elevation_crop_poudre <- crop(elevation, ext(poudre_hwy))

writeRaster(elevation_crop_poudre, "data/elevation_poudre.tif", overwrite=TRUE)

poudre_point_elevation <- extract(x= elevation_crop,
                                  y= poudre_points_prj,
                                  )

poudre_point_elevation

poudre_point_elevation %>%
  ggplot() +
  geom_col(mapping = aes(x=ID, y=Elevation, fill=ID))

# 4. Why are there 4 features in our Poudre Canyon Highway variable instead of 1?
















