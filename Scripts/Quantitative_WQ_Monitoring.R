#################################################
### Detecting and Quantifying Water Quality   ###
### Using Satellite Remote Sensing Algorithms ###
#################################################

#### Quantitative WQ Monitoring ####
# Author: Richard A. Johansen
# Date: March 25th 2021
# Source: https://github.com/ranghetti/sen2r
citation('sen2r')
citation('waterquality')

#Required R Packages
#Packages must be installed using install.packages("Package Name") or 
# devtools::install_github, if this is the first time you are using these packages.
library(tidyverse)
library(raster)
library(waterquality)
library(sf)
library(sp)
library(gdalUtils)
library(magrittr)
library(rgdal)
library(caret)

#install.packages("sen2r")
#install.packages(c("leafpm", "mapedit", "shiny", "shinyFiles", "shinydashboard", "shinyjs", "shinyWidgets",'geojsonlint'))
library(sen2r)

#### Use Sen2R for Sen2Cor Processing ------------------------------------------

#Check Dependencies
check_sen2r_deps()

#Sen2R can also be used to locally convert L1C Data to L2A Data
#Get a cup of coffee (~20-25 min run time for full scene)
sen2cor('Milford_20180709/S2B_MSIL1C_20180709T171159_N0206_R112_T14SPJ_20180709T204427.SAFE',
        outdir = "Milford_20180709")

#### Preprocess  Clipped & Masked Reflectance Imagery ####

#Image Directory Folder
# Bands 1 from 60m and Band 8 from 10m must be placed in the 20 m folder.
#Might want to create new folder for only bands 1-8A
Image_Directory = "Milford_Quantitative/S2B_MSIL2A_20180709T171159_N9999_R112_T14SPJ_20210319T155826.SAFE/GRANULE/L2A_T14SPJ_A007002_20180709T171818/IMG_DATA/R20m_WQ"

# Extacts all raster files from Image Directory with extention .jp2
Rasters = dir(Image_Directory, pattern = "*.jp2$", full.names = TRUE)

#Import Shapefile of Area of Interest
AOI = st_read("Milford_Quantitative/Milford_Reservoir.shp")

#Reproject AOI if needed (Example - UTM Zone 14)
Projection = "+proj=utm +zone=14 +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0"
AOI = st_transform(AOI, Projection)
AOI <- st_zm(AOI)
AOI = as(AOI, "Spatial")

# Create Template for Final Raster Stack
#For Sentinel-2 we want our final output to be in 20m so use Band 5
raster_B5 = raster(Rasters[[5]])
raster_template = raster(raster_B5)

# Resample, Crop, & Stack All Images
# **Takes a few minutes
raster_stack = Rasters %>% 
  lapply(raster) %>% 
  lapply(resample, raster_template) %>%
  lapply(crop, AOI) %>%
  stack()

# Mask cropped image for further reduction of the stacked image
Milford_L2A_Masked <- raster::mask(raster_stack,AOI)

#Plot Masked RGB Imagery
raster::plotRGB(Milford_L2A_Masked,
                r=4,
                g=3,
                b= 2,
                stretch='lin')

#Save Final Stacked Image as Tiff
writeRaster(x = Milford_L2A_Masked,
            filename= "Milford_Quantitative/Milford_L2A_Masked_20180709.tif", # save as a tif
            datatype ="FLT4S", # save as a float 4 significant digits
            overwrite = FALSE) #Overwrites same named file

#### Calculate Water Quality indices using waterquality ------------------------
Milford_2018 <- stack("Milford_Quantitative/Milford_L2A_Masked_20180709.tif") 

#Multiple Algorithms
Milford_2018_Chla <- wq_calc(raster_stack = Milford_2018, alg = "chlorophyll", sat = "sentinel2")

#Plot Chl-a indices
plot(Milford_2018_Chla)

#### Extract Values from Raster imagery from Shapefile--------------------------

#Create Spatial File from XY Coords
Milford_Data <- read.csv("Milford_Quantitative/Milford_WQ_Data_20180710.csv")

# Convert tabular data into spatial data
#Set sf with Long & Lat
Milford_Data <- Milford_Data %>% st_as_sf(coords = c("Long", "Lat"))

# Set CRS
Milford_Points = st_set_crs(Milford_Data, 4326)

#Check Lat Long
st_is_longlat(Milford_Points)

# Save Data as geopackage
st_write(Milford_Points, "Milford_Quantitative/Milford_GeoWQ_Data.gpkg")

# Extract values of raster using geospatial points
waterquality_data <- data.frame(Milford_Points, raster::extract(Milford_2018_Chla, Milford_Points))

# Map sampling points over a basemap
lake_extent <- st_read("Milford_Quantitative/Milford_Reservoir.shp")
Map_WQ_basemap(WQ_extent = lake_extent,
               sample_points = Milford_Points,
               WQ_parameter = "Chlorophyll_ugL",
               map_title= "Water Quality Map",
               points_style = "quantile",
               histogram = TRUE)

#Clean and Save Data
#Remove Geometry (Causes Problems)
waterquality_data <- subset(waterquality_data, select = -geometry)

# Remove two points with no data
waterquality_data <- waterquality_data[3:32,]

#Export Data
write.csv(waterquality_data, "Milford_Quantitative/Milford_GeoWQ_Data.csv")

#### Creating Statistical Models -----------------------------------------------

### THIS DATA IS SUBJECT TO ERRORS
# ****Imagery and in situ collect was done on two separate days!!!!
# Potential Mixed pixels
# Clouds/Poor AC methods

# Import Data 
waterquality_data <- read.csv("Milford_Quantitative/Milford_GeoWQ_Data.csv")
waterquality_data$Am092Bsub <- as.numeric(waterquality_data$Am092Bsub)
View(waterquality_data)

#### Basic Linear Model
#One Parameter & One Algorithm
extract_lm(parameter = "Chlorophyll_ugL", algorithm = "MM12NDCI", df = waterquality_data)

#### Linear Model with Cross-Validation
extract_lm_cv(parameter = "Chlorophyll_ugL", algorithm = "MM12NDCI",
              df = waterquality_data, train_method = "lm", control_method = "repeatedcv",
              folds = 3, nrepeats = 5)

#### Linear Model with Cross-Validation for multiple algorithms and parameters
# extract_lm_cv_all()
- `parameters` list of water quality parameters
- `df` data frame containing the values for parameter and algorithm arguments
- `train_method` A string specifying which classification or regression model to use (Default = "lm"). See ?caret::train for more details
- `control_method` A string specifying the resampling method (Default = "repeatedcv"). See ?caret::trainControl for more details
- `folds` the number of folds to be used in the cross validation model (Default = 3)
- `nrepeats` the number of iterations to be used in the cross validation model (Default = 5)

#Define Parameters
parameter <- c("Chlorophyll_ugL")

# Run Model
extract_lm_cv_all_results <- extract_lm_cv_all(parameters = parameter,
                                                 df = waterquality_data,
                                                 train_method = "lm",
                                                 control_method = "repeatedcv",
                                                 folds = 3,
                                                 nrepeats = 5)

#Show Results
View(extract_lm_cv_all_results)









