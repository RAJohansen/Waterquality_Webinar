#################################################
### Detecting and Quantifying Water Quality   ###
### Using Satellite Remote Sensing Algorithms ###
#################################################

#### Qualitative WQ Monitoring ####
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
library(sen2r)
#install.packages(c("leafpm", "mapedit", "shiny", "shinyFiles", "shinydashboard", "shinyjs", "shinyWidgets",'geojsonlint'))

#### Use Sen2R for data acquisition and Pre-processing -------------------------

#Check D
check_sen2r_deps()

#Explore GUI for Sentinel-2 Data and Pre-Processing 
sen2r()


#### Calculate Water Quality indices using waterquality ------------------------
Milford_raster <- stack("Milford_Qualitative/Milford_BOA_20200916.tif") 

#Plot RGB Image
raster::plotRGB(Milford_raster,
                r=4,
                g=3,
                b= 2,
                stretch='lin')

#Single Algorithm
Milford_NDCI <- wq_calc(raster_stack = Milford_raster, alg = "MM12NDCI", sat = "sentinel2")

#Multiple Algorithms
Milford_Chla <- wq_calc(raster_stack = Milford_raster, alg = "chlorophyll", sat = "sentinel2")

### Generate Map of NDCI Algorithm  --------------------------------------------

# Create Map of Milford Reservoir
Map_WQ_raster(WQ_raster = Milford_NDCI,
              map_title= "Water Quality Map",
              raster_style = "quantile",
              histogram = TRUE)

#Import Milford Lake Boundary Shapefile (Optional)
Milford_Extent <- read_sf("Milford_Qualitative/Milford_Reservoir.shp")
Projection = "+proj=utm +zone=14 +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0"
Milford_Extent = st_transform(Milford_Extent, Projection)
Milford_Extent <- st_zm(Milford_Extent)

#Create Mask
Milford_NDCI_Mask <- raster::mask(Milford_NDCI,Milford_Extent)

# Create Map of Milford Reservoir Mask
Map_WQ_raster(WQ_raster = Milford_NDCI_Mask,
              map_title= "Water Quality Map",
              raster_style = "quantile",
              histogram = FALSE)

#Export Relative Index/Indices as Tif
writeRaster(x = Milford_NDCI_Mask,
            filename= "Milford_Qualitative/Milford_NDCI_Mask_20200916.tif", # save as a tif
            datatype ="FLT4S", # save as a float 4 significant digits
            overwrite = FALSE) #Overwrites same named file
