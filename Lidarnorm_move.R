library(lidR)
library(rgdal)
library(rlas)
library(mapview)
library(progress)
library(future)
library(viridis)
library(dplyr)
library(stringr)
library(terra)
library(sf)
library(filesstrings)
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
index.dem<- readRDS('data/index.dem.RDS')
canopymodelpath <- 'output/thinned'
index <- st_read('data/Shapes/MI_LAS_Tiles.shp')
index <- subset(index, !is.na(newfolder))
index$f <- paste0(index$Name,'_Thinned.tif')
for (i in 1:nrow(index)){#i=1
  oldfolder <- index[i,]$Desc
  newfolder <- index[i,]$newfolder
f <- index[i,]$f
file.move(paste0('output/',oldfolder,'/',f), paste0('output/',newfolder))
}