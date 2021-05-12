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
wkt.new <- crs(rast(paste0(surfacemodelpath,'/','605_640','.img')))
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
demfolder <- 'data/Michigan/tif'
demlist <- list.files(demfolder)
index.dem<- readRDS('data/index.dem.RDS')
#surfacemodelpath <- 'E:/IMG'
#index <- st_read('E:/Shapes/MI_LAS_Tiles.shp')
surfacemodelpath <- 'data/IMG'
#index <- st_read('data/Shapes/MI_LAS_Tiles.shp')
#f <- '605_640'
#surface <- rast(paste0(surfacemodelpath,'/',f,'.img'))
index <- list.files(surfacemodelpath)
for (i in 1:length(index)){#i = 2
f <- index[i]
f.tif <- str_replace(f, '.img','.tif')
surface <- raster(paste0(surfacemodelpath,'/',f))
writeRaster(surface, 'tmp/surface.tif', overwrite=TRUE)
surface <- rast('tmp/surface.tif')
crs(surface) <- 'epsg:6497'
ext.surface <- ext(surface)

dem.select <- subset(index.dem,
                       xmin <= ext.surface[2] &
                       xmax >= ext.surface[1] &
                       ymin <= ext.surface[4] &
                       ymax >= ext.surface[3])

ground <- rast(paste0(demfolder,'/',dem.select$rname[1]))
ground <- project(ground, surface)
canopy <- surface - ground
crs(canopy) <- wkt.new
writeRaster(canopy, paste0('output/',f.tif), overwrite=TRUE)
}
plot(canopy)

