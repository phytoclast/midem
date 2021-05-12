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
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
demfolder <- 'data/Michigan/tif'
demlist <- list.files(demfolder)
index.dem<- readRDS('data/index.dem.RDS')
surfacemodelpath <- 'E:/IMG'
#surfacemodelpath <- 'data/IMG'
wkt.new <- crs(rast(paste0(surfacemodelpath,'/','605_640','.img')))
index <- st_read('data/Shapes/MI_LAS_Tiles.shp')
stands <- st_read('data/Shapes/stands.shp')
stands <- st_transform(stands, st_crs(index))
# plot(st_geometry(index))
# plot(st_geometry(stands), add=TRUE)
stands.folders <- unique(stands$folder)

for (i in 1:length(stands.folders)){#i = 1
  stand.select <- subset(stands, folder %in% stands.folders[i])
  index.select <- index[stand.select,] 
  if(nrow(index.select)>0){  
    f <- index.select$Name
    surface <- raster(paste0(surfacemodelpath,'/',f[1],'.IMG'))
    if(length(f) > 1){
      for (j in 2:length(f)) {
        surface1 <- raster(paste0(surfacemodelpath,'/',f[j],'.IMG'))
        surface <- merge(surface, surface1)
      }
    }
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
    if(length(dem.select$rname) > 1){
      for (j in 2:length(dem.select$rname)) {
        ground1 <- rast(paste0(demfolder,'/',dem.select$rname[j]))
        ground <- merge(ground, ground1)
      }
    }
    ground <- project(ground, surface)
    canopy <- surface - ground
    crs(canopy) <- wkt.new
    writeRaster(canopy, paste0('output/',stands.folders[i],'.tif'), overwrite=TRUE)
  }
}
plot(canopy)

