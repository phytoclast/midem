library(dplyr)
library(stringr)
library(terra)
library(raster)
library(sf)
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
olddem <- rast('data/Michigan/tif/oldUSGS_13_n45w085.tif')

huronpatch <- rast('D:/GIS/DEM/10mDEM/HMNF/hnf_patch')
huronpatch <- project(huronpatch, olddem)

newdem <- merge(huronpatch, olddem)
writeRaster(newdem, 'data/Michigan/tif/USGS_13_n45w085.tif', overwrite=TRUE)
plot(newdem)

#fix CHM ---- 
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
demfolder <- 'data/Michigan/tif'
demlist <- list.files(demfolder)
index.dem<- readRDS('data/index.dem.RDS')
surfacemodelpath <- 'E:/IMG_Thinned'
#index <- st_read('E:/Shapes/MI_LAS_Tiles.shp')
#surfacemodelpath <- 'data/IMG'
#index <- st_read('data/Shapes/MI_LAS_Tiles.shp')
#f <- '605_640'
#surface <- rast(paste0(surfacemodelpath,'/',f,'.img'))
index <- list.files(surfacemodelpath)
wkt.new <- crs(rast(paste0(surfacemodelpath,'/',index[1])))
#for (i in 1:length(index)){#i = 2
  f <- '645_450_Thinned.img'
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
  if(length(dem.select$rname) > 1){
    for (j in 2:length(dem.select$rname)) {
      ground1 <- rast(paste0(demfolder,'/',dem.select$rname[j]))
      ground <- merge(ground, ground1)
    }
  }
  ground <- project(ground, surface)
  canopy <- surface - ground 
  crs(canopy) <- wkt.new
  writeRaster(canopy, paste0('output/thinned/',f.tif), overwrite=TRUE)
#}
plot(canopy)

