library(sf)
library(terra)
library(raster)
library(minpack.lm)
library(growthmodels)
library(dplyr)

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
dem <- rast('D:/GIS/DEM/hmnfdem30.tif'); names(dem) = 'dem'

d3.max <- dem %>% aggregate(fact=3, fun='max', na.rm=T) %>% focal(fun='max', na.rm=T) %>% 
  focal(fun='mean', na.rm=T) %>% project(dem)
d3.min <- dem %>% aggregate(fact=3, fun='min', na.rm=T) %>% focal(fun='min', na.rm=T) %>% 
  focal(fun='mean', na.rm=T) %>% project(dem)


d10.max <- dem %>% aggregate(fact=10, fun='max', na.rm=T) %>% focal(fun='max', na.rm=T) %>% 
  focal(fun='mean', na.rm=T) %>% project(dem)
d10.min <- dem %>% aggregate(fact=10, fun='min', na.rm=T) %>% focal(fun='min', na.rm=T) %>% 
  focal(fun='mean', na.rm=T) %>% project(dem)


d33.max <- dem %>% aggregate(fact=33, fun='max', na.rm=T) %>% focal(fun='max', na.rm=T) %>% 
  focal(fun='mean', na.rm=T) %>% project(dem)
d33.min <- dem %>% aggregate(fact=33, fun='min', na.rm=T) %>% focal(fun='min', na.rm=T) %>% 
  focal(fun='mean', na.rm=T) %>% project(dem)


d3.mean <- dem %>% aggregate(fact=3, fun='mean', na.rm=T) %>% focal(fun='mean', na.rm=T) %>% 
  focal(fun='mean', na.rm=T) %>% project(dem)
d10.mean <- dem %>% aggregate(fact=10, fun='mean', na.rm=T) %>% focal(fun='mean', na.rm=T) %>% 
  focal(fun='mean', na.rm=T) %>% project(dem)
d33.mean <- dem %>% aggregate(fact=33, fun='mean', na.rm=T) %>% focal(fun='mean', na.rm=T) %>% 
  focal(fun='mean', na.rm=T) %>% project(dem)



tpi3l <- (dem - d3.mean)/(d3.mean - d3.min + 100)*(dem - d3.mean <0)
tpi10l <- (dem - d10.mean)/(d10.mean - d10.min + 100)*(dem - d10.mean <0)
tpi33l <- (dem - d33.mean)/(d33.mean - d33.min + 100)*(dem - d33.mean <0)

tpi3h <- (dem - d3.mean)/(d3.max - d3.mean + 100)*(dem - d3.mean >0)
tpi10h <- (dem - d10.mean)/(d10.max - d10.mean + 100)*(dem - d10.mean >0)
tpi33h <- (dem - d33.mean)/(d33.max - d33.mean + 100)*(dem - d33.mean >0)

tpi3 <- tpi3h+tpi3l
tpi10 <- tpi10h+tpi10l
tpi33 <- tpi33h+tpi33l




writeRaster(tpi3, 'output/tpi3.tif', overwrite=T)
writeRaster(tpi10, 'output/tpi10.tif', overwrite=T)
writeRaster(tpi33, 'output/tpi33.tif', overwrite=T)

tpi.mean <- mean(tpi3,tpi10,tpi33)
writeRaster(tpi.mean, 'output/tpi.mean.tif', overwrite=T)
10000*30
tpi.wtd <- 
  (tpi3*(d3.max - d3.min + .1)+
     tpi10*(d10.max - d10.min + .1)+
     tpi33*(d33.max - d33.min + .1))/
  ((d3.max - d3.min + .1)+(d10.max - d10.min + .1)+(d33.max - d33.min + .1))



writeRaster(tpi.wtd, 'output/tpi.wtd.tif', overwrite=T)






tpi.mean <- mean(tpi2, tpi5, tpi10, tpi30)
writeRaster(tpi.mean, 'output/tpi.mean.tif')
tpi.wtd <- 
  (tpi2*(d2.max - d2.min + 2)+
     tpi5*(d5.max - d5.min + 2)+
     tpi10*(d10.max - d10.min + 2)+
     tpi30*(d30.max - d30.min + 2))/
  ((d2.max - d2.min + 2)+(d5.max - d5.min + 2)+(d10.max - d10.min + 2)+(d30.max - d30.min + 2))
  


writeRaster(tpi.wtd, 'output/tpi.wtd.tif')



tpi2m <- (dem - d2.mean)/(d10.max - d10.min + 2)
tpi5m <- (dem - d5.mean)/(d10.max - d10.min + 2)
tpi10m <- (dem - d10.mean)/(d10.max - d10.min + 2)
tpi30m <- (dem - d30.mean)/(d10.max - d10.min + 2)
writeRaster(tpi2m, 'output/tpi2m.tif')
writeRaster(tpi5m, 'output/tpi5m.tif')
writeRaster(tpi10m, 'output/tpi10m.tif')
writeRaster(tpi30m, 'output/tpi30m.tif')

