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
#surfacemodelpath <- 'E:/IMG'
#index <- st_read('E:/Shapes/MI_LAS_Tiles.shp')
surfacemodelpath <- 'data/IMG'
#wkt.new <- crs(rast(paste0(surfacemodelpath,'/','605_640','.img')))
wkt.new <- 'epsg:6497'
index.dem <- NA

f <- '605_640'
surface <- rast(paste0(surfacemodelpath,'/',f,'.img'))
demfolder <- 'data/Michigan/tif'
demlist <- list.files(demfolder)

for(i in 1:length(demlist))
{#i=29
demname <- demlist[i]
ground <- rast(paste0(demfolder,'/',demname))

xmin = ext(ground)[1]
xmax = ext(ground)[2]
ymin = ext(ground)[3]
ymax = ext(ground)[4]
LR <- c(xmax, ymin)
LL <- c(xmin, ymin)
UR <- c(xmax, ymax)
UL <- c(xmin, ymax)

ex <- data.frame(rname=c('LR','LL','UR','UL'),
                 xcoord=c(LR[1],LL[1],UR[1],UL[1]),
                 ycoord=c(LR[2],LL[2],UR[2],UL[2])
)

ex <- sf::st_as_sf(as.data.frame(ex), coords = c("xcoord","ycoord"), crs=crs(ground))
ex.trans <- (st_transform(ex,crs=wkt.new))
ex.trans <- as.data.frame(st_coordinates(ex.trans))
exmnx <- data.frame(rname=c(demname),
                    xmin = min(ex.trans$X),
                    xmax = max(ex.trans$X),
                    ymin = min(ex.trans$Y),
                    ymax = max(ex.trans$Y)
)
if(is.na(index.dem)[1]){
  index.dem <- exmnx
}else{
  index.dem <- rbind(index.dem, exmnx)
}

}
saveRDS(index.dem,'data/index.dem.RDS')
