library(dplyr)
library(stringr)
library(terra)
library(raster)
library(sf)
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
MLRA93 <- rast('output/UP/MLRA93.tif')
MLRA94B <- rast('output/UP/MLRA94B.tif')
MLRA94AC <- rast('output/NL/MLRA94AC.tif')
MLRA96 <- rast('output/NL/MLRA96.tif')
MLRA97 <- rast('output/SL/MLRA97.tif')
MLRA98 <- rast('output/SL/MLRA98.tif')
MLRA99N <- rast('output/SL/MLRA99N.tif')
MLRA99S <- rast('output/SL/MLRA99S.tif')

CHM <- MLRA98
CHM[CHM > 50] <- NA
CHM[CHM < 0] <- 0
CHM.100 <- aggregate(CHM, fact = 20, fun='max', na.rm = TRUE)
writeRaster(CHM.100, 'output/CHM.100.tif', overwrite=T)
CHM <- MLRA97
CHM[CHM > 50] <- NA
CHM[CHM < 0] <- 0
CHM.100.2 <- aggregate(CHM, fact = 20, fun='max', na.rm = TRUE)
CHM.100 <- merge(CHM.100, CHM.100.2)
writeRaster(CHM.100, 'output/CHM.100.tif', overwrite=T)
CHM <- MLRA99N
CHM[CHM > 50] <- NA
CHM[CHM < 0] <- 0
CHM.100.2 <- aggregate(CHM, fact = 20, fun='max', na.rm = TRUE)
CHM.100 <- merge(CHM.100, CHM.100.2)
writeRaster(CHM.100, 'output/CHM.100.tif', overwrite=T)
CHM <- MLRA99S
CHM[CHM > 50] <- NA
CHM[CHM < 0] <- 0
CHM.100.2 <- aggregate(CHM, fact = 20, fun='max', na.rm = TRUE)
CHM.100 <- merge(CHM.100, CHM.100.2)
writeRaster(CHM.100, 'output/CHM.100.tif', overwrite=T)
CHM <- MLRA96
CHM[CHM > 50] <- NA
CHM[CHM < 0] <- 0
CHM.100.2 <- aggregate(CHM, fact = 20, fun='max', na.rm = TRUE)
CHM.100 <- merge(CHM.100, CHM.100.2)
writeRaster(CHM.100, 'output/CHM.100.tif', overwrite=T)

CHM.100 <- rast('output/CHM.100.tif')
CHM <- rast('output/NL/MLRA94AC.tif')
CHM[CHM > 50] <- NA
CHM[CHM < 0] <- 0
CHM.100.2 <- aggregate(CHM, fact = 20, fun='max', na.rm = TRUE)
CHM.100 <- merge(CHM.100, CHM.100.2)
writeRaster(CHM.100, 'output/CHM.100.tif', overwrite=T)

CHM.100 <- rast('output/CHM.100.tif')
CHM <- rast('output/UP/MLRA94B.tif')
CHM[CHM > 50] <- NA
CHM[CHM < 0] <- 0
CHM.100.2 <- aggregate(CHM, fact = 20, fun='max', na.rm = TRUE)
CHM.100 <- merge(CHM.100, CHM.100.2)
writeRaster(CHM.100, 'output/CHM.100.tif', overwrite=T)

CHM.100 <- rast('output/CHM.100.tif')
CHM <- rast('output/UP/MLRA93.tif')
CHM[CHM > 50] <- NA
CHM[CHM < 0] <- 0
CHM.100.2 <- aggregate(CHM, fact = 20, fun='max', na.rm = TRUE)
CHM.100 <- merge(CHM.100, CHM.100.2)
writeRaster(CHM.100, 'output/CHM.100.tif', overwrite=T)

