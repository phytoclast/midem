# ----
library(sf)
library(terra)
library(raster)
library(minpack.lm)
library(growthmodels)
library(dplyr)

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))


chm1 <- rast('output/NL/MLRA96.tif')
chm2 <- rast('output/NL/MLRA94AC.tif')

chm <- merge(chm1, chm2)
chm[chm<0] <- 0
chm[chm>50] <- 50

chm.30.max <- aggregate(chm, fact=30/5, na.rm=T, fun='max')
chm.30.mean <- aggregate(chm, fact=30/5, na.rm=T, fun='mean')
EV_Code <- read.csv('data/hmnfstandpoints/EV_Code.txt')
EV_Code$EV_COMMON_NAME <- tolower(EV_Code$EV_COMMON_NAME)
EV_Code <- unique(EV_Code[,c('EV_COMMON_NAME','First_EV_CODE')])
EV_Code <- EV_Code %>% group_by(First_EV_CODE) %>% summarise(type = first(EV_COMMON_NAME))
mukey <- rast('D:/GIS/SOIL/2021/northeast_gssurgo30m.tif'); names(mukey) = 'mukey'


# hmnf <- read.csv('data/hmnfstandpoints/SiteIndexRandomPoints_10perStand.txt')
# hmnf <- st_as_sf(hmnf, coords=c(x='Long', y='Lat'), crs='epsg:4326', remove = FALSE)
hmnf <- read_sf('data/hmnfstandpoints/hmnfstandpoints.join3.shp')
hmnf.mi <- sf::st_transform(hmnf, crs=crs(chm))
hmnf.albers <- sf::st_transform(hmnf, crs=crs(mukey))
hmnf.chm.max <- terra::extract(chm.30.max, vect(hmnf.mi))
hmnf.chm.mean <- terra::extract(chm.30.mean, vect(hmnf.mi))
hmnf.chm.mukey <- terra::extract(mukey, vect(hmnf.albers))
colnames(hmnf.chm.max) <- c('id','chm.max');colnames(hmnf.chm.mean) <- c('id','chm.mean');colnames(hmnf.chm.mukey) <- c('id','mukey')
hmnf.chm <- cbind(hmnf.mi, subset(hmnf.chm.max, chm=hmnf.chm[,2]))
hmnf.chm <- cbind(hmnf.chm, subset(hmnf.chm.mean, chm=hmnf.chm[,2]))
hmnf.chm <- cbind(hmnf.chm, subset(hmnf.chm.mukey, chm=hmnf.chm[,2]))
hmnf.chm$EV_CODE <- as.numeric(hmnf.chm$EV_CODE)
hmnf.chm <- hmnf.chm %>% left_join(EV_Code, by=c('EV_CODE'='First_EV_CODE'))
hmnf.chm$YEAR_OF_ORIGIN <- hmnf.chm$YEAR_OF_OR
hmnf.chm$SITE_INDEX_SPP <- hmnf.chm$SITE_IND_1
hmnf.chm$SITE_INDEX_REF <- hmnf.chm$SITE_IND_2
hmnf.chm$age <- 2021-hmnf.chm$YEAR_OF_ORIGIN
hmnf.chm$lmapunitiid <- hmnf.chm$mukey
colnames(hmnf.chm)
hmnf.2015 <- read_sf('data/hmnfstandpoints/hmnfstandpoints.join2.shp')
hmnf.2015$YEAR_OF_ORIGIN2 <- hmnf.2015$YEAR_OF_OR
hmnf.2015$EV_CODE2 <- hmnf.2015$EV_CODE
hmnf.2015 <- st_drop_geometry(hmnf.2015)
hmnf.2015 <- unique(subset(hmnf.2015, select=c(CID, YEAR_OF_ORIGIN2, EV_CODE2)))
hmnf.chm <- left_join(hmnf.chm, hmnf.2015)

hmnf.chm <-  subset(hmnf.chm, select=c("Lat","Long","EV_CODE","YEAR_OF_ORIGIN","SITE_INDEX","SITE_INDEX_SPP","SITE_INDEX_REF",
                                       "lmapunitiid", "chm.max", "chm.mean","type","age","YEAR_OF_ORIGIN2", "EV_CODE2")) %>% st_drop_geometry()
hmnf.chm <- subset(hmnf.chm, YEAR_OF_ORIGIN > 0)
hmnf.chm$age <- ifelse(hmnf.chm$YEAR_OF_ORIGIN > 2021 | (hmnf.chm$YEAR_OF_ORIGIN > 2020 & hmnf.chm$chm.mean > 5), 2021 - hmnf.chm$YEAR_OF_ORIGIN2, hmnf.chm$age)

write.csv(hmnf.chm, 'data/hmnfstandpoints/hmnf.chm.csv', row.names = F)
saveRDS(hmnf.chm, 'data/hmnfstandpoints/hmnf.chm.RDS')
#############################
#############################
library(sf)
library(terra)
library(raster)
library(minpack.lm)
library(growthmodels)
library(dplyr)
library(ggplot2)
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

hmnf.chm <- readRDS('data/hmnfstandpoints/hmnf.chm.RDS')

hmnf.chm$chm <- hmnf.chm$chm.max

ggplot(data=hmnf.chm)+
geom_point(aes(x=age, y=chm), col='red', alpha=0.02)

redpine <- subset(hmnf.chm, EV_CODE %in% 2 & !is.na(chm) & !is.na(age))
ggplot(data=redpine)+
  geom_point(aes(x=age, y=chm), col='red', alpha=0.02)

x=redpine$age
y=redpine$chm

mod1 <- nlsLM(y ~ SSgompertz(x, Asym, b2, b3))
x1 <- 0:200
pred <-  predict(mod1, list(x = x1))
ggplot()+
  geom_point(aes(x=age, y=chm),data=redpine, col='red', alpha=0.02)+
  geom_line(aes(x=x1,y=pred), col='black')

summary(mod)

redpine <- subset(hmnf.chm, EV_CODE %in% 2 & !is.na(chm) & !is.na(age) & age >= 5 & chm >=0.1)
x <- redpine$age
y <- redpine$chm
x=c(x,x*0-1)
y=c(y,y*0)
mod1 <- nlsLM(y ~ growthmodels::gompertz(x, alpha, beta, k), start = list(alpha = 30, beta = 1, k = 1))#can swap out start values for fixed model values
summary(mod1)
mod2 <- nlsLM(y ~ growthmodels::chapmanRichards(x, alpha, beta, k, m), start = list(alpha=30, beta=.1, k=.1, m=.1),
              control = list(maxiter = 150))
summary(mod2)
mod3 <- nlsLM(y ~ growthmodels::generalisedRichard(x, A, U, k, m, beta), start = list(A=0, U=30, k=.1, m=.1, beta=.1),
              control = list(maxiter = 150))
summary(mod3)


x1 <- 0:200+0
pred1 <-  predict(mod1, list(x = x1))
pred2 <-  predict(mod2, list(x = x1))
pred3 <-  predict(mod3, list(x = x1))
ggplot()+
  geom_point(aes(x=x, y=y), col='red', alpha=0.02)+
  geom_line(aes(x=x1,y=pred1), col='yellow', size=2)+
  geom_line(aes(x=x1,y=pred2), col='black')+
  geom_line(aes(x=x1,y=pred3), col='blue')


#------------------
jackpine <- subset(hmnf.chm, EV_CODE %in% 1 & !is.na(chm) & !is.na(age) & age >= 5 & chm >=0.1)
redpine <- subset(hmnf.chm, EV_CODE %in% 2 & !is.na(chm) & !is.na(age) & age >= 5 & chm >=0.1)
whitepine <- subset(hmnf.chm, EV_CODE %in% 3 & !is.na(chm) & !is.na(age) & age >= 5 & chm >=0.1)
aspen <- subset(hmnf.chm, EV_CODE %in% c(91,92,93) & !is.na(chm) & !is.na(age) & age >= 5 & chm >=0.1)
oak <- subset(hmnf.chm, EV_CODE %in% c(53,54,55,59,63) & !is.na(chm) & !is.na(age) & age >= 5 & chm >=0.1)
maple <- subset(hmnf.chm, EV_CODE %in% c(71,76,81,82,84,85,87) & !is.na(chm) & !is.na(age) & age >= 5 & chm >=0.1)

x1 <- jackpine$age
y1 <- jackpine$chm.max
x1=c(x1,x1*0-1)
y1=c(y1,y1*0)
x2 <- redpine$age
y2 <- redpine$chm.max
x2=c(x2,x2*0-1)
y2=c(y2,y2*0)

x3 <- whitepine$age
y3 <- whitepine$chm.max
x3=c(x3,x3*0-1)
y3=c(y3,y3*0)
x4 <- aspen$age
y4 <- aspen$chm.max
x4=c(x4,x4*0-1)
y4=c(y4,y4*0)

x5 <- oak$age
y5 <- oak$chm.max
x5=c(x5,x5*0-1)
y5=c(y5,y5*0)
x6 <- maple$age
y6 <- maple$chm.max
x6=c(x6,x6*0-1)
y6=c(y6,y6*0)


mod1 <- nlsLM(y1 ~ growthmodels::gompertz(x1, alpha, beta, k), start = list(alpha = 25, beta = 4, k = .05))
summary(mod1)
mod2 <- nlsLM(y2 ~ growthmodels::gompertz(x2, alpha, beta, k), start = list(alpha = 25, beta = 4, k = .05))
summary(mod2)
mod3 <- nlsLM(y3 ~ growthmodels::gompertz(x3, alpha, beta, k), start = list(alpha = 25, beta = 4, k = .05))
summary(mod1)
mod4 <- nlsLM(y4 ~ growthmodels::gompertz(x4, alpha, beta, k), start = list(alpha = 25, beta = 4, k = .05))
summary(mod2)
mod5 <- nlsLM(y5 ~ growthmodels::gompertz(x5, alpha, beta, k), start = list(alpha = 25, beta = 4, k = .05))
summary(mod1)
mod6 <- nlsLM(y6 ~ growthmodels::gompertz(x6, alpha, beta, k), start = list(alpha = 25, beta = 4, k = .05))
summary(mod2)

x0 <- 0:200+0
pred1 <-  predict(mod1, list(x1 = x0))
pred2 <-  predict(mod2, list(x2 = x0))
pred3 <-  predict(mod3, list(x3 = x0))
pred4 <-  predict(mod4, list(x4 = x0))
pred5 <-  predict(mod5, list(x5 = x0))
pred6 <-  predict(mod6, list(x6 = x0))

ggplot()+
  #geom_point(aes(x=x1, y=y1), col='red', alpha=0.02)+
  #geom_point(aes(x=x2, y=y2), col='blue', alpha=0.02)+
  geom_line(aes(x=x0,y=pred1), col='red', size=2)+
  geom_line(aes(x=x0,y=pred2), col='blue', size=2)
  

ggplot()+
  #geom_point(aes(x=x3, y=y3), col='red', alpha=0.02)+
  #geom_point(aes(x=x4, y=y4), col='blue', alpha=0.02)+
  geom_line(aes(x=x0,y=pred3), col='red', size=2)+
  geom_line(aes(x=x0,y=pred4), col='blue', size=2)

ggplot()+
  #geom_point(aes(x=x5, y=y5), col='red', alpha=0.02)+
  #geom_point(aes(x=x6, y=y6), col='blue', alpha=0.02)+
  geom_line(aes(x=x0,y=pred5), col='red', size=2)+
  geom_line(aes(x=x0,y=pred6), col='blue', size=2)

ggplot()+
  geom_line(aes(x=x0,y=pred1, col='jack pine'), size=1.5)+
  geom_line(aes(x=x0,y=pred2, col='red pine'), size=1.5)+
  geom_line(aes(x=x0,y=pred3, col='white pine'), size=1.5)+
  geom_line(aes(x=x0,y=pred4, col='aspen'), size=1.5)+
  geom_line(aes(x=x0,y=pred5, col='oak'), size=1.5)+
  geom_line(aes(x=x0,y=pred6, col='maple'), size=1.5)+
  scale_color_manual(   values = c('jack pine'='green','red pine'='darkgreen','white pine'='darkcyan',
                                   'aspen'='khaki4','oak'='orange','maple'='red'))+
  scale_x_continuous(name = 'age (years)')+
  scale_y_continuous(name = 'canopy height (meters)')
  
#############################
#acquire topographic covariates
#############################
library(sf)
library(terra)
library(raster)
library(minpack.lm)
library(growthmodels)
library(dplyr)
library(ggplot2)
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

hmnf.chm <- readRDS('data/hmnfstandpoints/hmnf.chm.RDS')

hmnf.chm$chm <- hmnf.chm$chm.max

solar <- rast('D:/GIS/DEM/hmnfsolarmean.tif'); names(solar) = 'solar'
toi <- rast('D:/GIS/DEM/hmnf-openess.tif'); names(toi) = 'toi'
toip <- rast('D:/GIS/DEM/hmnf-pos-openess.tif'); names(toip) = 'toip'
toin <- rast('D:/GIS/DEM/hmnf-neg-openess.tif'); names(toin) = 'toin'
toip100 <- rast('D:/GIS/DEM/hmnf-pos-openess100.tif'); names(toip100) = 'toip100'
toin100 <- rast('D:/GIS/DEM/hmnf-neg-openess100.tif'); names(toin100) = 'toin100'
twi <- rast('D:/GIS/DEM/hmnftwi.tif'); names(twi) = 'twi'
tpi <- rast('D:/GIS/DEM/hmnftopographicpositionindex.tif'); names(tpi) = 'tpi'
dem <- rast('D:/GIS/DEM/hmnfdem30.tif'); names(dem) = 'dem'
slope <- rast('D:/GIS/DEM/hmnfslope.tif'); names(slope) = 'slope'#radians
aspect <- rast('D:/GIS/DEM/hmnfaspect.tif'); names(aspect) = 'aspect'#radians
bt <- rast('D:/scripts/snow/output/bt.90.alt.tif'); names(bt) = 'bt'
tgs <- rast('D:/scripts/snow/output/tgs.90.alt.tif'); names(tgs) = 'tgs'
ppt <- rast('C:/a/Ecological_Sites/GIS/Climate/PRISM2010/P/p0112/w001001.adf'); names(ppt) = 'ppt'
p1 <- rast('C:/a/Ecological_Sites/GIS/Climate/PRISM2010/seasons/p1/w001001.adf'); names(p1) = 'p1'
p2 <- rast('C:/a/Ecological_Sites/GIS/Climate/PRISM2010/seasons/p2/w001001.adf'); names(p2) = 'p2'
p3 <- rast('C:/a/Ecological_Sites/GIS/Climate/PRISM2010/seasons/p3/w001001.adf'); names(p3) = 'p3'
p4 <- rast('C:/a/Ecological_Sites/GIS/Climate/PRISM2010/seasons/p4/w001001.adf'); names(p4) = 'p4'
slope.500 <- rast('D:/GIS/DEM/slope.500.tif'); names(slope.500) = 'slope.500'

bt <- project(bt, dem)
tgs <- project(tgs, dem)
ppt <- project(ppt, dem)
p1 <- project(p1, dem, method='bilinear')
p2 <- project(p2, dem, method='bilinear')
p3 <- project(p3, dem, method='bilinear')
p4 <- project(p4, dem, method='bilinear')

brk <- c(solar, toi, toip, toin, toip100, toin100, twi, tpi, dem, slope, aspect, slope.500, bt, tgs, ppt,p1,p2,p3,p4)

hmnf.pts <- subset(hmnf.chm, select=c(Lat, Long))

hmnf.pts <- st_as_sf(hmnf.pts, coords=c(x='Long', y='Lat'), crs='epsg:4326', remove = FALSE)
hmnf.pts <- sf::st_transform(hmnf.pts, crs=crs(dem))
hmnf.pts <- terra::extract(brk, vect(hmnf.pts))
hmnf.chm1 <- cbind(hmnf.chm, hmnf.pts)
hmnf.chm1$ID <- NULL

write.csv(hmnf.chm1, 'output/hmnf.chm1.csv', row.names = F)
#############################
#acquire soil covariates
#############################

library(sf)
library(terra)
library(raster)
library(minpack.lm)
library(growthmodels)
library(dplyr)
library(ggplot2)
library(Hmisc)
library(ranger)
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

hmnf.chm <- read.csv('output/hmnf.chm1.csv')

s <- read.csv('C:/workspace/SoilSorts/fy2021-refresh/s.csv')
colnames(s)

s$spodosols <- ifelse(s$taxorder %in% 'spodosols', 1,0)
s$spodic <- ifelse(grepl('spodic',s$taxsubgrp)|s$taxorder %in% 'spodosols', 1,0)
s$Bhs <- ifelse(s$Bhs %in% 'yes', 1,0)
s$spodosols.up <- s$spodosols*s$Water_Table>100
s$spodic.up <- s$spodic*s$Water_Table
s$Bhs.up <- s$Bhs*s$Water_Table

#fill missing awc data ---- 
s.train <- subset(s, !is.na(T150_AWC)&!is.na(T50_sand)&!is.na(T150_sand)&!is.na(T50_clay)&!is.na(T150_clay)&!is.na(T50_OM)&!is.na(T150_OM)&!is.na(T50_pH)&!is.na(Water_Table)&!is.na(spodic)&!is.na(spodosols)&!is.na(Bhs)&!is.na(carbdepth)&!is.na(humicdepth))
rf.awc <- ranger(T150_AWC ~ T50_sand+T150_sand+T50_clay+T150_clay+T50_OM+T150_OM+T50_pH+Water_Table+spodic+spodosols+Bhs+carbdepth+humicdepth+rockdepth,
                 data=s.train, num.trees=500, sample.fraction =0.02, 
                 write.forest = TRUE, importance = 'impurity')
s.predict <- subset(s, is.na(T150_AWC)&!is.na(T50_sand)&!is.na(T150_sand)&!is.na(T50_clay)&!is.na(T150_clay)&!is.na(T50_OM)&!is.na(T150_OM)&!is.na(T50_pH)&!is.na(Water_Table)&!is.na(spodic)&!is.na(spodosols)&!is.na(Bhs)&!is.na(carbdepth)&!is.na(humicdepth))

s.predict$new_AWC <- predictions(predict(rf.awc, data=s.predict))
s <- left_join(s, unique(s.predict[,c('coiid', 'new_AWC')]))
s$T150_AWC <- ifelse(is.na(s$T150_AWC),s$new_AWC,s$T150_AWC );s$new_AWC <- NULL
#fill missing pH data ----
s$taxsubgrp <- as.factor(s$taxsubgrp) 
s.train <- subset(s, !is.na(T50_sand)&!is.na(T150_sand)&!is.na(T50_clay)&!is.na(T150_clay)&!is.na(T50_OM)&!is.na(T150_OM)&!is.na(T50_pH)&!is.na(Water_Table)&!is.na(spodic)&!is.na(taxsubgrp)&!is.na(Bhs)&!is.na(carbdepth)&!is.na(humicdepth))
rf.ph <- ranger(T50_pH ~ T50_sand+T150_sand+T50_clay+T150_clay+T50_OM+T150_OM+Water_Table+spodic+taxsubgrp+Bhs+carbdepth+humicdepth+rockdepth,
                 data=s.train, num.trees=500, sample.fraction =0.02, 
                 write.forest = TRUE, importance = 'impurity')
s.predict <- subset(s, !is.na(T50_sand)&!is.na(T150_sand)&!is.na(T50_clay)&!is.na(T150_clay)&!is.na(T50_OM)&!is.na(T150_OM)&!is.na(Water_Table)&!is.na(spodic)&!is.na(taxsubgrp)&!is.na(Bhs)&!is.na(carbdepth)&!is.na(humicdepth))

s.predict$new_pH <- predictions(predict(rf.ph, data=s.predict))
s <- left_join(s, unique(s.predict[,c('coiid', 'new_pH')]))
s$T50_pH <- ifelse(is.na(s$T50_pH),s$new_pH,s$T50_pH );s$new_pH <- NULL
write.csv(s, 'C:/workspace/SoilSorts/fy2021-refresh/s.cleaned.csv', row.names = FALSE)
s <- read.csv('C:/workspace/SoilSorts/fy2021-refresh/s.cleaned.csv')
s.mu <- s %>% group_by(lmapunitiid=lmapunitiid) %>% dplyr::summarise(
  T150_AWC = wtd.mean(T150_AWC, weights=comppct_r, na.rm=T),
  T50_sand = wtd.mean(T50_sand, weights=comppct_r, na.rm=T),
  T150_sand = wtd.mean(T150_sand, weights=comppct_r, na.rm=T),
  T50_clay = wtd.mean(T50_clay, weights=comppct_r, na.rm=T),
  T150_clay = wtd.mean(T150_clay, weights=comppct_r, na.rm=T),
  T50_OM = wtd.mean(T50_OM, weights=comppct_r, na.rm=T),
  T150_OM = wtd.mean(T150_OM, weights=comppct_r, na.rm=T),
  T50_pH = wtd.mean(T50_pH, weights=comppct_r, na.rm=T),
  Water_Table = wtd.mean(Water_Table, weights=comppct_r, na.rm=T),
  spodic = wtd.mean(spodic, weights=comppct_r, na.rm=T),
  spodosols = wtd.mean(spodosols, weights=comppct_r, na.rm=T),
  Bhs = wtd.mean(Bhs, weights=comppct_r, na.rm=T),
  carbdepth = wtd.mean(carbdepth, weights=comppct_r, na.rm=T),
  humicdepth = wtd.mean(humicdepth, weights=comppct_r, na.rm=T)
)

colnames(s.mu)




mufilter <- unique(hmnf.chm$lmapunitiid)
hmnf.chm1 <- left_join(hmnf.chm, s.mu)
write.csv(hmnf.chm1, 'output/hmnf.chm2.csv', row.names = F)

#############################
#analysis
#############################

library(sf)
library(terra)
library(raster)
library(minpack.lm)
library(growthmodels)
library(plyr)
library(dplyr)
library(ggplot2)
library(Hmisc)
library(MASS)
library(ranger)
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

hmnf.chm <- read.csv('output/hmnf.chm2.csv')
hmnf.chm <- subset(hmnf.chm, !is.na(chm.max)&!is.na(T50_pH)&!is.na(T150_AWC))
hmnf.chm$jackpine <- ifelse(hmnf.chm$EV_CODE %in% 1,1,0)
hmnf.chm$redpine <- ifelse(hmnf.chm$EV_CODE %in% 2,1,0)
hmnf.chm$whitepine <- ifelse(hmnf.chm$EV_CODE %in% 3,1,0)
hmnf.chm$aspen <- ifelse(hmnf.chm$EV_CODE %in% c(91,92,93),1,0)
hmnf.chm$oak <- ifelse(hmnf.chm$EV_CODE %in% c(53,54,55,59,63),1,0)
hmnf.chm$maple <- ifelse(hmnf.chm$EV_CODE %in% c(71,76,81,82,84,85,87),1,0)
hmnf.chm$south <- hmnf.chm$slope^2+cos(hmnf.chm$aspect+3.141592)
hmnf.chm$west <- hmnf.chm$slope^2+cos(hmnf.chm$aspect+3.141592*1.5)
hmnf.chm$southwest <- hmnf.chm$slope^2+cos(hmnf.chm$aspect+3.141592*1.25)
hmnf.chm$southeast <- hmnf.chm$slope^2+cos(hmnf.chm$aspect+3.141592*0.75)
hmnf.chm$wet <- 25/(hmnf.chm$Water_Table+25)
hmnf.chm$hydric <- pmin(1,(pmax(50-hmnf.chm$Water_Table, 0)/25))
hmnf.chm$moist <- pmin(1,(pmax(200-hmnf.chm$Water_Table, 0)/100))
hmnf.chm$hilly <- (hmnf.chm$slope.500 > 0.05)*1
hmnf.chm$toip100

c.table <- as.data.frame(cor(hmnf.chm[,c("chm.max", "chm.mean","age","solar","wet","hilly",
                                         "toi","toip","toin","toip100","toin100","twi","tpi","dem","slope","south","west","southwest","southeast","bt","tgs",
"ppt","p1","p2","p3","p4","T150_AWC","T50_sand","T150_sand","T50_clay","T150_clay","T50_OM","spodic","spodosols","Bhs",
"T150_OM","T50_pH","Water_Table", "jackpine", "redpine", "whitepine", "aspen", "oak", "maple")], use='complete.obs'
))



#write.csv(c.table, 'output/c.table.csv', row.names = F)
model <- lm(chm.max~
              solar+
              #south+
              #west+
              #age+
              #I(age^0.5)+
              #twi+
              #toi+
              #slope+
              toip++
              #toin+
              #tpi+
              bt+
              tgs+
              ppt+
              T150_AWC+
              T50_sand+
              T150_sand+
              #T50_clay+
              #T150_clay+
              T50_OM+
              T150_OM+
              T50_pH+
              Water_Table+
              wet+
              #spodic+
              (spodosols*wet)+
              (humicdepth*wet)+
              #carbdepth+
              #Bhs+
              jackpine+
              redpine+
              whitepine+
              aspen+
              oak+
              maple:I(age^0.5)
            , data=hmnf.chm2)

stepAIC(model)
summary(model)
hmnf.chm2 <-  subset(hmnf.chm, Water_Table <=5000)
model <- lm(chm.max~
              solar:I(age^0.5)+
              #south:I(age^0.5)+
              #west:I(age^0.5)+
              #age:I(age^0.5)+
              #I(age^0.5)+
              #twi:I(age^0.5)+
              #toi:I(age^0.5)+
              hilly:I(age^0.5)+
              toip:I(age^0.5)+
              #toip100:I(age^0.5)+
              #toin:I(age^0.5)+
              #tpi:I(age^0.5)+
              #bt:I(age^0.5)+
              tgs:I(age^0.5)+
              ppt:I(age^0.5)+
              T150_AWC:I(age^0.5)+
              T50_sand:I(age^0.5)+
              #T150_sand:I(age^0.5)+
              #T50_clay:I(age^0.5)+
              #T150_clay:I(age^0.5)+
              T50_OM:I(age^0.5)+
              T150_OM:I(age^0.5)+
              T50_pH:I(age^0.5)+
              #Water_Table:I(age^0.5)+
              wet:I(age^0.5)+
              #spodic:I(age^0.5)+
              (spodosols*wet):I(age^0.5)+
              #(humicdepth*wet):I(age^0.5)+
              #carbdepth:I(age^0.5)+
              #Bhs:I(age^0.5)+
              jackpine:I(age^0.5)+
              redpine:I(age^0.5)+
              whitepine:I(age^0.5)+
              aspen:I(age^0.5)+
              oak:I(age^0.5)+
              maple:I(age^0.5)
            , data=hmnf.chm2)

stepAIC(model)
summary(model)

jackpine <- subset(hmnf.chm, jackpine %in% 1)
redpine <- subset(hmnf.chm, redpine %in% 1)
whitepine <- subset(hmnf.chm, whitepine %in% 1)
aspen <- subset(hmnf.chm, aspen %in% 1)
oak <- subset(hmnf.chm, oak %in% 1)
maple <- subset(hmnf.chm, maple %in% 1)






jackpine.model <- lm(chm.max~
                       solar+
                       age+twi+tpi+bt+tgs+
                       ppt+T150_AWC+T50_sand+T150_sand+T50_clay+T150_clay+T50_OM+
                       T150_OM+T50_pH+Water_Table+wet
                     , data=jackpine, weights = age)

redpine.model <- lm(chm.max~
                      solar+
                      age+tpi+bt+tgs+
                      ppt+T150_AWC+T50_sand+T150_sand+T50_clay+T150_clay+T50_OM+
                      T150_OM+T50_pH+Water_Table+wet
                     , data=redpine, weights = age)

whitepine.model <- lm(chm.max~
                        solar+
                        age+twi+tpi+bt+tgs+
                        ppt+T150_AWC+T50_sand+T150_sand+T50_clay+T150_clay+T50_OM+
                        T150_OM+T50_pH+Water_Table+wet
                     , data=whitepine, weights = age)

aspen.model <- lm(chm.max~
                    solar+
                    age+twi+tpi+bt+tgs+
                    ppt+T150_AWC+T50_sand+T150_sand+T50_clay+T150_clay+T50_OM+
                    T150_OM+T50_pH+Water_Table+wet
                     , data=aspen, weights = age)

oak.model <- lm(chm.max~
                  solar+
                  age+twi+tpi+bt+tgs+
                  ppt+T150_AWC+T50_sand+T150_sand+T50_clay+T150_clay+T50_OM+
                  T150_OM+T50_pH+Water_Table+wet
                     , data=oak, weights = age)

maple.model <- lm(chm.max~
                    solar+
                    age+tpi+bt+tgs+
                    ppt+T150_AWC+T50_sand+T150_sand+T50_clay+T150_clay+T50_OM+
                    T150_OM+T50_pH+Water_Table+wet
                     , data=maple, weights = age)

stepAIC(jackpine.model)
stepAIC(redpine.model)
stepAIC(whitepine.model)
stepAIC(aspen.model)
stepAIC(oak.model)
stepAIC(maple.model)



summary(jackpine.model)
summary(redpine.model)
summary(whitepine.model)
summary(aspen.model)
summary(oak.model)
summary(maple.model)

hmnf.chm$position <- (hmnf.chm$toi - median(hmnf.chm$toi))/(median(hmnf.chm$toi) - min(hmnf.chm$toi))
hmnf.chm$sunslope <- (hmnf.chm$solar - median(hmnf.chm$solar))/(median(hmnf.chm$solar) - min(hmnf.chm$solar))
hmnf.chm$elevation <- hmnf.chm$dem

soilsubset <- subset(hmnf.chm, T150_AWC < 15 & spodosols < 0.5 & Water_Table > 150 & age > 5 & chm > 1)
#soilsubset <- subset(soilsubset, slope.500 >= 0.05 & slope >= atan(0.0))
#soilsubset <- subset(soilsubset, slope.500 >= 0.05 & slope >= atan(0.0))
soilsubset <- hmnf.chm
jackpine <- subset(soilsubset, jackpine %in% 1)
redpine <- subset(soilsubset, redpine %in% 1)
whitepine <- subset(soilsubset, whitepine %in% 1)
aspen <- subset(soilsubset, aspen %in% 1)
oak <- subset(soilsubset, oak %in% 1)
maple <- subset(soilsubset, maple %in% 1)


ggplot()+
  geom_smooth(aes(x=sunslope, y=chm.max, weight = age, col='maple'), data = maple)+
  geom_smooth(aes(x=sunslope, y=chm.max, weight = age, col='oak'), data = oak)+
  geom_smooth(aes(x=sunslope, y=chm.max, weight = age, col='aspen'), data = aspen)+
  geom_smooth(aes(x=sunslope, y=chm.max, weight = age, col='white pine'), data = whitepine)+
  geom_smooth(aes(x=sunslope, y=chm.max, weight = age, col='red pine'), data = redpine)+
  geom_smooth(aes(x=sunslope, y=chm.max, weight = age, col='jack pine'), data = jackpine)+
  scale_color_manual(name='forest type',  values = c('jack pine'='green','red pine'='darkgreen','white pine'='darkcyan',
                                                     'aspen'='khaki4','oak'='orange','maple'='red'))
ggplot()+
  geom_smooth(aes(x=position, y=chm.max, weight = age, col='maple'), data = maple)+
  geom_smooth(aes(x=position, y=chm.max, weight = age, col='oak'), data = oak)+
  geom_smooth(aes(x=position, y=chm.max, weight = age, col='aspen'), data = aspen)+
  geom_smooth(aes(x=position, y=chm.max, weight = age, col='white pine'), data = whitepine)+
  geom_smooth(aes(x=position, y=chm.max, weight = age, col='red pine'), data = redpine)+
  geom_smooth(aes(x=position, y=chm.max, weight = age, col='jack pine'), data = jackpine)+
  scale_color_manual(name='forest type',  values = c('jack pine'='green','red pine'='darkgreen','white pine'='darkcyan',
                                                     'aspen'='khaki4','oak'='orange','maple'='red'))

ggplot()+
  geom_smooth(aes(x=Water_Table, y=chm.max, weight = age, col='maple'), data = maple)+
  geom_smooth(aes(x=Water_Table, y=chm.max, weight = age, col='oak'), data = oak)+
  geom_smooth(aes(x=Water_Table, y=chm.max, weight = age, col='aspen'), data = aspen)+
  geom_smooth(aes(x=Water_Table, y=chm.max, weight = age, col='white pine'), data = whitepine)+
  geom_smooth(aes(x=Water_Table, y=chm.max, weight = age, col='red pine'), data = redpine)+
  geom_smooth(aes(x=Water_Table, y=chm.max, weight = age, col='jack pine'), data = jackpine)+
  scale_color_manual(name='forest type',  values = c('jack pine'='green','red pine'='darkgreen','white pine'='darkcyan',
                                                     'aspen'='khaki4','oak'='orange','maple'='red'))


ggplot()+
  geom_smooth(aes(x=age, y=chm.max, col='maple'), data = maple)+
  geom_smooth(aes(x=age, y=chm.max, col='oak'), data = oak)+
  geom_smooth(aes(x=age, y=chm.max, col='aspen'), data = aspen)+
  geom_smooth(aes(x=age, y=chm.max, col='white pine'), data = whitepine)+
  geom_smooth(aes(x=age, y=chm.max, col='red pine'), data = redpine)+
  geom_smooth(aes(x=age, y=chm.max, col='jack pine'), data = jackpine)+
  scale_color_manual(name='forest type',  values = c('jack pine'='green','red pine'='darkgreen','white pine'='darkcyan',
                                                     'aspen'='khaki4','oak'='orange','maple'='red'))

exposedridge <- subset(oak, position >=.50 & sunslope >=.50)
medianslope <- subset(oak,  position >= -.50 & position < .50 & sunslope >= -.50 & sunslope < .50)
shadycove <- subset(oak, position < -.50 & sunslope < -.50) 
cove <- subset(oak, position < -.50 )#& sunslope >= -.50 & sunslope < .50) 
ridge <- subset(oak, position >= .50 )# & sunslope >= -.50 & sunslope < .50) 
shadyslope <- subset(oak, position >= -.50 & position < .50 & sunslope < -.50) 
sunnyslope <- subset(oak, position >= -.50 & position < .50 & sunslope >= .50)

ggplot()+
  #geom_smooth(aes(x=age, y=chm, col='sunny ridge'), data= exposedridge)+
  geom_smooth(aes(x=age, y=chm, col='ridge'), data= ridge)+
  #geom_smooth(aes(x=age, y=chm, col='sunny slope'), data= sunnyslope)+
  geom_smooth(aes(x=age, y=chm, col='median slope'), data= medianslope)+
  geom_smooth(aes(x=age, y=chm, col='shady slope'), data= shadyslope)+
  geom_smooth(aes(x=age, y=chm, col='cove'), data= cove)+
  geom_smooth(aes(x=age, y=chm, col='shady cove'), data= shadycove)+
  geom_smooth(aes(x=age, y=chm, col='all slopes'), data= oak)+
  scale_color_manual(name='site',  values = c('all slopes'='black','sunny ridge'='red','ridge'='orange','sunny slope'='yellow',
                                              'median slope'='green',
                                              'shady slope'='darkgreen','cove'='cyan','shady cove'='blue'
  ))+
  scale_y_continuous(name='Canopy Height (m)', breaks=c((0:10)*5))+
  scale_x_continuous(name='age (years)', breaks=c((0:50)*5))+
  coord_fixed(ratio=2/1,ylim=c(0,40))+
  labs(title = 'Oak')


exposedridge <- subset(whitepine, position >=.50 & sunslope >=.50)
medianslope <- subset(whitepine,  position >= -.50 & position < .50 & sunslope >= -.50 & sunslope < .50)
shadycove <- subset(whitepine, position < -.50 & sunslope < -.50) 
cove <- subset(whitepine, position < -.50  )#& sunslope >= -.50 & sunslope < .50) 
ridge <- subset(whitepine, position >= .50  )#& sunslope >= -.50 & sunslope < .50) 
shadyslope <- subset(whitepine, position >= -.50 & position < .50 & sunslope < -.50) 
sunnyslope <- subset(whitepine, position >= -.50 & position < .50 & sunslope >= .50)

ggplot()+
  #geom_smooth(aes(x=age, y=chm, col='sunny ridge'), data= exposedridge)+
  #geom_smooth(aes(x=age, y=chm, col='ridge'), data= ridge)+
  #geom_smooth(aes(x=age, y=chm, col='sunny slope'), data= sunnyslope)+
  geom_smooth(aes(x=age, y=chm, col='median slope'), data= medianslope)+
  geom_smooth(aes(x=age, y=chm, col='shady slope'), data= shadyslope)+
  geom_smooth(aes(x=age, y=chm, col='cove'), data= cove)+
  geom_smooth(aes(x=age, y=chm, col='shady cove'), data= shadycove)+
  geom_smooth(aes(x=age, y=chm, col='all slopes'), data= whitepine)+
  scale_color_manual(name='site',  values = c('all slope'='black','sunny ridge'='red','ridge'='orange','sunny slope'='yellow',
                                              'median slope'='green',
                                              'shady slope'='darkgreen','cove'='cyan','shady cove'='blue'
  ))+
  scale_y_continuous(name='Canopy Height (m)', breaks=c((0:10)*5))+
  scale_x_continuous(name='age (years)', breaks=c((0:50)*5))+
  coord_fixed(ratio=2/1,ylim=c(0,40))+
  labs(title = 'white pine')


exposedridge <- subset(maple, position >=.50 & sunslope >=.50)
medianslope <- subset(maple,  position >= -.50 & position < .50 & sunslope >= -.50 & sunslope < .50)
shadycove <- subset(maple, position < -.50 & sunslope < -.50) 
cove <- subset(maple, position < -.50  )#& sunslope >= -.50 & sunslope < .50) 
ridge <- subset(maple, position >= .50  )#& sunslope >= -.50 & sunslope < .50) 
shadyslope <- subset(maple, position >= -.50 & position < .50 & sunslope < -.50) 
sunnyslope <- subset(maple, position >= -.50 & position < .50 & sunslope >= .50)

ggplot()+
  #geom_smooth(aes(x=age, y=chm, col='sunny ridge'), data= exposedridge)+
  geom_smooth(aes(x=age, y=chm, col='ridge'), data= ridge)+
  #geom_smooth(aes(x=age, y=chm, col='sunny slope'), data= sunnyslope)+
  geom_smooth(aes(x=age, y=chm, col='median slope'), data= medianslope)+
  geom_smooth(aes(x=age, y=chm, col='shady slope'), data= shadyslope)+
  geom_smooth(aes(x=age, y=chm, col='cove'), data= cove)+
  #geom_smooth(aes(x=age, y=chm, col='shady cove'), data= shadycove)+
  geom_smooth(aes(x=age, y=chm, col='all slopes'), data= maple)+
  scale_color_manual(name='site',  values = c('all slope'='black','sunny ridge'='red','ridge'='orange','sunny slope'='yellow',
                                              'median slope'='green',
                                              'shady slope'='darkgreen','cove'='cyan','shady cove'='blue'
  ))+
  scale_y_continuous(name='Canopy Height (m)', breaks=c((0:10)*5))+
  scale_x_continuous(name='age (years)', breaks=c((0:50)*5))+
  coord_fixed(ratio=2/1,ylim=c(0,40))+
  labs(title = 'maple')

exposedridge <- subset(aspen, position >=.50 & sunslope >=.50)
medianslope <- subset(aspen,  position >= -.50 & position < .50 & sunslope >= -.50 & sunslope < .50)
shadycove <- subset(aspen, position < -.50 & sunslope < -.50) 
cove <- subset(aspen, position < -.50  )#& sunslope >= -.50 & sunslope < .50) 
ridge <- subset(aspen, position >= .50  )#& sunslope >= -.50 & sunslope < .50) 
shadyslope <- subset(aspen, position >= -.50 & position < .50 & sunslope < -.50) 
sunnyslope <- subset(aspen, position >= -.50 & position < .50 & sunslope >= .50)

ggplot()+
  #geom_smooth(aes(x=age, y=chm, col='sunny ridge'), data= exposedridge)+
  geom_smooth(aes(x=age, y=chm, col='ridge'), data= ridge)+
  #geom_smooth(aes(x=age, y=chm, col='sunny slope'), data= sunnyslope)+
  geom_smooth(aes(x=age, y=chm, col='median slope'), data= medianslope)+
  geom_smooth(aes(x=age, y=chm, col='shady slope'), data= shadyslope)+
  geom_smooth(aes(x=age, y=chm, col='cove'), data= cove)+
  #geom_smooth(aes(x=age, y=chm, col='shady cove'), data= shadycove)+
  geom_smooth(aes(x=age, y=chm, col='all slopes'), data= aspen)+
  scale_color_manual(name='site',  values = c('all slope'='black','sunny ridge'='red','ridge'='orange','sunny slope'='yellow',
                                              'median slope'='green',
                                              'shady slope'='darkgreen','cove'='cyan','shady cove'='blue'
  ))+
  scale_y_continuous(name='Canopy Height (m)', breaks=c((0:10)*5))+
  scale_x_continuous(name='age (years)', breaks=c((0:50)*5))+
  coord_fixed(ratio=2/1,ylim=c(0,40))+
  labs(title = 'Aspen')


spodic <- subset(hmnf.chm, T150_AWC < 15 & spodosols > 0.5 & Water_Table > 150 & age > 5 & chm > 1 & aspen %in% 1)
entic <- subset(hmnf.chm, T150_AWC < 15 & spodosols < 0.5 & Water_Table > 150 & age > 5 & chm > 1 & aspen %in% 1)
mwdsand <- subset(hmnf.chm, T150_AWC < 15 & spodosols > -1 & Water_Table > 75 & Water_Table < 150 & age > 5 & chm > 1 & aspen %in% 1)
spdsand <- subset(hmnf.chm, T150_AWC < 15 & spodosols > -1 & Water_Table > 25 & Water_Table < 75 & age > 5 & chm > 1 & aspen %in% 1)
loamy <- subset(hmnf.chm, T150_AWC >= 15 & spodosols > -1 & Water_Table > 25 & age > 5 & chm > 1 & aspen %in% 1)
acidwet <- subset(hmnf.chm,  (spodosols > 0.5 | T50_pH < 5.5) & Water_Table < 25 & age > 5 & chm > 1 & aspen %in% 1)
basicwet <- subset(hmnf.chm, (spodosols < 0.5 & T50_pH >= 5.5) & Water_Table > 150 & age > 5 & chm > 1 & aspen %in% 1)
aspen <- subset(hmnf.chm, age > 5 & chm > 1 & aspen %in% 1)

ggplot()+
  geom_smooth(aes(x=age, y=chm, col='spodic'), data= spodic)+
  geom_smooth(aes(x=age, y=chm, col='entic'), data= entic)+
  geom_smooth(aes(x=age, y=chm, col='mwdsand'), data= mwdsand)+
  geom_smooth(aes(x=age, y=chm, col='spdsand'), data= spdsand)+
  geom_smooth(aes(x=age, y=chm, col='loamy'), data= loamy)+
  geom_smooth(aes(x=age, y=chm, col='acidwet'), data= acidwet)+
  geom_smooth(aes(x=age, y=chm, col='basicwet'), data= basicwet)+
  geom_smooth(aes(x=age, y=chm, col='all soils'), data= aspen)+
  scale_color_manual(name='site',  values = c('all soils'='black','entic'='red','spodic'='orange','mwdsand'='yellow',
                                              'spdsand'='green',
                                              'loamy'='darkgreen','acidwet'='cyan','basicwet'='blue'
  ))+
  scale_y_continuous(name='Canopy Height (m)', breaks=c((0:10)*5))+
  scale_x_continuous(name='age (years)', breaks=c((0:50)*5))+
  coord_fixed(ratio=2/1,ylim=c(0,40))+
  labs(title = 'Aspen')

spodic <- subset(hmnf.chm, T150_AWC < 15 & spodosols > 0.5 & Water_Table > 150 & age > 5 & chm > 1 & oak %in% 1)
entic <- subset(hmnf.chm, T150_AWC < 15 & spodosols < 0.5 & Water_Table > 150 & age > 5 & chm > 1 & oak %in% 1)
mwdsand <- subset(hmnf.chm, T150_AWC < 15 & spodosols > -1 & Water_Table > 75 & Water_Table < 150 & age > 5 & chm > 1 & oak %in% 1)
spdsand <- subset(hmnf.chm, T150_AWC < 15 & spodosols > -1 & Water_Table > 25 & Water_Table < 75 & age > 5 & chm > 1 & oak %in% 1)
loamy <- subset(hmnf.chm, T150_AWC >= 15 & spodosols > -1 & Water_Table > 25 & age > 5 & chm > 1 & oak %in% 1)
acidwet <- subset(hmnf.chm,  (spodosols > 0.5 | T50_pH < 5.5) & Water_Table < 25 & age > 5 & chm > 1 & oak %in% 1)
basicwet <- subset(hmnf.chm, (spodosols < 0.5 & T50_pH >= 5.5) & Water_Table > 150 & age > 5 & chm > 1 & oak %in% 1)
oak <- subset(hmnf.chm, age > 5 & chm > 1 & oak %in% 1)

ggplot()+
  geom_smooth(aes(x=age, y=chm, col='spodic'), data= spodic)+
  geom_smooth(aes(x=age, y=chm, col='entic'), data= entic)+
  geom_smooth(aes(x=age, y=chm, col='mwdsand'), data= mwdsand)+
  geom_smooth(aes(x=age, y=chm, col='spdsand'), data= spdsand)+
  geom_smooth(aes(x=age, y=chm, col='loamy'), data= loamy)+
  geom_smooth(aes(x=age, y=chm, col='acidwet'), data= acidwet)+
  geom_smooth(aes(x=age, y=chm, col='basicwet'), data= basicwet)+
  geom_smooth(aes(x=age, y=chm, col='all soils'), data= oak)+
  scale_color_manual(name='site',  values = c('all soils'='black','entic'='red','spodic'='orange','mwdsand'='yellow',
                                              'spdsand'='green',
                                              'loamy'='darkgreen','acidwet'='cyan','basicwet'='blue'
  ))+
  scale_y_continuous(name='Canopy Height (m)', breaks=c((0:10)*5))+
  scale_x_continuous(name='age (years)', breaks=c((0:50)*5))+
  coord_fixed(ratio=2/1,ylim=c(0,40))+
  labs(title = 'Oak')

#randomforest----
#
hmnf.chm2 <- hmnf.chm
hmnf.chm2$age <- 0
hmnf.chm2$chm.max <- 0
hmnf.chm2 <- rbind(hmnf.chm, hmnf.chm2)


x1 <- 0:500
pred1 <-  predict(mod1, list(age = x1))
pred2 <-  predict(mod2, list(age = x1))
pred3 <-  predict(mod3, list(age = x1))
pred4 <-  predict(mod4, list(age = x1))
pred5 <-  predict(mod5, list(age = x1))

rpmod <- function(age){
b0= 2.0434
b1= 0.9978
b2= -0.0147
b3= 1.0937
b4= -0.0035
S= 58#subset(redpine,SITE_INDEX_SPP %in%  'PIRE') $ SITE_INDEX %>% median()

(0+b0*(S^b1)*(1-exp(b2*age))^(b3*S^b4)   )*0.3048}

asmod <- function(age){
b0= 5.2188
b1= 0.6855
b2= -0.0301
b3= 50.0071
b4= -0.8695
S= 63#subset(aspen,SITE_INDEX_SPP %in%  c('POGR4','POTR5')) $ SITE_INDEX %>% median()
(0+b0*(S^b1)*(1-exp(b2*age))^(b3*S^b4)   )*0.3048}


# treecurve <- function(age, b1,b2,b3){
#  
#   b1*(1-exp(b2*age))^b3 #essentially equivalent to Chapman-Richards Model with additional factors to connect to site index
# }
# 
# mod1 <- nlsLM(chm.max ~ treecurve(age, b1,b2,b3), start = list(
#   b1= 27.2291 ,b2= -0.0301  ,b3= 1.3630
# ), data = aspen)
# fake4 <- predict(mod1, list(age=x1))


hmnf.chm2$pred1 <-  asmod(hmnf.chm2$age)
hmnf.chm2$pred2 <-  rpmod(age = hmnf.chm2$age)

mod <- lm(chm.max~
               solar+
               hilly+
               age+
               pred1+
               pred2+
               toip+
               toin+
               bt+
               tgs+
               ppt+
               T150_AWC+
               T50_sand+
               T150_sand+
               T50_clay+
               T150_clay+
               T50_OM+
               T150_OM+
               T50_pH+
               Water_Table+
               wet+
               spodic+
               carbdepth+
               Bhs+
               jackpine+
               redpine+
               whitepine+
               aspen+
               oak+
               maple
             , data=hmnf.chm2)
summary(mod)

hmnf.chm2$wts <-  ifelse(hmnf.chm2$redpine %in% 1, 1, 1)
rf <- ranger(chm.max~
           solar+
           hilly+
           #age+
           pred1+
           pred2+
           toip+
           toin+
           bt+
           tgs+
           ppt+
           T150_AWC+
           T50_sand+
           T150_sand+
           T50_clay+
           T150_clay+
           T50_OM+
           T150_OM+
           T50_pH+
           Water_Table+
           wet+
           spodic+
           carbdepth+
           Bhs+
           jackpine+
           redpine+
           whitepine+
           aspen+
           oak+
           maple
         , data=hmnf.chm2,
         num.trees=500, sample.fraction =0.002, always.split.variables = c('pred1'),
           write.forest = TRUE, importance = 'impurity', case.weights= hmnf.chm2$wts
         )
importance(rf)


hmnf.chm.0 <- hmnf.chm2 %>% mutate(age=0, pred1 = 0, pred2 = 0, chm=0, aspen=0, oak=0, maple=0, jackpine=0, redpine=1, whitepine=0)
hmnf.chm.10 <- hmnf.chm2 %>% mutate(age=10, pred1 = asmod(10), pred2 = rpmod(10), aspen=0, oak=0, maple=0, jackpine=0, redpine=1, whitepine=0)
hmnf.chm.25 <- hmnf.chm2 %>% mutate(age=25, pred1 = asmod(25), pred2 = rpmod(25), aspen=0, oak=0, maple=0, jackpine=0, redpine=1, whitepine=0)
hmnf.chm.50 <- hmnf.chm2 %>% mutate(age=50, pred1 = asmod(50), pred2 = rpmod(50), aspen=0, oak=0, maple=0, jackpine=0, redpine=1, whitepine=0)
hmnf.chm.75 <- hmnf.chm2 %>% mutate(age=75, pred1 = asmod(75), pred2 = rpmod(75), aspen=0, oak=0, maple=0, jackpine=0, redpine=1, whitepine=0)
hmnf.chm.100 <- hmnf.chm2 %>% mutate(age=100, pred1 = asmod(100), pred2 = rpmod(100), aspen=0, oak=0, maple=0, jackpine=0, redpine=1, whitepine=0)
hmnf.chm.150 <- hmnf.chm2 %>% mutate(age=150, pred1 = asmod(150), pred2 = rpmod(150), aspen=0, oak=0, maple=0, jackpine=0, redpine=1, whitepine=0)
hmnf.chm.200 <- hmnf.chm2 %>% mutate(age=200, pred1 = asmod(200), pred2 = rpmod(200), aspen=0, oak=0, maple=0, jackpine=0, redpine=1, whitepine=0)
hmnf.chm.10 <- hmnf.chm.10 %>% mutate(chm = predictions(predict(rf, data=hmnf.chm.10)))
hmnf.chm.25 <- hmnf.chm.25 %>% mutate(chm = predictions(predict(rf, data=hmnf.chm.25)))
hmnf.chm.50 <- hmnf.chm.50 %>% mutate(chm = predictions(predict(rf, data=hmnf.chm.50)))
hmnf.chm.75 <- hmnf.chm.75 %>% mutate(chm = predictions(predict(rf, data=hmnf.chm.75)))
hmnf.chm.100 <- hmnf.chm.100 %>% mutate(chm = predictions(predict(rf, data=hmnf.chm.100)))
hmnf.chm.150 <- hmnf.chm.150 %>% mutate(chm = predictions(predict(rf, data=hmnf.chm.150)))
hmnf.chm.200 <- hmnf.chm.200 %>% mutate(chm = predictions(predict(rf, data=hmnf.chm.200)))

x1<-c(0,10,25,50,75,100,150,200)
aspen = asmod(x1)
redpine = rpmod(x1)

hmnf.chm.alltime <-  rbind(hmnf.chm.0,hmnf.chm.10,hmnf.chm.25, hmnf.chm.50,hmnf.chm.75,hmnf.chm.100,hmnf.chm.150,hmnf.chm.200)
ggplot()+
  geom_boxplot(aes(x=as.factor(age),y=chm), data=hmnf.chm.alltime)+
  geom_point(aes(x=as.factor(x1),y=aspen,col='aspen'))+
  geom_point(aes(x=as.factor(x1),y=redpine,col='red pine'))+
  scale_color_manual(name='curve',  values = c('aspen'='blue','red pine'='red'))+
  scale_y_continuous(name='Canopy Height (m)', breaks=c((0:10)*5))+
  scale_x_discrete(name='Age')

### 90 m predicted SI  ----
### 
solar <- rast('D:/GIS/DEM/hmnfsolarmean.tif'); names(solar) = 'solar'
toip <- rast('D:/GIS/DEM/hmnf-pos-openess.tif'); names(toip) = 'toip'
toin <- rast('D:/GIS/DEM/hmnf-neg-openess.tif'); names(toin) = 'toin'
slope <- rast('D:/GIS/DEM/hmnfslope.tif'); names(slope) = 'slope'#radians
dem <- rast('D:/GIS/DEM/hmnfdem30.tif'); names(dem) = 'dem'
bt <- rast('D:/scripts/snow/output/bt.90.alt.tif'); names(bt) = 'bt'
tgs <- rast('D:/scripts/snow/output/tgs.90.alt.tif'); names(tgs) = 'tgs'
ppt <- rast('C:/a/Ecological_Sites/GIS/Climate/PRISM2010/P/p0112/w001001.adf'); names(ppt) = 'ppt'
p1 <- rast('C:/a/Ecological_Sites/GIS/Climate/PRISM2010/seasons/p1/w001001.adf'); names(p1) = 'p1'
p2 <- rast('C:/a/Ecological_Sites/GIS/Climate/PRISM2010/seasons/p2/w001001.adf'); names(p2) = 'p2'
p3 <- rast('C:/a/Ecological_Sites/GIS/Climate/PRISM2010/seasons/p3/w001001.adf'); names(p3) = 'p3'
p4 <- rast('C:/a/Ecological_Sites/GIS/Climate/PRISM2010/seasons/p4/w001001.adf'); names(p4) = 'p4'

slope.500 <- rast('D:/GIS/DEM/slope.500.tif'); names(slope.500) = 'slope.500'

mukey <- rast('D:/GIS/SOIL/2021/northeast_gssurgo90m.tif'); names(mukey) = 'mukey'

solar <- aggregate(solar, fact=3, fun='mean')
toip <- aggregate(toip, fact=3, fun='mean')
toin <- aggregate(toin, fact=3, fun='mean')
slope <- aggregate(slope, fact=3, fun='mean')
dem <- aggregate(dem, fact=3, fun='mean')
slope.500 <- aggregate(slope.500, fact=3, fun='mean')
hilly <- (slope.500 > 0.05)*1
hilly <-  focal(hilly, w=focalMat(hilly, 500, type=c('circle')))
hilly <-  focal(hilly, w=focalMat(hilly, 500, type=c('circle')))
names(hilly) <- 'hilly'
bt <- project(bt, dem, method='bilinear')
tgs <- project(tgs, dem, method='bilinear')
ppt <- project(ppt, dem, method='bilinear')
p1 <- project(p1, dem, method='bilinear')
p2 <- project(p2, dem, method='bilinear')
p3 <- project(p3, dem, method='bilinear')
p4 <- project(p4, dem, method='bilinear')
mukey <- project(mukey, dem, method='near')

brk <- c(mukey, solar, toip, toin, dem, slope, hilly, bt, tgs, ppt, p1,p2,p3,p4)
brk.pts <- as.data.frame(brk, xy=TRUE)

#soils processing
s <- read.csv('C:/workspace/SoilSorts/fy2021-refresh/s.cleaned.csv')
s <- subset(s, majcompflag %in% 'TRUE' & !is.na(T150_AWC)&!is.na(T50_sand)&!is.na(T150_sand)&!is.na(T50_clay)&!is.na(T150_clay)&!is.na(T50_OM)&!is.na(T150_OM)&!is.na(T50_pH)&!is.na(Water_Table)&!is.na(carbdepth)&!is.na(Bhs))
colnames(s)

s$spodosols <- ifelse(s$taxorder %in% 'spodosols', 1,0)
s$spodic <- ifelse(grepl('spodic',s$taxsubgrp)|s$taxorder %in% 'spodosols', 1,0)
s$Bhs <- ifelse(s$Bhs %in% c('yes', 1), 1,0)
s$spodosols.up <- s$spodosols*s$Water_Table>100
s$spodic.up <- s$spodic*s$Water_Table
s$Bhs.up <- s$Bhs*s$Water_Table
#---

s.mu <- s %>% group_by(lmapunitiid=lmapunitiid) %>% dplyr::summarise(
  T150_AWC = wtd.mean(T150_AWC, weights=comppct_r, na.rm=T),
  T50_sand = wtd.mean(T50_sand, weights=comppct_r, na.rm=T),
  T150_sand = wtd.mean(T150_sand, weights=comppct_r, na.rm=T),
  T50_clay = wtd.mean(T50_clay, weights=comppct_r, na.rm=T),
  T150_clay = wtd.mean(T150_clay, weights=comppct_r, na.rm=T),
  T50_OM = wtd.mean(T50_OM, weights=comppct_r, na.rm=T),
  T150_OM = wtd.mean(T150_OM, weights=comppct_r, na.rm=T),
  T50_pH = wtd.mean(T50_pH, weights=comppct_r, na.rm=T),
  Water_Table = wtd.mean(Water_Table, weights=comppct_r, na.rm=T),
  spodic = wtd.mean(spodic, weights=comppct_r,  na.rm=T),
  spodosols = wtd.mean(spodosols, weights=comppct_r,  na.rm=T),
  Bhs = wtd.mean(Bhs, weights=comppct_r,  na.rm=T),
  carbdepth = wtd.mean(carbdepth, weights=comppct_r,  na.rm=T),
  humicdepth = wtd.mean(humicdepth, weights=comppct_r, na.rm=T)
)
s.mu$wet <- 25/(s.mu$Water_Table+25)
s.mu$hydric <- pmin(1,(pmax(50-s.mu$Water_Table, 0)/25))
s.mu$moist <- pmin(1,(pmax(200-s.mu$Water_Table, 0)/100))




brk.pts <- left_join(brk.pts, s.mu, by=c('mukey'='lmapunitiid'))
brk.pts <- subset(brk.pts, !is.na(T150_AWC)&!is.na(T50_sand)&!is.na(T150_sand)&!is.na(T50_clay)&!is.na(T150_clay)&!is.na(T50_OM)&!is.na(T150_OM)&!is.na(T50_pH)&!is.na(Water_Table)&!is.na(wet)&!is.na(spodic)&!is.na(carbdepth)&!is.na(Bhs))

saveRDS(brk.pts, 'output/brk.pts.RDS')
saveRDS(hmnf.chm2, 'output/hmnf.chm3.RDS')



#model predictions ----
#
library(sf)
library(terra)
library(raster)
library(minpack.lm)
library(growthmodels)
library(dplyr)
library(ggplot2)
library(Hmisc)
library(MASS)
library(ranger)

library(partykit)
#partykit::lmtree()

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
rpmod <- function(age){
  b0= 2.0434
  b1= 0.9978
  b2= -0.0147
  b3= 1.0937
  b4= -0.0035
  S= 58#subset(redpine,SITE_INDEX_SPP %in%  'PIRE') $ SITE_INDEX %>% median()
  
  (0+b0*(S^b1)*(1-exp(b2*age))^(b3*S^b4)   )*0.3048}

asmod <- function(age){
  b0= 5.2188
  b1= 0.6855
  b2= -0.0301
  b3= 50.0071
  b4= -0.8695
  S= 63#subset(aspen,SITE_INDEX_SPP %in%  c('POGR4','POTR5')) $ SITE_INDEX %>% median()
  (0+b0*(S^b1)*(1-exp(b2*age))^(b3*S^b4)   )*0.3048}

brk.pts <- readRDS('output/brk.pts.RDS')
hmnf.chm2 <- readRDS('output/hmnf.chm3.RDS')
dem <- rast('D:/GIS/DEM/hmnfdem30.tif'); names(dem) = 'dem'
brk.pts$toi <- brk.pts$toip-brk.pts$toin
hmnf.chm2$toi <- hmnf.chm2$toip-hmnf.chm2$toin

1/mean(hmnf.chm2$jackpine)/2
1/mean(hmnf.chm2$redpine)/2
1/mean(hmnf.chm2$whitepine)/2
1/mean(hmnf.chm2$aspen)/2
1/mean(hmnf.chm2$oak)/2
1/mean(hmnf.chm2$maple)/2

#age+taxon+soil+topo+clim
#texture+drainage
#temp+ppt
#
#generic linear model
# hmnf.chm2$wts <-   ((hmnf.chm2$hilly >= 0.5)*3+1)
# 
# linear.mod <- lm(chm.max~
#                    pred1+
#                    pred1:solar+
#                    pred1:toip+
#                    # pred1:toin+
#                    # pred1:bt+
#                    pred1:tgs+
#                    pred1:ppt+
#                    # pred1:p1+
#                    # pred1:p2+
#                    # pred1:p3+
#                    # pred1:p4+
#                    pred1:T150_AWC+
#                    pred1:T50_sand+
#                    pred1:T150_sand+
#                    pred1:T50_clay+
#                    pred1:T150_clay+
#                    pred1:T50_OM+
#                    pred1:T150_OM+
#                    pred1:T50_pH+
#                    pred1:Water_Table+
#                    pred1:wet+
#                    pred1:hydric+
#                    pred1:moist+
#                    pred1:spodosols+
#                    pred1:Bhs+
#                    pred1:spodosols:Water_Table+
#                    # pred1:Bhs:Water_Table+
#                    
#                    jackpine:pred1+
#                    redpine:pred1+
#                    whitepine:pred1+
#                    aspen:pred1+
#                    oak:pred1+
#                    maple:pred1#+
#                  
#                  # jackpine:pred1:hydric+
#                  # redpine:pred1:hydric+
#                  # whitepine:pred1:hydric+
#                  # aspen:pred1:hydric+
#                  # oak:pred1:hydric+
#                  # maple:pred1:hydric+
#                  # 
#                  # jackpine:pred1:T50_pH+
#                  # redpine:pred1:T50_pH+
#                  # whitepine:pred1:T50_pH+
#                  # aspen:pred1:T50_pH+
#                  # oak:pred1:T50_pH+
#                  # maple:pred1:T50_pH
#                  , data=hmnf.chm2)
# 
# #stepAIC(linear.mod)
# summary(linear.mod)

#step 1 rf model 
#step 2 substitute with median climate run model calculate residuals
#step 3 lm model with climate run model calculate residuals
#step 4 final rf model
# alternatively run lm model first?
# 
# 
 


lmod <- function(train){
  linear.mod <- lm(chm.max~
                     pred1+
                     pred1:solar+
                     pred1:toip+
                     pred1:tgs+
                     pred1:ppt+
                     pred1:T150_AWC+
                     pred1:T50_sand+
                     pred1:T150_sand+
                     pred1:T50_clay+
                     pred1:T150_clay+
                     pred1:T50_OM+
                     pred1:T150_OM+
                     pred1:T50_pH+
                     pred1:Water_Table+
                     pred1:wet+
                     pred1:spodosols+
                     pred1:Bhs+
                     pred1:spodosols:Water_Table+
                     
                     jackpine:pred1+
                     redpine:pred1+
                     whitepine:pred1+
                     aspen:pred1+
                     oak:pred1+
                     maple:pred1
                   , data=train, weights = train$wts)
  
  summary(linear.mod)
    return(linear.mod)
}

rf.model <- function(train){
  splitwts <- c(
    solar=1/5/3,
    #hilly=1/5,
    #age=1/5,
    pred1=1/5/2,
    pred2=1/5/2,
    toip=1/5/3,
    toin=1/5/3,
    bt=1/5/2/2,
    tgs=1/5/2/2,
    ppt=1/5/2/5,
    p1=1/5/2/5,
    p2=1/5/2/5,
    p3=1/5/2/5,
    p4=1/5/2/5,
    T150_AWC=1/5/2/11,
    T50_sand=1/5/2/11,
    T150_sand=1/5/2/11,
    T50_clay=1/5/2/11,
    T150_clay=1/5/2/11,
    T50_OM=1/5/2/11,
    T150_OM=1/5/2/11,
    T50_pH=1/5/2/11,
    Water_Table=1/5/2/2,
    wet=1/5/2/2,
    spodic=1/5/2/11,
    carbdepth=1/5/2/11,
    Bhs=1/5/2/11,
    jackpine=1/5/6,
    redpine=1/5/6,
    whitepine=1/5/6,
    aspen=1/5/6,
    oak=1/5/6,
    maple=1/5/6
  )
  rf <- ranger(chm.max~
                 solar+
                 #hilly+
                 #age+
                 pred1+
                 pred2+
                 toip+
                 toin+
                 bt+
                 tgs+
                 ppt+
                 p1+
                 p2+
                 p3+
                 p4+
                 T150_AWC+
                 T50_sand+
                 T150_sand+
                 T50_clay+
                 T150_clay+
                 T50_OM+
                 T150_OM+
                 T50_pH+
                 Water_Table+
                 wet+
                 spodic+
                 carbdepth+
                 Bhs+
                 jackpine+
                 redpine+
                 whitepine+
                 aspen+
                 oak+
                 maple
               , data=train,
               num.trees=500, sample.fraction =.2,
               
               split.select.weights = splitwts,
               # max.depth = 5,  
               # mtry=7,
               write.forest = TRUE, importance = 'impurity', case.weights= train$wts
  )
  return(rf)
}



hybrid.model <- function(train){
  splitwts <- c(
    solar=1/5/3,
    #hilly=1/5,
    #age=1/5,
    pred1=1/5/2,
    pred2=1/5/2,
    toip=1/5/3,
    toin=1/5/3,
    bt=1/5/2/2,
    tgs=1/5/2/2,
    ppt=1/5/2/5,
    p1=1/5/2/5,
    p2=1/5/2/5,
    p3=1/5/2/5,
    p4=1/5/2/5,
    T150_AWC=1/5/2/11,
    T50_sand=1/5/2/11,
    T150_sand=1/5/2/11,
    T50_clay=1/5/2/11,
    T150_clay=1/5/2/11,
    T50_OM=1/5/2/11,
    T150_OM=1/5/2/11,
    T50_pH=1/5/2/11,
    Water_Table=1/5/2/2,
    wet=1/5/2/2,
    spodic=1/5/2/11,
    carbdepth=1/5/2/11,
    Bhs=1/5/2/11,
    jackpine=1/5/6,
    redpine=1/5/6,
    whitepine=1/5/6,
    aspen=1/5/6,
    oak=1/5/6,
    maple=1/5/6
  )
  
  
  train$resid <- train$chm.max - predict(linear.mod, train)
  
  
  rf <- ranger(resid~
                 solar+
                 #hilly+
                 #age+
                 pred1+
                 pred2+
                 toip+
                 toin+
                 bt+
                 tgs+
                 ppt+
                 p1+
                 p2+
                 p3+
                 p4+
                 T150_AWC+
                 T50_sand+
                 T150_sand+
                 T50_clay+
                 T150_clay+
                 T50_OM+
                 T150_OM+
                 T50_pH+
                 Water_Table+
                 wet+
                 spodic+
                 carbdepth+
                 Bhs+
                 jackpine+
                 redpine+
                 whitepine+
                 aspen+
                 oak+
                 maple
               , data=train,
               num.trees=500, sample.fraction =.2,
               
               split.select.weights = splitwts,
               # max.depth = 5,  
               # mtry=7,
               write.forest = TRUE, importance = 'impurity', case.weights= train$wts
  )
     return(rf)
}

brk.pts.50 <- brk.pts %>% mutate(age=50, pred1 = asmod(50), pred2 = rpmod(50), aspen=0, oak=0, maple=1, jackpine=0, redpine=0, whitepine=0
                                 # ,
                                 # mi=ppt/(bt*58.93+ppt),
                                 # solar=median(solar),
                                 # toip=median(toip),
                                 # toip=median(toin),
                                 # bt=median(bt),
                                 # tgs=median(tgs),
                                 # ppt=median(ppt),
                                 # p1=median(p1),
                                 # p2=median(p2),
                                 # p3=median(p3),
                                 # p4=median(p4),
                                 # T150_AWC=median(T150_AWC),
                                 # T50_sand=median(T50_sand),
                                 # T150_sand=median(T150_sand),
                                 # T50_clay=median(T50_clay),
                                 # T150_clay=median(T150_clay),
                                 # T50_OM=median(T50_OM),
                                 # T150_OM=median(T150_OM),
                                 # T50_pH=median(T50_pH),
                                 # Water_Table=median(Water_Table),
                                 # wet=median(wet),
                                 # hydric=median(hydric),
                                 # moist=median(moist),
                                 # spodosols=median(spodosols),
                                 # Bhs=median(Bhs)
                                 
                              
)


brk.pts.50 <- brk.pts.50 %>% mutate(chm = predict(linear.mod, newdata=brk.pts.50))

generic.si <- rast(cbind(x=brk.pts.50$x, y=brk.pts.50$y, generic.si=brk.pts.50$chm), type="xyz", crs=crs(dem))
plot(generic.si)
writeRaster(generic.si, 'output/generic.si.lm.tif', overwrite=T)










#jack pine si ----
testingrows <- sample(rownames(hmnf.chm2), floor(nrow(hmnf.chm2)*0.1))

train <- hmnf.chm2[!rownames(hmnf.chm2) %in% testingrows,]
test <- hmnf.chm2[rownames(hmnf.chm2) %in% testingrows,]

train$wts <-   ((train$jackpine %in% 1)*1/mean(train$jackpine)/1+1)*((train$hilly >= 0.5)*3+1)
linear.mod <- lmod(train)
test <- test %>% mutate(chm.mod = predict(linear.mod, newdata=test))
Metrics::rmse(test$chm.max,test$chm.mod)

hy.mod <- hybrid.model(train)
test <- test %>% mutate(resid = predictions(predict(hy.mod, data=test)), chm.lm = predict(linear.mod, newdata=test), chm.mod = chm.lm+resid)
Metrics::rmse(test$chm.max,test$chm.mod)

brk.pts.50 <- brk.pts %>% mutate(age=50, pred1 = asmod(50), pred2 = rpmod(50), aspen=0, oak=0, maple=0, jackpine=1, redpine=0, whitepine=0)
brk.pts.50 <- brk.pts.50 %>% mutate(resid = predictions(predict(hy.mod, data=brk.pts.50)), chm.lm = predict(linear.mod, newdata=brk.pts.50), chm=chm.lm+resid)

jackpine.si <- rast(cbind(x=brk.pts.50$x, y=brk.pts.50$y, jackpine.si=brk.pts.50$chm), type="xyz", crs=crs(dem))
writeRaster(jackpine.si, 'output/jackpine.si.tif', overwrite=T)

#red pine si  ----
testingrows <- sample(rownames(hmnf.chm2), floor(nrow(hmnf.chm2)*0.1))

train <- hmnf.chm2[!rownames(hmnf.chm2) %in% testingrows,]
test <- hmnf.chm2[rownames(hmnf.chm2) %in% testingrows,]

train$wts <-   ((train$redpine %in% 1)*1/mean(train$redpine)/1+1)*((train$hilly >= 0.5)*3+1)
linear.mod <- lmod(train)
test <- test %>% mutate(chm.mod = predict(linear.mod, newdata=test))
Metrics::rmse(test$chm.max,test$chm.mod)

hy.mod <- hybrid.model(train)
test <- test %>% mutate(resid = predictions(predict(hy.mod, data=test)), chm.lm = predict(linear.mod, newdata=test), chm.mod = chm.lm+resid)
Metrics::rmse(test$chm.max,test$chm.mod)

brk.pts.50 <- brk.pts %>% mutate(age=50, pred1 = asmod(50), pred2 = rpmod(50), aspen=0, oak=0, maple=0, jackpine=0, redpine=1, whitepine=0)
brk.pts.50 <- brk.pts.50 %>% mutate(resid = predictions(predict(hy.mod, data=brk.pts.50)), chm.lm = predict(linear.mod, newdata=brk.pts.50), chm=chm.lm+resid)
redpine.si <- rast(cbind(x=brk.pts.50$x, y=brk.pts.50$y, redpine.si=brk.pts.50$chm), type="xyz", crs=crs(dem))

writeRaster(redpine.si, 'output/redpine.si.tif', overwrite=T)
#whitepine si  ----
testingrows <- sample(rownames(hmnf.chm2), floor(nrow(hmnf.chm2)*0.1))

train <- hmnf.chm2[!rownames(hmnf.chm2) %in% testingrows,]
test <- hmnf.chm2[rownames(hmnf.chm2) %in% testingrows,]

train$wts <-   ((train$whitepine %in% 1)*1/mean(train$whitepine)/1+1)*((train$hilly >= 0.5)*3+1)
linear.mod <- lmod(train)
test <- test %>% mutate(chm.mod = predict(linear.mod, newdata=test))
Metrics::rmse(test$chm.max,test$chm.mod)

hy.mod <- hybrid.model(train)
test <- test %>% mutate(resid = predictions(predict(hy.mod, data=test)), chm.lm = predict(linear.mod, newdata=test), chm.mod = chm.lm+resid)
Metrics::rmse(test$chm.max,test$chm.mod)

brk.pts.50 <- brk.pts %>% mutate(age=50, pred1 = asmod(50), pred2 = rpmod(50), aspen=0, oak=0, maple=0, jackpine=0, redpine=0, whitepine=1)
brk.pts.50 <- brk.pts.50 %>% mutate(resid = predictions(predict(hy.mod, data=brk.pts.50)), chm.lm = predict(linear.mod, newdata=brk.pts.50), chm=chm.lm+resid)

whitepine.si <- rast(cbind(x=brk.pts.50$x, y=brk.pts.50$y, whitepine.si=brk.pts.50$chm), type="xyz", crs=crs(dem))

writeRaster(whitepine.si, 'output/whitepine.si.tif', overwrite=T)

#aspen si  ----
testingrows <- sample(rownames(hmnf.chm2), floor(nrow(hmnf.chm2)*0.1))

train <- hmnf.chm2[!rownames(hmnf.chm2) %in% testingrows,]
test <- hmnf.chm2[rownames(hmnf.chm2) %in% testingrows,]

train$wts <-   ((train$aspen %in% 1)*1/mean(train$aspen)/1+1)*((train$hilly >= 0.5)*3+1)
linear.mod <- lmod(train)
test <- test %>% mutate(chm.mod = predict(linear.mod, newdata=test))
Metrics::rmse(test$chm.max,test$chm.mod)

hy.mod <- hybrid.model(train)
test <- test %>% mutate(resid = predictions(predict(hy.mod, data=test)), chm.lm = predict(linear.mod, newdata=test), chm.mod = chm.lm+resid)
Metrics::rmse(test$chm.max,test$chm.mod)

brk.pts.50 <- brk.pts %>% mutate(age=50, pred1 = asmod(50), pred2 = rpmod(50), aspen=1, oak=0, maple=0, jackpine=0, redpine=0, whitepine=0)
brk.pts.50 <- brk.pts.50 %>% mutate(resid = predictions(predict(hy.mod, data=brk.pts.50)), chm.lm = predict(linear.mod, newdata=brk.pts.50), chm=chm.lm+resid)

aspen.si <- rast(cbind(x=brk.pts.50$x, y=brk.pts.50$y, aspen.si=brk.pts.50$chm), type="xyz", crs=crs(dem))

writeRaster(aspen.si, 'output/aspen.si.tif', overwrite=T)

#oak si  ----
testingrows <- sample(rownames(hmnf.chm2), floor(nrow(hmnf.chm2)*0.1))

train <- hmnf.chm2[!rownames(hmnf.chm2) %in% testingrows,]
test <- hmnf.chm2[rownames(hmnf.chm2) %in% testingrows,]

train$wts <-   ((train$oak %in% 1)*1/mean(train$oak)/1+1)*((train$hilly >= 0.5)*3+1)
linear.mod <- lmod(train)
test <- test %>% mutate(chm.mod = predict(linear.mod, newdata=test))
Metrics::rmse(test$chm.max,test$chm.mod)

hy.mod <- hybrid.model(train)
test <- test %>% mutate(resid = predictions(predict(hy.mod, data=test)), chm.lm = predict(linear.mod, newdata=test), chm.mod = chm.lm+resid)
Metrics::rmse(test$chm.max,test$chm.mod)

brk.pts.50 <- brk.pts %>% mutate(age=50, pred1 = asmod(50), pred2 = rpmod(50), aspen=0, oak=1, maple=0, jackpine=0, redpine=0, whitepine=0)
brk.pts.50 <- brk.pts.50 %>% mutate(resid = predictions(predict(hy.mod, data=brk.pts.50)), chm.lm = predict(linear.mod, newdata=brk.pts.50), chm=chm.lm+resid)

oak.si <- rast(cbind(x=brk.pts.50$x, y=brk.pts.50$y, oak.si=brk.pts.50$chm), type="xyz", crs=crs(dem))

writeRaster(oak.si, 'output/oak.si.tif', overwrite=T)


#maple si  ----
testingrows <- sample(rownames(hmnf.chm2), floor(nrow(hmnf.chm2)*0.1))

train <- hmnf.chm2[!rownames(hmnf.chm2) %in% testingrows,]
test <- hmnf.chm2[rownames(hmnf.chm2) %in% testingrows,]

train$wts <-   ((train$maple %in% 1)*1/mean(train$maple)/1+1)*((train$hilly >= 0.5)*3+1)
linear.mod <- lmod(train)
test <- test %>% mutate(chm.mod = predict(linear.mod, newdata=test))
Metrics::rmse(test$chm.max,test$chm.mod)

hy.mod <- hybrid.model(train)
test <- test %>% mutate(resid = predictions(predict(hy.mod, data=test)), chm.lm = predict(linear.mod, newdata=test), chm.mod = chm.lm+resid)
Metrics::rmse(test$chm.max,test$chm.mod)

rf.mod <- rf.model(train)
test <- test %>% mutate(chm.mod = predictions(predict(rf.mod, data=test)))
Metrics::rmse(test$chm.max,test$chm.mod)

brk.pts.50 <- brk.pts %>% mutate(age=50, pred1 = asmod(50), pred2 = rpmod(50), aspen=0, oak=0, maple=1, jackpine=0, redpine=0, whitepine=0)
brk.pts.50 <- brk.pts.50 %>% mutate(chm.rf = predictions(predict(rf.mod, data=brk.pts.50)), resid = predictions(predict(hy.mod, data=brk.pts.50)), chm.lm = predict(linear.mod, newdata=brk.pts.50), chm=chm.lm+resid)

maple.si <- rast(cbind(x=brk.pts.50$x, y=brk.pts.50$y, chm=brk.pts.50$chm), type="xyz", crs=crs(dem))
plot(maple.si)
writeRaster(maple.si, 'output/maple.si.tif', overwrite=T)

maple.si.rf <- rast(cbind(x=brk.pts.50$x, y=brk.pts.50$y, chm=brk.pts.50$chm.rf), type="xyz", crs=crs(dem))
plot(maple.si.rf)
writeRaster(maple.si.rf, 'output/maple.si.rf.tif', overwrite=T)

maple.si.lm <- rast(cbind(x=brk.pts.50$x, y=brk.pts.50$y, chm=brk.pts.50$chm.lm), type="xyz", crs=crs(dem))
plot(maple.si.lm)
writeRaster(maple.si.lm, 'output/maple.si.lm.tif', overwrite=T)
#----

focalsoil <- subset(hmnf.chm2, Lat > 44.2 & Lat < 44.4 & Long > -85.75 & Long < -85.6)
brk.pts.50 <- brk.pts %>% mutate(age=50, pred1 = asmod(50), pred2 = rpmod(50), aspen=0, oak=0, maple=1, jackpine=0, redpine=0, whitepine=0,
                                 # solar=median(focalsoil$solar),
                                 # hilly=median(focalsoil$hilly),
                                 # toip=median(focalsoil$toip),
                                 # toin=median(focalsoil$toin),
                                 bt=median(focalsoil$bt),
                                 tgs=median(focalsoil$tgs),
                                 ppt=median(focalsoil$ppt),
                                 p1=median(focalsoil$p1),
                                 p2=median(focalsoil$p2),
                                 p3=median(focalsoil$p3),
                                 p4=median(focalsoil$p4),
                                 T150_AWC=median(focalsoil$T150_AWC),
                                 T50_sand=median(focalsoil$T50_sand),
                                 T150_sand=median(focalsoil$T150_sand),
                                 T50_clay=median(focalsoil$T50_clay),
                                 T150_clay=median(focalsoil$T150_clay),
                                 T50_OM=median(focalsoil$T50_OM),
                                 T150_OM=median(focalsoil$T150_OM),
                                 T50_pH=median(focalsoil$T50_pH),
                                 Water_Table=median(focalsoil$Water_Table),
                                 wet=median(focalsoil$wet),
                                 spodic=median(focalsoil$spodic),
                                 spodosols=median(focalsoil$spodosols),
                                 carbdepth=median(focalsoil$carbdepth),
                                 hydric=median(focalsoil$hydric),
                                 moist=median(focalsoil$moist),
                                 Bhs=median(focalsoil$Bhs))

brk.pts.50 <- brk.pts.50 %>% mutate(resid = predictions(predict(hy.mod, data=brk.pts.50)), chm.lm = predict(linear.mod, newdata=brk.pts.50), chm=chm.lm+resid)

maple.si.topo <- rast(cbind(x=brk.pts.50$x, y=brk.pts.50$y, maple.si=brk.pts.50$chm), type="xyz", crs=crs(dem))
plot(maple.si.topo)
writeRaster(maple.si.topo, 'output/maple.si.topo.tif', overwrite=T)
saveRDS(brk.pts.50, 'tmp/brk.pts.50.topo.RDS')
#ggplot of model ----
focalsoil <- subset(hmnf.chm2, Lat > 44.2 & Lat < 44.4 & Long > -85.75 & Long < -85.6)
brk.pts.50 <- hmnf.chm2 %>% mutate(age=50, pred1 = asmod(50), pred2 = rpmod(50), aspen=0, oak=0, maple=1, jackpine=0, redpine=0, whitepine=0,
                                 # solar=median(focalsoil$solar),
                                 # hilly=median(focalsoil$hilly),
                                 # toip=median(focalsoil$toip),
                                 # toin=median(focalsoil$toin),
                                 bt=median(focalsoil$bt),
                                 tgs=median(focalsoil$tgs),
                                 ppt=median(focalsoil$ppt),
                                 p1=median(focalsoil$p1),
                                 p2=median(focalsoil$p2),
                                 p3=median(focalsoil$p3),
                                 p4=median(focalsoil$p4),
                                 T150_AWC=median(focalsoil$T150_AWC),
                                 T50_sand=median(focalsoil$T50_sand),
                                 T150_sand=median(focalsoil$T150_sand),
                                 T50_clay=median(focalsoil$T50_clay),
                                 T150_clay=median(focalsoil$T150_clay),
                                 T50_OM=median(focalsoil$T50_OM),
                                 T150_OM=median(focalsoil$T150_OM),
                                 T50_pH=median(focalsoil$T50_pH),
                                 Water_Table=median(focalsoil$Water_Table),
                                 wet=median(focalsoil$wet),
                                 spodic=median(focalsoil$spodic),
                                 spodosols=median(focalsoil$spodosols),
                                 carbdepth=median(focalsoil$carbdepth),
                                 hydric=median(focalsoil$hydric),
                                 moist=median(focalsoil$moist),
                                 Bhs=median(focalsoil$Bhs))

brk.pts.50 <- brk.pts.50 %>% mutate(resid = predictions(predict(hy.mod, data=brk.pts.50)), chm.lm = predict(linear.mod, newdata=brk.pts.50), chm=chm.lm+resid)





brk.pts.50$toi <- brk.pts.50$toip -brk.pts.50$toin
brk.pts.50$toi <- (brk.pts.50$toi - min(brk.pts.50$toi))/(max(brk.pts.50$toi)-min(brk.pts.50$toi))
brk.pts.50$normsolar <- (brk.pts.50$solar - mean(brk.pts.50$solar))/(max(brk.pts.50$solar)-min(brk.pts.50$solar))
mean(brk.pts.50$solar)
# 
# 
# breaks.solar <- (max(brk.pts.50$solar)-min(brk.pts.50$solar))/20*c(0:20)+min(brk.pts.50$solar)
# breaks.toin <- (max(brk.pts.50$toin)-min(brk.pts.50$toin))/20*c(0:20)+min(brk.pts.50$toin)
# breaks.toip <- (max(brk.pts.50$toip)-min(brk.pts.50$toip))/20*c(0:20)+min(brk.pts.50$toip)
# i=1
# for(i in 1:20){
#   for(j in 1:20){  
#     for(k in 1:20){
# chm <- mean(subset(brk.pts.50, toin >= breaks.toin[i] & toin >= breaks.toin[i+1] &
#                      toip >= breaks.toip[j] & toip >= breaks.toip[j+1] &
#                      solar >= breaks.solar[j] & solar >= breaks.solar[j+1] )$chm)
# toin <- mean(breaks.toin[i],breaks.toin[i+1])
# toip <- mean(breaks.toip[j],breaks.toip[j+1])
# solar <- mean(breaks.solar[k],breaks.solar[k+1])
# 
# if(i==1&j==1){simple.df<- as.data.frame(cbind(chm,toin,toip,solar))}else{simple.df<- rbind(simple.df,as.data.frame(cbind(chm,toin,toip,solar)))}
# }}}
# saveRDS(simple.df, 'tmp/simple.df.RDS')
# 
# 
# Get percentiles for each slope position and sun position and then plot all three factors
# 
ggplot()+
  geom_smooth(aes(x=toi, y=chm, col='Maple'), data= brk.pts.50)+
  
  scale_color_manual(name='site',  values = c('model'='black'
  ))+
  scale_y_continuous(name='Canopy Height (m)')+
  scale_x_continuous(name='position')+
  #coord_fixed(ratio=2/1,ylim=c(0,40))+
  labs(title = 'Maple SI by slope position')

ggplot()+
  geom_smooth(aes(x=normsolar, y=chm, col='Maple'), data= brk.pts.50)+
  
  scale_color_manual(name='site',  values = c('model'='black'
  ))+
  scale_y_continuous(name='Canopy Height (m)')+
  scale_x_continuous(name='solar')+
  #coord_fixed(ratio=2/1,ylim=c(0,40))+
  labs(title = 'Maple SI by solar exposure')

ggplot()+
  geom_smooth(aes(x=ppt, y=chm, col='Maple'), data= brk.pts.50)+
  
  scale_color_manual(name='site',  values = c('model'='black'
  ))+
  scale_y_continuous(name='Canopy Height (m)')+
  scale_x_continuous(name='MAP')+
  #coord_fixed(ratio=2/1,ylim=c(0,40))+
  labs(title = 'Maple SI by Mean Annual Precipitation')




#boosted ----
# require(gbm)
# require(MASS)#package with the boston housing dataset

# rf.maple <- gbm(chm.max~
#                      solar+
#                      #hilly+
#                      #age+
#                      pred1+
#                      pred2+
#                      toip+
#                      toin+
#                      bt+
#                      tgs+
#                      ppt+
#                      p1+
#                      p2+
#                      p3+
#                      p4+
#                      T150_AWC+
#                      T50_sand+
#                      T150_sand+
#                      T50_clay+
#                      T150_clay+
#                      T50_OM+
#                      T150_OM+
#                      T50_pH+
#                      Water_Table+
#                      wet+
#                      spodic+
#                      carbdepth+
#                      Bhs+
#                      jackpine+
#                      redpine+
#                      whitepine+
#                      aspen+
#                      oak+
#                      maple
#                    , data=train,
#                 distribution = "gaussian",
#                 n.trees = 500,
#                 train.fraction = 1
#                 
# )
# 
# 
# test <- test %>% mutate(chm.mod = predict.gbm(rf.maple, newdata=test))
# Metrics::rmse(test$chm.max,test$chm.mod)

#linear tree ----
# library(partykit)
# lm.maple <- lmtree(chm.max~
#                      pred1+
#                      pred1:ppt+
#                      pred1:tgs+
#                      pred1:toip+
#                      pred1:solar+
#                      pred1:Water_Table+
#                      pred1:T150_OM+
#                      pred1:T150_AWC|
#                      
#                      
#                      pred2+
#                      toin+
#                      bt+
#                      p1+
#                      p2+
#                      p3+
#                      p4+                     
#                      T50_sand+
#                      T150_sand+
#                      T50_clay+
#                      T150_clay+
#                      T50_OM+
#                      T50_pH+
#                      wet+
#                      spodic+
#                      carbdepth+
#                      Bhs+
#                      jackpine+
#                      redpine+
#                      whitepine+
#                      aspen+
#                      oak+
#                      maple
#                    , data=hmnf.chm2
# )
# 

#jack pine prob ----
hmnf.chm2$wts <-  1
rf.jackpine <- ranger(jackpine ~
                        solar+
                        hilly+
                        toip+
                        toin+
                        bt+
                        tgs+
                        ppt+
                        T150_AWC+
                        T50_sand+
                        T150_sand+
                        T50_clay+
                        T150_clay+
                        T50_OM+
                        T150_OM+
                        T50_pH+
                        Water_Table+
                        wet+
                        spodic+
                        carbdepth+
                        Bhs
                      , data=hmnf.chm2,
                      num.trees=500, sample.fraction =0.02,
                      write.forest = TRUE, importance = 'impurity'
)

brk.pts.50 <- brk.pts
brk.pts.50 <- brk.pts.50 %>% mutate(prob = predictions(predict(rf.jackpine, data=brk.pts.50)))

jackpine.prob <- rast(cbind(x=brk.pts.50$x, y=brk.pts.50$y, jackpine.prob=brk.pts.50$prob), type="xyz", crs=crs(dem))

writeRaster(jackpine.prob, 'output/jackpine.prob.tif', overwrite=T)
#red pine prob ----
hmnf.chm2$wts <-  1
rf.redpine <- ranger(redpine ~
                        solar+
                        hilly+
                        toip+
                        toin+
                        bt+
                        tgs+
                        ppt+
                        T150_AWC+
                        T50_sand+
                        T150_sand+
                        T50_clay+
                        T150_clay+
                        T50_OM+
                        T150_OM+
                        T50_pH+
                        Water_Table+
                        wet+
                        spodic+
                        carbdepth+
                        Bhs
                      , data=hmnf.chm2,
                      num.trees=500, sample.fraction =0.02,
                      write.forest = TRUE, importance = 'impurity'
)

brk.pts.50 <- brk.pts
brk.pts.50 <- brk.pts.50 %>% mutate(prob = predictions(predict(rf.redpine, data=brk.pts.50)))

redpine.prob <- rast(cbind(x=brk.pts.50$x, y=brk.pts.50$y, redpine.prob=brk.pts.50$prob), type="xyz", crs=crs(dem))

writeRaster(redpine.prob, 'output/redpine.prob.tif', overwrite=T)
#whitepine prob ----
hmnf.chm2$wts <-  1
rf.whitepine <- ranger(whitepine ~
                        solar+
                        hilly+
                        toip+
                        toin+
                        bt+
                        tgs+
                        ppt+
                        T150_AWC+
                        T50_sand+
                        T150_sand+
                        T50_clay+
                        T150_clay+
                        T50_OM+
                        T150_OM+
                        T50_pH+
                        Water_Table+
                        wet+
                        spodic+
                        carbdepth+
                        Bhs
                      , data=hmnf.chm2,
                      num.trees=500, sample.fraction =0.02,
                      write.forest = TRUE, importance = 'impurity'
)

brk.pts.50 <- brk.pts
brk.pts.50 <- brk.pts.50 %>% mutate(prob = predictions(predict(rf.whitepine, data=brk.pts.50)))

whitepine.prob <- rast(cbind(x=brk.pts.50$x, y=brk.pts.50$y, whitepine.prob=brk.pts.50$prob), type="xyz", crs=crs(dem))

writeRaster(whitepine.prob, 'output/whitepine.prob.tif', overwrite=T)
#aspen prob ----
hmnf.chm2$wts <-  1
rf.aspen <- ranger(aspen ~
                        solar+
                        hilly+
                        toip+
                        toin+
                        bt+
                        tgs+
                        ppt+
                        T150_AWC+
                        T50_sand+
                        T150_sand+
                        T50_clay+
                        T150_clay+
                        T50_OM+
                        T150_OM+
                        T50_pH+
                        Water_Table+
                        wet+
                        spodic+
                        carbdepth+
                        Bhs
                      , data=hmnf.chm2,
                      num.trees=500, sample.fraction =0.02,
                      write.forest = TRUE, importance = 'impurity'
)

brk.pts.50 <- brk.pts
brk.pts.50 <- brk.pts.50 %>% mutate(prob = predictions(predict(rf.aspen, data=brk.pts.50)))

aspen.prob <- rast(cbind(x=brk.pts.50$x, y=brk.pts.50$y, aspen.prob=brk.pts.50$prob), type="xyz", crs=crs(dem))

writeRaster(aspen.prob, 'output/aspen.prob.tif', overwrite=T)
#oak prob ----
hmnf.chm2$wts <-  1
rf.oak <- ranger(oak ~
                        solar+
                        hilly+
                        toip+
                        toin+
                        bt+
                        tgs+
                        ppt+
                        T150_AWC+
                        T50_sand+
                        T150_sand+
                        T50_clay+
                        T150_clay+
                        T50_OM+
                        T150_OM+
                        T50_pH+
                        Water_Table+
                        wet+
                        spodic+
                        carbdepth+
                        Bhs
                      , data=hmnf.chm2,
                      num.trees=500, sample.fraction =0.02,
                      write.forest = TRUE, importance = 'impurity'
)

brk.pts.50 <- brk.pts
brk.pts.50 <- brk.pts.50 %>% mutate(prob = predictions(predict(rf.oak, data=brk.pts.50)))

oak.prob <- rast(cbind(x=brk.pts.50$x, y=brk.pts.50$y, oak.prob=brk.pts.50$prob), type="xyz", crs=crs(dem))

writeRaster(oak.prob, 'output/oak.prob.tif', overwrite=T)
#maple prob ----
hmnf.chm2$wts <-  1
rf.maple <- ranger(maple ~
                        solar+
                        hilly+
                        toip+
                        toin+
                        bt+
                        tgs+
                        ppt+
                        T150_AWC+
                        T50_sand+
                        T150_sand+
                        T50_clay+
                        T150_clay+
                        T50_OM+
                        T150_OM+
                        T50_pH+
                        Water_Table+
                        wet+
                        spodic+
                        carbdepth+
                        Bhs
                      , data=hmnf.chm2,
                      num.trees=500, sample.fraction =0.02,
                      write.forest = TRUE, importance = 'impurity'
)

brk.pts.50 <- brk.pts
brk.pts.50 <- brk.pts.50 %>% mutate(prob = predictions(predict(rf.maple, data=brk.pts.50)))

maple.prob <- rast(cbind(x=brk.pts.50$x, y=brk.pts.50$y, maple.prob=brk.pts.50$prob), type="xyz", crs=crs(dem))

writeRaster(maple.prob, 'output/maple.prob.tif', overwrite=T)

plot(maple.prob)

v1 <- (jackpine.prob >1/6)*(jackpine.prob >= max(jackpine.prob, whitepine.prob, oak.prob, maple.prob, aspen.prob, redpine.prob))
v3 <- (whitepine.prob >1/6)*(whitepine.prob >= max(jackpine.prob,  whitepine.prob, oak.prob, maple.prob, aspen.prob, redpine.prob))*3
v5 <- (oak.prob >1/6)*(oak.prob >= max(jackpine.prob,  whitepine.prob, oak.prob, maple.prob, aspen.prob, redpine.prob))*5
v6 <- (maple.prob >1/6)*(maple.prob >= max(jackpine.prob,  whitepine.prob, oak.prob, maple.prob, aspen.prob, redpine.prob))*6
v4 <- (aspen.prob >1/6)*(aspen.prob >= max(jackpine.prob,  whitepine.prob, oak.prob, maple.prob, aspen.prob, redpine.prob))*4#*(max(v1,v3,v5,v6)==0)
v2 <- (redpine.prob >1/6)*(redpine.prob >= max(jackpine.prob,  whitepine.prob, oak.prob, maple.prob, aspen.prob, redpine.prob))#*2*(max(v1,v3,v4,v5,v6)==0)
vegbest <- max(v1,v2, v3, v4, v5, v6)
plot(vegbest)

writeRaster(vegbest, 'output/vegbest.tif', overwrite=T)
#fixed soil slope relations ---
hmnf.chm2$wts <-  (hmnf.chm2$maple %in% 1)*3+(hmnf.chm2$hilly >= 0.5)*3+1
#To maximize topographic difference in model, I tried more or less weighting of the hilly parameter, and using it as a model input, but it is best to have moderate weight of 3 than no weight or weight of 6, and better not to have as a model input. Also not as good to include toi as an input.
rf.maple <- ranger(chm.max~
                     solar+
                     #hilly+
                     #age+
                     pred1+
                     pred2+
                     toip+
                     toin+
                     bt+
                     tgs+
                     ppt+
                     T150_AWC+
                     T50_sand+
                     T150_sand+
                     T50_clay+
                     T150_clay+
                     T50_OM+
                     T150_OM+
                     T50_pH+
                     Water_Table+
                     wet+
                     spodic+
                     carbdepth+
                     Bhs+
                     jackpine+
                     redpine+
                     whitepine+
                     aspen+
                     oak+
                     maple
                   , data=hmnf.chm2,
                   num.trees=500, sample.fraction =0.02, always.split.variables = c('pred1', 'maple'),
                   write.forest = TRUE, importance = 'impurity', case.weights= hmnf.chm2$wts
)

brk.pts.sand <- subset(brk.pts, Bhs >0.5 & Water_Table > 150 & T150_sand > 0.75)
brk.pts.sand <- brk.pts.sand %>% mutate(T150_AWC = mean(T150_AWC),
                                          T50_sand = mean(T50_sand),
                                          T150_sand = mean(T150_sand),
                                          T50_clay = mean(T50_clay),
                                          T150_clay = mean(T150_clay),
                                          T50_OM = mean(T50_OM),
                                          T150_OM = mean(T150_OM),
                                          T50_pH = mean(T50_pH)
                                          )
  
  
brk.pts.50 <- brk.pts %>% mutate(age=50, pred1 = asmod(50), pred2 = rpmod(50), aspen=0, oak=0, maple=1, jackpine=0, redpine=0, whitepine=0,
                                 
                                 
                                 
                                 
                                 
                                 
                                 )
brk.pts.50 <- brk.pts.50 %>% mutate(chm = predictions(predict(rf.maple, data=brk.pts.50)))





### other mountains  ----
### 

library(sf)
library(terra)
library(raster)
library(minpack.lm)
library(growthmodels)
library(dplyr)
library(ggplot2)
library(Hmisc)
library(MASS)
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

katahdin <- rast('data/lidar/katahdin.tif')
bullhead <- rast('data/lidar/bullrun.tif')
porters <- rast('data/lidar/porters.tif')
port.solar <- rast('data/lidar/porters_solar.tif')
port.solar.annual <- mean(port.solar)
plot(port.solar.annual)
names(porters) <- c("dem","slope","aspect","toip","toin","twi","chm")
names(port.solar) <- c('d0101','d0127','d0222','d0319','d0414','d0510','d0605','d0701','d0727','d0822','d0917','d1013','d1108','d1204','d1230')

sveg <- st_read('D:/scripts/lidar/data/vegpolygons/GRSM_veg.shp')
sveg <- st_transform(sveg, crs(porters))
sveg.r <- vect(sveg)
sveg.r <- rasterize(sveg.r, porters$dem, field = 'GRSM_code')
sveg.tab <- sveg[,c('GRSM_code', 'Short_name','Ecogroup')] %>% st_drop_geometry() %>% unique()
  
  
  
names(katahdin) <- c("dem", "toip","toin","tpi","slope","aspect","chm")
kat <- as.data.frame(katahdin, xy=T)

names(bullhead) <- c("dem", "toip","toin","tpi","slope","aspect","chm")
bull <- as.data.frame(bullhead, xy=T)


port <- as.data.frame(c(porters, port.solar,sveg.r), xy=T)
port <- merge(port, sveg.tab, by='GRSM_code', all.x=T)

bull$toi <- (bull$toip-bull$toin)/2
kat$toi <- (kat$toip-kat$toin)/2

bull$south <- bull$slope^2+cos(bull$aspect+3.141592)
bull$west <- bull$slope^2+cos(bull$aspect+3.141592*1.5)
bull$southwest <- bull$slope^2+cos(bull$aspect+3.141592*1.25)
bull$southeast <- bull$slope^2+cos(bull$aspect+3.141592*0.75)
kat$south <- kat$slope^2+cos(kat$aspect+3.141592)
kat$west <- kat$slope^2+cos(kat$aspect+3.141592*1.5)
kat$southwest <- kat$slope^2+cos(kat$aspect+3.141592*1.25)
kat$southeast <- kat$slope^2+cos(kat$aspect+3.141592*0.75)


bull <- subset(bull, !is.na(toip)& !is.na(chm))
kat <- subset(kat, !is.na(toip)& !is.na(chm))
ggplot()+
  geom_smooth(aes(x=toi, y=chm, col='GSMNP'), data= bull)+
  geom_smooth(aes(x=toi, y=chm, col='Katahdin'), data= kat)+
  scale_color_manual(name='site',  values = c('GSMNP'='green','Katahdin'='red'))

c.table <- as.data.frame(cor(bull[,c("chm", "toi","toip","toin","tpi","dem","slope","south","west","southwest","southeast")], use='complete.obs'
))

port$saspect <- ifelse(port$slope > 0.25, (port$aspect)/2/3.141592*360, NA)
port$toi <- port$toip-port$toin
port$solar <- apply(port[,c('d0101','d0127','d0222','d0319','d0414','d0510','d0605','d0701','d0727','d0822','d0917','d1013','d1108','d1204','d1230')], MARGIN=1, FUN='mean')

port$toi2 <- (65.81*port$toin + -103.32*port$toip)/-(65.81+103.32)

port2<- subset(port, !GRSM_code %in% c(1015, 1014, 1003, 1019, 4048, 1006))


saveRDS(port, 'output/port.RDS')
saveRDS(kat, 'output/kat.RDS')
saveRDS(bull, 'output/bull.RDS')
port <- readRDS('output/port.RDS')
kat <- readRDS('output/kat.RDS')
bull <- readRDS('output/bull.RDS')


ggplot()+
  geom_smooth(aes(x=saspect, y=chm, col='GSMNP'), data= port)+
  scale_color_manual(name='site',  values = c('GSMNP'='green'))+
  scale_y_continuous(name='Canopy Height (m)', breaks=c((0:10)*5))+
  scale_x_continuous(name='Aspect', breaks=c((0:8)*45), labels=c("N","NE","E","SE","S","SW","W","NW","N"))
ggplot()+
  geom_smooth(aes(x=twi, y=chm, col='GSMNP'), data= port)+
  scale_color_manual(name='site',  values = c('GSMNP'='green'))+
  scale_y_continuous(name='Canopy Height (m)', breaks=c((0:10)*5))+
  scale_x_continuous(name='twi')



c.table <- as.data.frame(cor(port[,c('chm','d0101','d0127','d0222','d0319','d0414','d0510','d0605','d0701','d0727','d0822','d0917','d1013','d1108','d1204','d1230')], use='complete.obs'
))

model <- lm(chm~
              toip+tpi+dem+south
            , data=bull)
stepAIC(model)
summary(model)

model <- lm(chm~
              tpi+dem+south
            , data=kat)
stepAIC(model)
summary(model)


model <- lm(chm~
              toip+tpi+dem+south
            , data=bull)
stepAIC(model)
summary(model)

model <- lm(chm~
              toip+toin+dem+d1230
            , data=port2)
stepAIC(model)
summary(model)

ggplot()+
  geom_smooth(aes(x=toi2, y=chm, col='GSMNP'), data= port)+
  scale_color_manual(name='site',  values = c('GSMNP'='green'))+
  scale_y_continuous(name='Canopy Height (m)', breaks=c((0:10)*5))+
  scale_x_continuous(name='topographic position')

library(rpart)
port$position <- (port$toi2 - median(port$toi2))/(median(port$toi2)-min(port$toi2))*100
port$sunslope <- (port$d1230 - median(port$d1230))/(median(port$d1230)-min(port$d1230))*100
port$elevation <- port$dem
rp <- rpart(chm~
              position+sunslope+elevation, data = port)
rpart.plot::rpart.plot(rp)


ggplot()+
  geom_smooth(aes(x=position, y=chm, col='GSMNP'), data= port)+
  scale_color_manual(name='site',  values = c('GSMNP'='green'))+
  scale_y_continuous(name='Canopy Height (m)', breaks=c((0:10)*5))+
  scale_x_continuous(name='topographic position')
ggplot()+
  geom_smooth(aes(x=sunslope, y=chm, col='GSMNP'), data= port)+
  scale_color_manual(name='site',  values = c('GSMNP'='green'))+
  scale_y_continuous(name='Canopy Height (m)', breaks=c((0:10)*5))+
  scale_x_continuous(name='sunslope')

exposedridge <- subset(port, position >=50 & sunslope >=50)
medianslope <- subset(port,  position >= -50 & position < 50 & sunslope >= -50 & sunslope < 50)
shadycove <- subset(port, position < -50 & sunslope < -50) 
cove <- subset(port, position < -50 & sunslope >= -50 & sunslope < 50) 
ridge <- subset(port, position >= 50 & sunslope >= -50 & sunslope < 50) 
shadyslope <- subset(port, position >= -50 & position < 50 & sunslope < -50) 
sunnyslope <- subset(port, position >= -50 & position < 50 & sunslope >= 50)

ggplot()+
  geom_smooth(aes(x=elevation, y=chm, col='sunny ridge'), data= exposedridge)+
  geom_smooth(aes(x=elevation, y=chm, col='ridge'), data= ridge)+
  geom_smooth(aes(x=elevation, y=chm, col='sunny slope'), data= sunnyslope)+
  geom_smooth(aes(x=elevation, y=chm, col='median slope'), data= medianslope)+
  geom_smooth(aes(x=elevation, y=chm, col='shady slope'), data= shadyslope)+
  geom_smooth(aes(x=elevation, y=chm, col='cove'), data= cove)+
  geom_smooth(aes(x=elevation, y=chm, col='shady cove'), data= shadycove)+
  #geom_smooth(aes(x=elevation, y=chm, col='average slope'), data= port)+
  scale_color_manual(name='site',  values = c('sunny ridge'='red','ridge'='orange','sunny slope'='yellow',
                                              'median slope'='green',
                                              'shady slope'='darkgreen','cove'='cyan','shady cove'='blue'
                                              ))+
  scale_y_continuous(name='Canopy Height (m)', breaks=c((0:10)*5))+
  scale_x_continuous(name='elevation', breaks=c((0:5)*500))

kat$position <- (kat$toi - median(kat$toi))/(median(kat$toi)-min(kat$toi))*100
kat$sunslope <- (kat$south - median(kat$south))/(median(kat$south)-min(kat$south))*100
kat$elevation <- kat$dem
rp <- rpart(chm~
              position+sunslope+elevation, data = kat)
rpart.plot::rpart.plot(rp)

#mtlaff ----
mtlaff <- rast('data/lidar/mtlaff-slopeaspecpositivenegative.tif')
mtlaff.solar <- rast('data/lidar/mtlaff-insolation.tif')
mtlaff.solar.annual <- mean(port.solar)
plot(mtlaff.solar.annual)
names(mtlaff) <- c("dem","slope","aspect","toip","toin","chm")
names(mtlaff.solar) <- c('d0101','d0127','d0222','d0319','d0414','d0510','d0605','d0701','d0727','d0822','d0917','d1013','d1108','d1204','d1230')
lafcomb <- c(mtlaff, mtlaff.solar)
laf <- as.data.frame(lafcomb, xy=T)

model <- lm(chm ~ dem+toip+toin+d1230+d0701+d0414, data=laf)
summary(model)

laf$chm.p <- predict(model, laf)
laf$chm.r <- laf$chm - laf$chm.p
laf <- laf %>% mutate(p.rank = percent_rank(chm.r), chm.rank = percent_rank(chm))
laf.s <- subset(laf, p.rank > 0.05 & p.rank < 0.95 & p.rank < 0.97)
ggplot()+
  geom_point(aes(x=dem, y=chm), data=laf.s, alpha=0.05)

model <- lm(chm ~ dem+toip+toin+d1230+d0701, data=laf.s)
summary(model)
laf$chm.p2 <- predict(model, laf)
pchm.lm <- rast(cbind(x=laf$x,y=laf$y,z=laf$chm.p2), type="xyz", crs=crs(mtlaff$dem))
names(pchm.lm) <- 'pchm.lm'
plot(pchm.lm)
writeRaster(pchm.lm, 'data/lidar/mtlaff-pchm.lm.tif')
library(ranger)
rf <- ranger(chm ~ slope+dem+toip+toin+d1230+d0701, data=laf.s, 
             num.trees=500, sample.fraction =0.01, 
             max.depth = 15,  write.forest = TRUE)
laf$chm.p3 <-predictions(predict(rf, data=laf))
pchm.rf <- rast(cbind(x=laf$x,y=laf$y,z=laf$chm.p3), type="xyz", crs=crs(mtlaff$dem))
names(pchm.rf) <- 'pchm.rf'
plot(pchm.rf)
writeRaster(pchm.rf, 'data/lidar/mtlaff-pchm.rf.tif')
chm.x <- rast(cbind(x=laf.s$x,y=laf.s$y,z=laf.s$chm), type="xyz", crs=crs(mtlaff$dem))
names(chm.x) <- 'chm.x'
plot(chm.x)
writeRaster(chm.x, 'data/lidar/mtlaff-chm.x.tif')

