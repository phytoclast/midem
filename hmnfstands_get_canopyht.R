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


hmnf <- read.csv('data/hmnfstandpoints/SiteIndexRandomPoints_10perStand.txt')
hmnf <- st_as_sf(hmnf, coords=c(x='Long', y='Lat'), crs='epsg:4326', remove = FALSE)
hmnf.mi <- sf::st_transform(hmnf, crs=crs(chm))
hmnf.chm.max <- terra::extract(chm.30.max, vect(hmnf.mi))
hmnf.chm.mean <- terra::extract(chm.30.mean, vect(hmnf.mi))
colnames(hmnf.chm.max) <- c('id','chm.max');colnames(hmnf.chm.mean) <- c('id','chm.mean')
hmnf.chm <- cbind(hmnf.mi, subset(hmnf.chm.max, chm=hmnf.chm[,2]))
hmnf.chm <- cbind(hmnf.chm, subset(hmnf.chm.mean, chm=hmnf.chm[,2]))
hmnf.chm <- hmnf.chm %>% left_join(EV_Code, by=c('EV_CODE'='First_EV_CODE'))
hmnf.chm$age <- 2021-hmnf.chm$YEAR_OF_ORIGIN
hmnf.chm$lmapunitiid <- hmnf.chm$HMNF_MapunitRaster_10m
colnames(hmnf.chm)
hmnf.chm <-  subset(hmnf.chm, select=c("Lat","Long","EV_CODE","YEAR_OF_ORIGIN","SITE_INDEX","SITE_INDEX_SPP","SITE_INDEX_REF",
                                       "lmapunitiid", "chm.max", "chm.mean","type","age")) %>% st_drop_geometry()
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
# slope.10 <- aggregate(slope, fact= 10, fun='mean', na.rm=T)
# slope.focal <- focal(slope.10, focalWeight(slope.10, 500, "Gauss"))
# slope.500 <- resample(slope.focal, slope)
# writeRaster(slope.500, 'D:/GIS/DEM/slope.500.tif', overwrite=T)
slope.500 <- rast('D:/GIS/DEM/slope.500.tif'); names(slope.500) = 'slope.500'

bt <- project(bt, dem)
tgs <- project(tgs, dem)
ppt <- project(ppt, dem)

brk <- c(solar, toi, toip, toin, toip100, toin100, twi, tpi, dem, slope, aspect, slope.500, bt, tgs, ppt)

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
s.mu <- s %>% group_by(lmapunitiid=lmapunitiid) %>% summarise(
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
library(dplyr)
library(ggplot2)
library(Hmisc)
library(MASS)
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
hmnf.chm$hilly <- (hmnf.chm$slope.500 > 0.05)*1
hmnf.chm$toip100

c.table <- as.data.frame(cor(hmnf.chm[,c("chm.max", "chm.mean","age","solar","wet","hilly",
                                         "toi","toip","toin","toip100","toin100","twi","tpi","dem","slope","south","west","southwest","southeast","bt","tgs",
"ppt","T150_AWC","T50_sand","T150_sand","T50_clay","T150_clay","T50_OM","spodic","spodosols","Bhs",
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


soilsubset <- subset(hmnf.chm, T150_AWC < 15 & spodosols > 0.5 & Water_Table > 150)
soilsubset <- subset(soilsubset, slope.500 >= 0.05 & slope >= atan(0.0))

jackpine <- subset(soilsubset, jackpine %in% 1)
redpine <- subset(soilsubset, redpine %in% 1)
whitepine <- subset(soilsubset, whitepine %in% 1)
aspen <- subset(soilsubset, aspen %in% 1)
oak <- subset(soilsubset, oak %in% 1)
maple <- subset(soilsubset, maple %in% 1)




ggplot()+
  geom_smooth(aes(x=toip, y=chm.max, weight = age, col='maple'), data = maple)+
  geom_smooth(aes(x=toip, y=chm.max, weight = age, col='oak'), data = oak)+
  geom_smooth(aes(x=toip, y=chm.max, weight = age, col='aspen'), data = aspen)+
  geom_smooth(aes(x=toip, y=chm.max, weight = age, col='white pine'), data = whitepine)+
  geom_smooth(aes(x=toip, y=chm.max, weight = age, col='red pine'), data = redpine)+
  geom_smooth(aes(x=toip, y=chm.max, weight = age, col='jack pine'), data = jackpine)+
  scale_color_manual(name='forest type',  values = c('jack pine'='green','red pine'='darkgreen','white pine'='darkcyan',
                                   'aspen'='khaki4','oak'='orange','maple'='red'))

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

rp <- rpart(chm~
        toi2+toi+toip+toin+dem+d1230, data = port)
rpart.plot::rpart.plot(rp)


rp <- rpart(chm~
              dem+toi+toip+toin+tpi+slope+aspect, data = kat)
rpart.plot::rpart.plot(rp)