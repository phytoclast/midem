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
twi <- rast('D:/GIS/DEM/hmnftwi.tif'); names(twi) = 'twi'
tpi <- rast('D:/GIS/DEM/hmnftopographicpositionindex.tif'); names(tpi) = 'tpi'
dem <- rast('D:/GIS/DEM/hmnfdem30.tif'); names(dem) = 'dem'
slope <- rast('D:/GIS/DEM/hmnfslope.tif'); names(slope) = 'slope'
aspect <- rast('D:/GIS/DEM/hmnfaspect.tif'); names(aspect) = 'aspect'
bt <- rast('D:/scripts/snow/output/bt.90.alt.tif'); names(bt) = 'bt'
tgs <- rast('D:/scripts/snow/output/tgs.90.alt.tif'); names(tgs) = 'tgs'
ppt <- rast('C:/a/Ecological_Sites/GIS/Climate/PRISM2010/P/p0112/w001001.adf'); names(ppt) = 'ppt'

bt <- project(bt, dem)
tgs <- project(tgs, dem)
ppt <- project(ppt, dem)

brk <- c(solar, twi, tpi, dem, slope, aspect, bt, tgs, ppt)

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
s.mu <- s %>% group_by(lmapunitiid=lmapunitiid) %>% summarise(
  T150_AWC = wtd.mean(T150_AWC, weights=comppct_r, na.rm=T),
  T50_sand = wtd.mean(T50_sand, weights=comppct_r, na.rm=T),
  T150_sand = wtd.mean(T150_sand, weights=comppct_r, na.rm=T),
  T50_clay = wtd.mean(T50_clay, weights=comppct_r, na.rm=T),
  T150_clay = wtd.mean(T150_clay, weights=comppct_r, na.rm=T),
  T50_OM = wtd.mean(T50_OM, weights=comppct_r, na.rm=T),
  T150_OM = wtd.mean(T150_OM, weights=comppct_r, na.rm=T),
  T50_pH = wtd.mean(T50_pH, weights=comppct_r, na.rm=T),
  Water_Table = wtd.mean(Water_Table, weights=comppct_r, na.rm=T)
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
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

hmnf.chm <- read.csv('output/hmnf.chm2.csv')

hmnf.chm <- subset()
c.table <- as.data.frame(cor(hmnf.chm[,c("chm.max", "chm.mean","age","solar",
"twi","tpi","dem","slope","aspect","bt","tgs",
"ppt","T150_AWC","T50_sand","T150_sand","T50_clay","T150_clay","T50_OM",
"T150_OM","T50_pH","Water_Table")], use='complete.obs'
))

write.csv(c.table, 'output/c.table.csv', row.names = F)

model <- lm(chm.max~age+solar+
            twi+tpi+dem+slope+aspect+bt+tgs+
            ppt+T150_AWC+T50_sand+T150_sand+T50_clay+T150_clay+T50_OM+
            T150_OM+T50_pH+Water_Table, data=hmnf.chm)

summary(model)