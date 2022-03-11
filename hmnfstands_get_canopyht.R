library(sf)
library(terra)
library(raster)
library(minpack.lm)
library(growthmodels)

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

hmnf <- read.csv('data/hmnfstandpoints/SiteIndexRandomPoints_10perStand.txt')


chm1 <- rast('output/NL/MLRA96.tif')
chm2 <- rast('output/NL/MLRA94AC.tif')

chm <- merge(chm1, chm2)
chm[chm<0] <- 0
chm[chm>50] <- 50

chm.30 <- aggregate(chm, fact=30/5, na.rm=T, fun='max')


hmnf <- st_as_sf(hmnf, coords=c(x='Long', y='Lat'), crs='epsg:4326')
hmnf.mi <- sf::st_transform(hmnf, crs=crs(chm))
hmnf.chm <- terra::extract(chm.30, vect(hmnf.mi))
colnames(hmnf.chm) <- c('id','chm')
hmnf.chm1 <- cbind(hmnf.mi, subset(hmnf.chm, chm=hmnf.chm[,2]))

hmnf.chm1$age <- 2021-hmnf.chm1$YEAR_OF_ORIGIN

library(ggplot2)
ggplot(data=hmnf.chm1)+
geom_point(aes(x=age, y=chm), col='red', alpha=0.02)

redpine <- subset(hmnf.chm1, EV_CODE %in% 2 & !is.na(chm) & !is.na(age))
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

redpine <- subset(hmnf.chm1, EV_CODE %in% 2 & !is.na(chm) & !is.na(age) & age >= 5 & chm >=0.1)
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
redpine1 <- subset(hmnf.chm1, EV_CODE %in% 1 & !is.na(chm) & !is.na(age) & age >= 5 & chm >=0.1 & 
                     SITE_INDEX >= mean(redpine$SITE_INDEX+1))
redpine2 <- subset(hmnf.chm1, EV_CODE %in% 1 & !is.na(chm) & !is.na(age) & age >= 5 & chm >=0.1 & 
                     SITE_INDEX < mean(redpine$SITE_INDEX)-1)
x1 <- redpine1$age
y1 <- redpine1$chm
x1=c(x1,x1*0-1)
y1=c(y1,y1*0)
x2 <- redpine2$age
y2 <- redpine2$chm
x2=c(x2,x2*0-1)
y2=c(y2,y2*0)
mod1 <- nlsLM(y1 ~ growthmodels::gompertz(x1, alpha, beta, k), start = list(alpha = 30, beta = .1, k = .1))
summary(mod1)
mod2 <- nlsLM(y2 ~ growthmodels::gompertz(x2, alpha, beta, k), start = list(alpha = 30, beta = .1, k = .1))
summary(mod2)

x0 <- 0:200+0
pred1 <-  predict(mod1, list(x1 = x0))
pred2 <-  predict(mod2, list(x2 = x0))

ggplot()+
  geom_point(aes(x=x1, y=y1), col='red', alpha=0.02)+
  geom_point(aes(x=x2, y=y2), col='blue', alpha=0.02)+
  geom_line(aes(x=x0,y=pred1), col='red', size=2)+
  geom_line(aes(x=x0,y=pred2), col='blue', size=2)
  