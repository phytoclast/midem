setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
library(neuralnet)
library(dplyr)
library(ranger)

brk.pts <- readRDS('output/brk.pts.RDS')
hmnf.chm2 <- readRDS('output/hmnf.chm3.RDS')
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

brk.pts <- brk.pts %>% mutate(age=50, pred1 = asmod(age = 50), pred2 =  rpmod(age = 50), jackpine = 0, redpine = 0, whitepine = 0, aspen = 0, maple = 1, oak = 0, status = 'raster')

hmnf.chm2$status <- 'training'

mergedset <- dplyr::bind_rows(brk.pts, hmnf.chm2)


inputcols <-   c("chm.max","solar","hilly","age","pred1","pred2","toip","toin","bt","tgs","ppt","p1","p2","p3","p4","T150_AWC","T50_sand","T150_sand","T50_clay","T150_clay","T50_OM","T150_OM","T50_pH","Water_Table","wet","spodic","carbdepth","Bhs","jackpine","redpine","whitepine","aspen","oak","maple")

# inputcols2 <- c(solar,hilly,age,pred1,pred2,toip,toin,bt,tgs,ppt,p1,p2,p3,p4,T150_AWC,T50_sand,T150_sand,T50_clay,T150_clay,T50_OM,T150_OM,T50_pH,Water_Table,wet,spodic,carbdepth,Bhs,jackpine,redpine,whitepine,aspen,oak,maple)

mergedset.norm <- mergedset %>% mutate(across(inputcols, ~ scale(.)[, 1]))
# mergedset.norm2 <- mergedset %>% select(all_of(inputcols)) %>% scale(.)[, 1]
meanchm <-  mean(hmnf.chm2$chm.max)
sdchm <-  sd(hmnf.chm2$chm.max)

trainset <- subset(mergedset.norm, status  %in% 'training')
trainset <- trainset %>% mutate(prob = ifelse(hilly > 0, 4,1)) 
trainset = trainset[sample(rownames(trainset),20000, prob=trainset$prob),]
testingrows <- sample(rownames(trainset), floor(nrow(trainset)*0.1))
train <- trainset[!rownames(trainset) %in% testingrows,]
test <- trainset[rownames(trainset) %in% testingrows,]
formular <- chm.max~
  solar+
  #hilly+
  age+
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

rf <- ranger(formular
             , data=train,
             num.trees=500, sample.fraction = 0.2,
             # max.depth = 5,  
             # mtry=7,
             write.forest = TRUE, importance = 'impurity'
)


nn <- neuralnet(formular
               , data=train,  hidden=c(2,2), threshold = 0.2, stepmax = 2e+05, act.fct	= 'tanh')

plot(nn)
nn$result.matrix





rf1 <- ranger(formular
              , data=train,
              num.trees=500, sample.fraction = 0.2,
              # max.depth = 5,  
              # mtry=7,
              write.forest = TRUE, importance = 'impurity'
)
rf2 <- ranger(formular
              , data=train,
              num.trees=500, sample.fraction = 0.5,
              # max.depth = 5,  
              # mtry=7,
              write.forest = TRUE, importance = 'impurity'
)
rf3 <- ranger(formular
              , data=train,
              num.trees=500, sample.fraction = 0.75,
              # max.depth = 5,  
              # mtry=7,
              write.forest = TRUE, importance = 'impurity'
)
rf4 <- ranger(formular
              , data=train,
              num.trees=500, sample.fraction = 1,
              # max.depth = 5,  
              # mtry=7,
              write.forest = TRUE, importance = 'impurity'
)


test <- test %>% mutate(rfpred1 = predictions(predict(rf1, data=test)),
                        rfpred2 = predictions(predict(rf2, data=test)),
                        rfpred3 = predictions(predict(rf3, data=test)),
                        rfpred4 = predictions(predict(rf4, data=test)),
                        nnpred = predict(nn, newdata=test))

resids <- test %>% mutate(res1 = (rfpred1 - chm.max)^2,
                          res2 = (rfpred2 - chm.max)^2,
                          res3 = (rfpred3 - chm.max)^2,
                          res4 = (rfpred4 - chm.max)^2,
                          resn = (nnpred - chm.max)^2)

mean(resids$res1)^0.5
mean(resids$res2)^0.5
mean(resids$res3)^0.5
mean(resids$res4)^0.5
mean(resids$resn)^0.5


mergedset.pred <- mergedset.norm %>% mutate(resn = (predict(nn, newdata=mergedset.norm))) %>% subset(status %in% 'raster') 


mergedset.pred <- mergedset.pred %>% mutate(newchm = resn*sdchm + meanchm)

library(terra)
dem <- rast('D:/GIS/DEM/hmnfdem30.tif'); names(dem) = 'dem'

nnmap <- rast(cbind(x=mergedset.pred$x, y=mergedset.pred$y, resn=mergedset.pred$newchm), type="xyz", crs=crs(dem))

plot(nnmap)
writeRaster(nnmap, 'output/nnmaple.tif', overwrite=T)
#regression-enhanced random forest and ensemble models incorporating linear models?
#

