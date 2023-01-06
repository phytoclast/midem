setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
library(neuralnet)
library(dplyr)
library(ranger)

brk.pts <- readRDS('output/brk.pts.RDS')
hmnf.chm2 <- readRDS('output/hmnf.chm3.RDS')
inputcols <-   c("chm.max","solar","hilly","age","pred1","pred2","toip","toin","bt","tgs","ppt","p1","p2","p3","p4","T150_AWC","T50_sand","T150_sand","T50_clay","T150_clay","T50_OM","T150_OM","T50_pH","Water_Table","wet","spodic","carbdepth","Bhs","jackpine","redpine","whitepine","aspen","oak","maple")

# inputcols2 <- c(solar,hilly,age,pred1,pred2,toip,toin,bt,tgs,ppt,p1,p2,p3,p4,T150_AWC,T50_sand,T150_sand,T50_clay,T150_clay,T50_OM,T150_OM,T50_pH,Water_Table,wet,spodic,carbdepth,Bhs,jackpine,redpine,whitepine,aspen,oak,maple)

hmnf.chm2.norm <- hmnf.chm2 %>% mutate(across(inputcols, ~ scale(.)[, 1]))

trainset = hmnf.chm2.norm[sample(rownames(hmnf.chm2.norm),10000),]
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
               , data=train,  hidden=c(2,2), threshold = 0.2, act.fct	= 'tanh')

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


#regression-enhanced random forest and ensemble models incorporating linear models?
