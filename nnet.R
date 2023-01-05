setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
library(neuralnet)
library(dplyr)
library(ranger)

brk.pts <- readRDS('output/brk.pts.RDS')
hmnf.chm2 <- readRDS('output/hmnf.chm3.RDS')


colnames(hmnf.chm2)
trainset = hmnf.chm2[sample(rownames(hmnf.chm2),10000),]
testingrows <- sample(rownames(trainset), floor(nrow(trainset)*0.1))

train <- trainset[!rownames(trainset) %in% testingrows,]
test <- trainset[rownames(trainset) %in% testingrows,]

rf <- ranger(chm.max~
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
             , data=train,
             num.trees=500, sample.fraction = 0.2,
             # max.depth = 5,  
             # mtry=7,
             write.forest = TRUE, importance = 'impurity'
)


nn <- neuralnet(chm.max~
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
               , data=train)



test <- test %>% mutate(rfpred = predictions(predict(rf, data=test)),
                        nnpred = predict(nn, newdata=test)[1])
nnpred = (predict(nn, newdata=test))