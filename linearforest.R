library(grf)
library(glmnet)
library(ggplot2)
library(dplyr)
library(ranger)
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

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

mergedset.norm <- mergedset# %>% mutate(across(inputcols, ~ scale(.)[, 1]))
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
Y <- train$chm.max
X <- train %>% select(solar,
                          age,
                          pred1,
                          pred2,
                          toip,
                          toin,
                          bt,
                          tgs,
                          ppt,
                          p1,
                          p2,
                          p3,
                          p4,
                          T150_AWC,
                          T50_sand,
                          T150_sand,
                          T50_clay,
                          T150_clay,
                          T50_OM,
                          T150_OM,
                          T50_pH,
                          Water_Table,
                          wet,
                          spodic,
                          carbdepth,
                          Bhs,
                          jackpine,
                          redpine,
                          whitepine,
                          aspen,
                          oak,
                          maple)

test.X <- test %>% select(solar,
                           age,
                           pred1,
                           pred2,
                           toip,
                           toin,
                           bt,
                           tgs,
                           ppt,
                           p1,
                           p2,
                           p3,
                           p4,
                           T150_AWC,
                           T50_sand,
                           T150_sand,
                           T50_clay,
                           T150_clay,
                           T50_OM,
                           T150_OM,
                           T50_pH,
                           Water_Table,
                           wet,
                           spodic,
                           carbdepth,
                           Bhs,
                           jackpine,
                           redpine,
                           whitepine,
                           aspen,
                           oak,
                           maple)



rf <- ranger(formular
             , data=train,
             num.trees=100, sample.fraction = 0.5,
             # max.depth = 5,  
             # mtry=7,
             write.forest = TRUE, importance = 'impurity'
)

lf <- grf::ll_regression_forest(X,Y, num.trees = 500, enable.ll.split = T)

results.llf <- predict(lf, test.X)
preds.llf <- results.llf$predictions

test <- test %>% mutate(pred.rf = predictions(predict(rf, test)),
                        pred.lf = preds.llf
                        )


mean((test$pred.rf - test$chm.max)^2)^0.5
mean((test$pred.lf - test$chm.max)^2)^0.5




####
####
####
n <- 1000
p <- 20
X <- matrix(runif(n * p, 0, 1), nrow = n)


X[,2] <- sin(X[,1]*50)
X[,3] <- 1/((X[,1])+1)
X[,4] <- tan(X[,1])
X[,5] <- (X[,1])^3
X[,6] <- (X[,1])^0.5
X[,7] <- cos(X[,1]*20)
Y <- X[,1]*0.1+X[,2]*.5+X[,3]*10+X[,4]*.5+X[,5]*1+X[,6]*1+X[,7]*0.6+X[,8]*0.2

X[,2] <- (X[,2]*5 + runif(n, 0, 1)*1)/2
X[,3] <- (X[,3]*5 + runif(n, 0, 1)*1)/2
X[,4] <- (X[,4]*5 + runif(n, 0, 1)*1)/2
X[,5] <- (X[,5]*5 + runif(n, 0, 1)*1)/2
X[,6] <- (X[,6]*5 + runif(n, 0, 1)*1)/2
X[,7] <- (X[,7]*5 + runif(n, 0, 1)*1)/2





trainset <- sample(1:nrow(X), 200)
testset <- sample(1:nrow(X), 100)
X.train <- X[trainset,]
Y.train <- Y[trainset]
removeX <- which(X.train[,1] < 0.5 | X.train[,1] > 0.7)
X.train <- X.train[removeX,]
Y.train <- Y.train[removeX]
X.test <- X[testset,]
Y.test <- Y[testset]
g1 <- ggplot() +
  geom_point(aes(x=X[,1], y = Y)) +
  geom_point(aes(x=X.train[,1], y = Y.train), color='red') +
  xlab("x") + ylab("y") + theme_bw()
g1


forest <- regression_forest(X.train, Y.train, honesty = F)
mod.lm <- lm(Y.train ~ .,  data=as.data.frame(X.train))
summary(mod.lm)
mod.rf <- ranger(Y.train ~ .,  data=as.data.frame(X.train),
             sample.fraction = 1, num.trees = 2000)
ll.forest <- ll_regression_forest(X.train, Y.train, enable.ll.split = F, honesty = F)
b.forest <- boosted_regression_forest(X.train, Y.train, honesty = F)

preds.forest <- predict(forest, X.test)$predictions
preds.rf <- predict(mod.rf, as.data.frame(X.test))$predictions
preds.lm <- predict(mod.lm, as.data.frame(X.test))
preds.llf <- predict(ll.forest, X.test,
                     linear.correction.variables = 1)$predictions
preds.bf <- predict(b.forest, X.test)$predictions

g1 <- ggplot() +
  geom_line(aes(x=X[,1], y = Y), color='black') +
  #geom_point(aes(x=X.train[,1], y = Y.train), color='red') +
  geom_point(aes(x=X.test[,1], y = preds.lm), color='red') +
  geom_point(aes(x=X.test[,1], y = preds.forest), color='blue') +
  geom_point(aes(x=X.test[,1], y = preds.bf), color='yellow') +
  geom_point(aes(x=X.test[,1], y = preds.llf), color='green') 
# +
#   xlab("x") + ylab("y") + theme_bw()
g1

# mean((Y.test - preds.forest)^2)
# mean((Y.test - preds.rf)^2)
# mean((Y.test - preds.lm)^2)
# mean((Y.test - preds.llf)^2)
# mean((Y.test - mean(Y.test))^2)


#preds.forest
mean((Y.test - preds.forest)^2)/mean((Y.test - mean(Y.test))^2)*100
#preds.rf
mean((Y.test - preds.rf)^2)/mean((Y.test - mean(Y.test))^2)*100
#preds.lm
mean((Y.test - preds.lm)^2)/mean((Y.test - mean(Y.test))^2)*100
#preds.llf
mean((Y.test - preds.llf)^2)/mean((Y.test - mean(Y.test))^2)*100
#preds.bf
mean((Y.test - preds.bf)^2)/mean((Y.test - mean(Y.test))^2)*100


