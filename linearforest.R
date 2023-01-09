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
lf <- grf::lm_forest(formular
             , data=train,
             num.trees=500, sample.fraction = 0.2,
             # max.depth = 5,  
             # mtry=7,
             write.forest = TRUE, importance = 'impurity'
)


####
####
####
p <- 20
n <- 1000
sigma <- sqrt(20)

mu <- function(x){ log(1 + exp(6 * x)) }
X <- matrix(runif(n * p, -1, 1), nrow = n)
Y <- mu(X[,1]) + sigma * rnorm(n)

X.test <- matrix(runif(n * p, -1, 1), nrow = n)
ticks <- seq(-1, 1, length = n)
X.test[,1] <- ticks
truth <- mu(ticks)

forest <- regression_forest(X, Y)
preds.forest <- predict(forest, X.test)$predictions

df <- data.frame(cbind(ticks, truth, preds.forest))
g1 <- ggplot(df, aes(ticks)) +
  geom_point(aes(y = preds.forest, color = "Regression Forest"), show.legend = F, size = 0.6) +
  geom_line(aes(y = truth)) +
  xlab("x") + ylab("y") + theme_bw()
g1

ll.forest <- ll_regression_forest(X, Y, enable.ll.split = TRUE)
preds.llf <- predict(ll.forest, X.test,
                     linear.correction.variables = 1)$predictions

df.llf <- data.frame(cbind(ticks, truth, preds.llf))
g2 <- ggplot(df.llf, aes(ticks)) +
  geom_point(aes(y = preds.llf, color = "Local Linear Forest"), show.legend = F, size = 0.6) +
  geom_line(aes(y = truth)) +
  xlab("x") + ylab("y") + theme_bw()
g2



