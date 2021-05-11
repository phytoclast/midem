library(stringr)
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
wd <- 'Michigan'
fn <- 'data.txt'
tif.path <- paste0('data/', wd,'/tif')
filelist <- paste0('data/', wd,'/',fn)
urls <- read.delim(filelist, header = F)
urls <- subset(urls, grepl('.tif$|.TIF$',V1))
urls$fcount <- str_count(urls$V1, '/')
urls$fn <-  str_split_fixed(urls$V1, '/', urls$fcount+1)[,urls$fcount+1]
if(!dir.exists(tif.path)){dir.create(tif.path)}
for (i in 1:nrow(urls)){
  download.file(urls$V1[i], paste0('data/', wd,'/tif/', urls$fn[i]), method = 'curl')
}

