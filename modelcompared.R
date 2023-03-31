
library(ggplot2)
library(terra)

maple.soil <- rast('D:/scripts/midem/output/maple.si.soil.tif') |> as.data.frame() |> data.frame(model='soil')

maple.si  <- rast('D:/scripts/midem/output/maple.si.tif') |> as.data.frame()  |> data.frame(model='hybrid')

maple.si.topo  <- rast('D:/scripts/midem/output/maple.si.topo.tif') |> as.data.frame()  |> data.frame(model='topography')

maple.si.lm  <- rast('D:/scripts/midem/output/maple.si.lm.tif') |> as.data.frame()  |> data.frame(model='lm')

maple.si.rf  <- rast('D:/scripts/midem/output/maple.si.rf.tif') |> as.data.frame()  |> data.frame(model='rf')

maple.si.topo.clim  <- rast('D:/scripts/midem/output/maple.si.topo.clim.tif') |> as.data.frame()  |> data.frame(model='no soil')


colnames(maple.soil) <- c('chm','model')
colnames(maple.si) <- c('chm','model')
colnames(maple.si.topo) <- c('chm','model')
colnames(maple.si.lm) <- c('chm','model')
colnames(maple.si.rf) <- c('chm','model')
colnames(maple.si.topo.clim) <- c('chm','model')

df <- rbind(maple.soil, maple.si, maple.si.topo, maple.si.lm, maple.si.rf, maple.si.topo.clim)
df <- df |> mutate(model = factor(model, levels = c("lm","rf","hybrid","soil","no soil","topography")))
ggplot(df, aes(x=model, y=chm))+
  geom_boxplot()


df.sd <- df |> group_by(model) |> summarise(sd = sd(chm), mean = mean(chm))
write.csv(df.sd, 'D:/scripts/midem/output/df.sd.csv', row.names = F)
