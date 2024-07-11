
library(terra)
library(data.table)
library(dplyr)
library(stringr)
library(tidyterra)

basePath <- '/root/autodl-tmp/humPoulResult/data/'
#basePath <-'/root/autodl-tmp/result2'
world.map <- rnaturalearth::ne_countries(returnclass = "sf") |>filter(continent != "Antarctica")
globalCountry <- vect(world.map) 
coast <- ne_coastline(scale = "small", returnclass = "sf")
crs <- '+proj=longlat +datum=WGS84'

##########AE calculation###############
# allDf <- fread(paste0(basePath,'allDf786_reclass.csv'))
# speciesPixelNumPath <- list.files('/root/autodl-tmp/humPoulResult/data/single_model',pattern = '.tif',full.names = T)
# spName <- basename(speciesPixelNumPath) %>% str_sub(.,1,-5)                      
# speciesPixelNumPath2 <- speciesPixelNumPath[spName%in%allDf$LatName]
# length(speciesPixelNumPath2)
# 
# #######1.Cumulative species months#########
# speciesPixelNum2 <- rast(speciesPixelNumPath2)
# allMonth <- sum(speciesPixelNum2,na.rm=T)%>% mask(globalCountry)
# #writeRaster(allMonth,paste0(basePath,'AE_data/allMonth.tif'), overwrite=T)
# 
# #######2.AE##########
# allMonth <- rast(paste0(basePath,'AE_data/allMonth.tif'))
# calEntropy <- lapply(speciesPixelNumPath2, function(x){
#   r <- rast(x) %>% sum(.,na.rm=T)
#   pi <- r/allMonth
#   y <- -pi*log(pi)
#   names(y) <- str_sub(basename(x),1,-5)
#   return(y)
# })
# calEntropy <- rast(calEntropy)
# AE<-sum(calEntropy,na.rm = T)%>% mask(globalCountry)
# plot(AE)
#writeRaster(AE,paste0(basePath,'AE_data/AE.tif'), overwrite=T)



#new WAE 0618按照分月份thin###########
# speciesPre<-rast('/root/autodl-tmp/result2/presenceRast.tif')
# 
# allDf <- fread(paste0(basePath,'allDf786_reclass.csv'))
# 
# matching_indices <- which(names(speciesPre) %in% allDf$LatName)
# length(matching_indices)
# 
# layers_list <- lapply(1:nlyr(speciesPre), function(i) speciesPre[[i]])
# 
# layers_list786 <- layers_list[matching_indices]
# combined_raster <- rast(layers_list786)
# #allMonth <- sum(combined_raster,na.rm=T)%>% mask(globalCountry)


allDf <- fread(paste0(basePath,'allDf786_reclass.csv'))
speciesPixelNumPath <- list.files('/root/autodl-tmp/SDMsingle',pattern = '.tif',full.names = T)
spName <- basename(speciesPixelNumPath) %>% str_sub(.,1,-5)                      
speciesPixelNumPath2 <- speciesPixelNumPath[spName%in%allDf$LatName]
length(speciesPixelNumPath2)

#######1.Cumulative species months#########
speciesPixelNum2 <- rast(speciesPixelNumPath2)
allMonth <- sum(speciesPixelNum2,na.rm=T)%>% mask(globalCountry)
#writeRaster(allMonth,paste0(basePath,'AE_data/allMonth.tif'), overwrite=T)

#######2.AE##########
#allMonth <- rast(paste0(basePath,'AE_data/allMonth.tif'))
calEntropy <- lapply(speciesPixelNumPath2, function(x){
  r <- rast(x) %>% sum(.,na.rm=T)
  pi <- r/allMonth
  y <- -pi*log(pi)
  names(y) <- str_sub(basename(x),1,-5)
  return(y)
})
calEntropy <- rast(calEntropy)
#calEntropy <- clamp(calEntropy, lower=0)# 将所有小于0的值设置为0: AE0619先给分熵的负值赋予0
AE<-sum(calEntropy,na.rm = T)%>% mask(globalCountry); plot(AE)
#AE0620 <- clamp(AE, lower=0)# 将所有小于0的值设置为0: AE0620先给分熵的负值赋予0
AE0620<-AE
plot(AE0620)
#writeRaster(AE, '/root/autodl-tmp/result2/AE0619.tif', overwrite=T)
#writeRaster(AE0620, '/root/autodl-tmp/result2/AE0620.tif', overwrite=T)

#对比
AE0620<-rast('/root/autodl-tmp/result2/AE0620.tif')
plot(AE0620)
AE_used<-rast(paste0(basePath,'AE_data/AE.tif'))
plot(AE_used)


AE_com<-AE0620-AE_used
AE_com1 <- ifel(AE0620>4.107,4,0)

plot(AE_com1)

raster_values <- values(calEntropy)
negative_layers <- which(apply(raster_values, 2, function(x) any(x < 0)))# 检查哪些栅格层包含负值
print(negative_layers)# 打印出包含负值的栅格层索引

AE_A<-calEntropy$`Actitis hypoleucos`
AE_A1 <- ifel(calEntropy$`Actitis hypoleucos`>0,1,0)
plot(AE_A1)

library(terra)
library(data.table)
library(dplyr)
library(stringr)
library(tidyterra)
#######Fig2. WAE validation#########################
`%notin%` <- Negate(`%in%`)
crs <- '+proj=longlat +datum=WGS84'
globalRaster <- rast(vals=1:259200,nrows=360, ncols=720,xmin=-180, xmax=180,ymin=-90, ymax=90,crs=crs)
globalSHP <- vect('/root/autodl-tmp/worldBorder/ne_10m_admin_0_countries.shp')
# calculate  threshold
outBreakData <- fread('/root/autodl-tmp/YANZHENG/point/allData.csv')%>% subset(Diagnosis.status=='Confirmed')
outBreakData2 <- na.omit(outBreakData[,c('Country','Longitude','Latitude')])
addGeom <- cellFromXY(globalRaster,outBreakData2[,c('Longitude','Latitude')]) %>% 
  cbind(id=.,outBreakData2)
thinData <- unique(data.table(addGeom),by='id') %>% dplyr::select(.,-id)
breakVect <- vect(thinData,geom=c('Longitude','Latitude'),crs=crs)
breakRast <- rasterize(breakVect,globalRaster)

TP_points <- vect(outBreakData2,geom=c('Longitude','Latitude'),crs=crs)
AE0619<-rast('/root/autodl-tmp/result2/AE0619.tif')
sample_TP <- terra::extract(AE0619,TP_points,mean,bind=T) %>% as.data.frame()
sample_TP$label <- 'good'

trueNgRast <- mask(AE0619,breakRast,inverse=T)

set.seed(1234)
sample_TN<- terra::spatSample(trueNgRast,dim(sample_TP)[1],method='random',exhaustive=T,na.rm=T)
sample_TN$label <- 'poor'
sampleData <- rbind(sample_TN,sample_TP[,c('sum','label')])
dfroc1<- roc(sampleData$label, sampleData$sum)
dfroc1$auc

#litter
plot(dfroc1,col="red",
     legacy.axes=T,
     print.auc=TRUE,
     expand=c(0,0),
     print.thres=TRUE,
     grid=c(0.2,0.2),grid.col=c("lightgrey","lightgrey"),
     cex.main=1.5, 
     cex.sub=1.5,  
     cex.axis=1.5, 
     cex.lab=1.5)


#a##########
threshold <- 4.2

WAE<- rast(paste0(basePath,'AE_data/AE.tif'))%>% mask(globalCountry)#used
#WAE<-rast('/root/autodl-tmp/result2/AE0620.tif')

plotWAE <- ifel(WAE>=threshold,1,0)   #4.19

coast <- rnaturalearth::ne_coastline(scale = "small", returnclass = "sf")
crs <- '+proj=longlat +datum=WGS84'
ggplot() +
  geom_spatraster(data = plotWAE) +
  geom_spatvector(data=coast,fill=NA)+
  #coord_sf(crs = crs,xlim=c(-160,180),ylim=c(-56,90))+
  geom_spatvector(data=TP_points,size=1, shape = 4,color='#ff4f4c',fill='#f26c6a', alpha = 0.5,stroke = 0.3)+   #shape = 21是圆圈，shape = 4是×
  theme_bw()+
  #scale_color_manual(values='#fc4e4e', labels='Outbreak of\navian influenza') +
  scale_fill_gradient(low = "lightgrey",high = "yellow" ,space = "Lab",n.break=2,
                      labels=c('Low risk','High risk'),na.value='white')+
  guides(
    fill=guide_legend()  )+
  labs(title = NULL)+
  theme(
    text = element_text(size = 15),
    #axis.text = element_text(size = 15),
    plot.title = element_text(hjust=0.5),
    # axis.line = element_line(color = 'black'),
    legend.title = element_blank(),
    legend.text = element_text(size = 15),
    legend.position = "none",
    legend.direction='horizontal',
    legend.key.width = unit(0.4,'cm'), 
    legend.key.height = unit(3,'cm')
  )







# #AE validation############
# outBreak <- fread('/root/autodl-tmp/YANZHENG/point/allData.csv')
# #outBreak$Latitude <- as.numeric(outBreak$Latitude)
# outBreak$Longitude <- as.numeric(outBreak$Longitude)
# crs <- '+proj=longlat +datum=WGS84'
# globalRaster <- rast(vals=1:259200,nrows=360, ncols=720,xmin=-180, xmax=180,ymin=-90, ymax=90,crs=crs)
# addGeom <- cellFromXY(globalRaster,outBreak[,c('Longitude','Latitude')]) %>% 
#   cbind(id=.,outBreak)
# thinData <- unique(data.table(addGeom),by='id') %>% dplyr::select(.,-id)
# outBreak2 <- vect(thinData,geom=c('Longitude','Latitude'),crs=crs)
# plot(outBreak2)
# 
# world.map <- rnaturalearth::ne_countries(returnclass = "sf") |>filter(continent != "Antarctica")
# globalCountry <- vect(world.map) 
# #overEntropyNum2 <- rast('/root/autodl-tmp/zyresult/EntropySum(868版).tif')%>% mask(globalCountry)
# overEntropyNum2 <- rast('/root/autodl-tmp/zyresult/Entropy/Global_birds_Entropy.tif')%>% mask(globalCountry)
# global(overEntropyNum2,quantile,probs=seq(0, 1, 0.05),na.rm=T)
# 
# overEntropyNum3 <- ifel(overEntropyNum2>4,overEntropyNum2,NA)
# 
# # ggplot()+
# #   geom_spatraster(data=overEntropyNum2,na.rm=T)+
# #   scale_fill_gradientn(colors=paletteer_c("grDevices::Zissou 1",30))+
# #   geom_spatvector(data=outBreak2,fill=NA,na.rm=T)
# 
# df <- terra::extract(overEntropyNum2,outBreak2,mean,bind=T,exact=T,touches=T,na.rm=T,cells=T)%>% as.data.frame()
# allDf <- terra::as.data.frame(overEntropyNum2,cells=T)
# threshold <- 4.107
# 
# TP <- sum(df$sum>=threshold,na.rm = T)  ;TP
# FN <- sum(df$sum<threshold,na.rm = T)   ;FN
# FP <- sum(allDf$sum>=threshold,na.rm = T)-TP  ;FP
# TN <- nrow(allDf)-TP-FN-FP  ;TN
# 
# accuracy <- (TP+TN)/(TP+TN+FN+FP)  ;accuracy
# precise <- (TP)/(TP+FP)  ;precise
# recall <- (TP)/(TP+FN)  ;recall
# FI <- 2*((recall*precise)/(recall+precise))  ;FI
# 
# library(pROC)
# 
# df$label <- 'label'
# rocDf <- left_join(allDf,df[,c('cell','label')])
# rocDf$label <- ifelse(is.na(rocDf$label),'poor',rocDf$label)
# dfroc1<- roc(rocDf$label, rocDf$sum)
# plot(dfroc1,col="red",#颜色
#      legacy.axes=T,#y轴格式更改
#      print.auc=TRUE,#显示AUC面积
#      expand=c(0,0),
#      print.thres=TRUE,#添加截点和95%CI
#      grid=c(0.2,0.2),grid.col=c("lightgrey","lightgrey"),
#      cex.main=1.5,  # 主标题字体放大1.5倍
#      cex.sub=1.5,   # 副标题字体放大1.5倍
#      cex.axis=1.5,  # 坐标轴刻度字体放大1.5倍
#      cex.lab=1.5)#网格线设置
# 
# 
# 






