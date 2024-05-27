

basePath <- '/root/autodl-tmp/humPoulResult/data/'

world.map <- rnaturalearth::ne_countries(returnclass = "sf") |>filter(continent != "Antarctica")
globalCountry <- vect(world.map) 
coast <- ne_coastline(scale = "small", returnclass = "sf")
crs <- '+proj=longlat +datum=WGS84'

##########AE calculation###############
allDf <- fread(paste0(basePath,'allDf786_reclass.csv'))
speciesPixelNumPath <- list.files('/root/autodl-tmp/humPoulResult/data/single_model',pattern = '.tif',full.names = T)
spName <- basename(speciesPixelNumPath) %>% str_sub(.,1,-5)                      
speciesPixelNumPath2 <- speciesPixelNumPath[spName%in%allDf$LatName]


#######1.Cumulative species months#########
speciesPixelNum2 <- rast(speciesPixelNumPath2)
allMonth <- sum(speciesPixelNum2,na.rm=T)%>% mask(globalCountry)
#writeRaster(allMonth,paste0(basePath,'AE_data/allMonth.tif'), overwrite=T)

#######2.AE##########
allMonth <- rast(paste0(basePath,'AE_data/allMonth.tif'))
calEntropy <- lapply(speciesPixelNumPath2, function(x){
  r <- rast(x) %>% sum(.,na.rm=T)
  pi <- r/allMonth
  y <- -pi*log(pi)
  names(y) <- str_sub(basename(x),1,-5)
  return(y)
})
calEntropy <- rast(calEntropy)
AE<-sum(calEntropy,na.rm = T)%>% mask(globalCountry)
plot(AE)
#writeRaster(AE,paste0(basePath,'AE_data/AE.tif'), overwrite=T)



#AE-Functional Group#############



#AE validation############
outBreak <- fread('/root/autodl-tmp/YANZHENG/point/allData.csv')
#outBreak$Latitude <- as.numeric(outBreak$Latitude)
outBreak$Longitude <- as.numeric(outBreak$Longitude)
crs <- '+proj=longlat +datum=WGS84'
globalRaster <- rast(vals=1:259200,nrows=360, ncols=720,xmin=-180, xmax=180,ymin=-90, ymax=90,crs=crs)
addGeom <- cellFromXY(globalRaster,outBreak[,c('Longitude','Latitude')]) %>% 
  cbind(id=.,outBreak)
thinData <- unique(data.table(addGeom),by='id') %>% dplyr::select(.,-id)
outBreak2 <- vect(thinData,geom=c('Longitude','Latitude'),crs=crs)
plot(outBreak2)

world.map <- rnaturalearth::ne_countries(returnclass = "sf") |>filter(continent != "Antarctica")
globalCountry <- vect(world.map) 
#overEntropyNum2 <- rast('/root/autodl-tmp/zyresult/EntropySum(868版).tif')%>% mask(globalCountry)
overEntropyNum2 <- rast('/root/autodl-tmp/zyresult/Entropy/Global_birds_Entropy.tif')%>% mask(globalCountry)
global(overEntropyNum2,quantile,probs=seq(0, 1, 0.05),na.rm=T)

overEntropyNum3 <- ifel(overEntropyNum2>4,overEntropyNum2,NA)

# ggplot()+
#   geom_spatraster(data=overEntropyNum2,na.rm=T)+
#   scale_fill_gradientn(colors=paletteer_c("grDevices::Zissou 1",30))+
#   geom_spatvector(data=outBreak2,fill=NA,na.rm=T)

df <- terra::extract(overEntropyNum2,outBreak2,mean,bind=T,exact=T,touches=T,na.rm=T,cells=T)%>% as.data.frame()
allDf <- terra::as.data.frame(overEntropyNum2,cells=T)
threshold <- 4.107

TP <- sum(df$sum>=threshold,na.rm = T)  ;TP
FN <- sum(df$sum<threshold,na.rm = T)   ;FN
FP <- sum(allDf$sum>=threshold,na.rm = T)-TP  ;FP
TN <- nrow(allDf)-TP-FN-FP  ;TN

accuracy <- (TP+TN)/(TP+TN+FN+FP)  ;accuracy
precise <- (TP)/(TP+FP)  ;precise
recall <- (TP)/(TP+FN)  ;recall
FI <- 2*((recall*precise)/(recall+precise))  ;FI

library(pROC)

df$label <- 'label'
rocDf <- left_join(allDf,df[,c('cell','label')])
rocDf$label <- ifelse(is.na(rocDf$label),'poor',rocDf$label)
dfroc1<- roc(rocDf$label, rocDf$sum)
plot(dfroc1,col="red",#颜色
     legacy.axes=T,#y轴格式更改
     print.auc=TRUE,#显示AUC面积
     expand=c(0,0),
     print.thres=TRUE,#添加截点和95%CI
     grid=c(0.2,0.2),grid.col=c("lightgrey","lightgrey"),
     cex.main=1.5,  # 主标题字体放大1.5倍
     cex.sub=1.5,   # 副标题字体放大1.5倍
     cex.axis=1.5,  # 坐标轴刻度字体放大1.5倍
     cex.lab=1.5)#网格线设置









