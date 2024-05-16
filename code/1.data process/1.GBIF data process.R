#####GBIF data process####
#####2024/5/1##########
######author Qiaoyuxin Zhanyue########
library(data.table)
library(terra)
library(tidyterra)
library(stringr)
`%notin%` <- Negate(`%in%`)
#####global0.5 ° raster definition#############
crs <- '+proj=longlat +datum=WGS84'
globalRaster <- rast(vals=1:259200,nrows=360, ncols=720,xmin=-180, xmax=180,ymin=-90, ymax=90,crs=crs)

#####GBIF data###########
basePath <- '/root/autodl-tmp/GBIF_Data/'
allPath <- list.files(paste0(basePath,''),pattern='.csv',full.names=T)
####IUCN border######
borderList <- list.files('/root/autodl-tmp/IUCN_Data/IUCNsp/',pattern = '.shp',full.names = T)

########切割gbif数据,目转种并筛去小于108点的物种################
for (path in allPath) {
  birdData <- fread(path,select = c("decimalLatitude","decimalLongitude",'species','month','year','basisOfRecord','order',
                                    'family','genus'),quote="") %>% 
    subset(.,basisOfRecord%in%c('HUMAN_OBSERVATION','OBSERVATION','OCCURRENCE','MACHINE_OBSERVATION')) %>% 
    .[,basisOfRecord:=NULL]
  birdData$species <- ifelse(birdData$species=='',NA,birdData$species)
  birdData <- na.omit(birdData)
  colnames(birdData)[1:2] <- c('lat','lon')
  ####对每个物种进行清洗####
  speciesList <- unique(birdData$species)
  for (sp in speciesList) {
    # skip_to_next <- FALSE
    # tryCatch({
    spData <- subset(na.omit(birdData),species==sp)
    if(dim(spData)[1]>=108){
      print(paste(sp,'has been done!'))
      print(path)
      fwrite(spData,paste0('/root/autodl-tmp/GBIF_Data/data/',sp,'.csv'))
    }
  }
}


########去掉IUCN外的点，并thin raster#################
#IUCN boder
allPath2 <- list.files('/root/autodl-tmp/GBIF_Data/data',pattern='.csv',full.names=T)
for (path in allPath2) {
  spData <- fread(path)
  sp <- basename(path) %>% str_sub(.,1,-5)
  
  spBorderName <- grep(sp,borderList,value=T)
  if(length(spBorderName)>0){
    spBorder <- vect(spBorderName)
    spVect <- vect(spData,geom=c('lon','lat'),crs=crs)
    #筛选在IUCN边界内的观测点
    containPoints <- relate(spBorder,spVect,'contains',pairs=T)
    containData <- spData[containPoints[,2],]
    #观测记录数大于108
    # monthLength <- spData[,.N,keyby='month']
    #thin raster
    addGeom <- cellFromXY(globalRaster,containData[,c('lon','lat')]) %>% 
      cbind(id=.,containData)
    thinData <- unique(addGeom,by=1) %>% dplyr::select(.,-id)
    print(paste(sp,'has been done!'))
    fwrite(thinData,paste0('/root/autodl-tmp/GBIF_Data/filterData/',sp,'.csv'))
  }
}

#valis border
newValisBorder <- vect('/root/autodl-tmp/Wallace_zoogeographic/newValisBorder.shp')
valisName <- fread('/root/autodl-tmp/result/valisName.csv',header = T) %>% .[,1:2]
valisName$valisBorder <- ifelse(valisName$valisBorder=='',NA,valisName$valisBorder)
valisName2 <- na.omit(valisName)
# fwrite(valisName2,'/root/autodl-tmp/result/valisName2.csv')
for (i in 1:dim(valisName2)[1]) {
  valisPath <- grep(valisName2[i,1],allPath2,value = T,ignore.case=F)
  if(length(valisPath)>1){
    #有杂交情况
    valisPath <-valisPath[which.min(str_count(valisPath))]
  }
  if(length(valisPath)>0){
    valisDf <- fread(valisPath)
    spValis <- str_split(valisName2[i,2],',')
    valisSubBorder <- subset(newValisBorder,newValisBorder$name%in%spValis[[1]])
    spVect <- vect(valisDf,geom=c('lon','lat'),crs=crs)
    #筛选在华莱士界内的观测点
    containPoints <- relate(valisSubBorder,spVect,'contains',pairs=T)
    containData <- valisDf[containPoints[,2],]
    #thin raster
    addGeom <- cellFromXY(globalRaster,containData[,c('lon','lat')]) %>% 
      cbind(id=.,containData)
    thinData <- unique(addGeom,by=1) %>% dplyr::select(.,-id)
    print(paste(valisName2[i,1],'has been done!'))
    fwrite(thinData,paste0('/root/autodl-tmp/GBIF_Data/filterData/',valisName2[i,1],'.csv'))
  }
}

##########模型数据背景点采样################
climPath <- list.files("/root/autodl-tmp/全球基准时期BIO10m", pattern = ".tif", 
                       full.names = T)
crs <- '+proj=longlat +datum=WGS84'
globalRaster <- rast(vals=1:259200,nrows=360, ncols=720,xmin=-180, xmax=180,ymin=-90, ymax=90,crs=crs)
clim <- rast(climPath) %>% resample(.,globalRaster)
# 海拔数据
ele <- rast('/root/autodl-tmp/worldClim2dem数据10min/wc2.1_10m_elev.tif')%>% resample(.,globalRaster)
names(ele) <- 'ele'
#水体百分比
waterPer <-  rast('/root/autodl-tmp/waterPer/waterPer.tif')%>% resample(.,globalRaster)
names(waterPer) <- 'waterPer'

predictors <- rast()
`add<-`(predictors,clim) %>% `add<-`(.,ele)%>%  `add<-`(.,waterPer)

####相关性计算排除r>0.8#######
climateNames <- c('BIO_10','BIO_12','BIO_04','BIO_19','ele','waterPer')
allClimData <- predictors[[climateNames]]

#########按月切割##################
filterPath <- list.files('/root/autodl-tmp/GBIF_Data/filterData',full.names = T)
for(path in filterPath){
  df <- fread(path)
  name <- basename(path) %>% str_sub(.,1,-5)
  # polygon <- grep(name,polygonPath,value = T)
  # polyVect <- vect(polygon)
  for(m in unique(df$month)){
    dfMonth <- df[df$month==m,]
    ########背景值采样################
    background <- spatSample(x = allClimData,
                             size = dim(dfMonth)[1],    # generate 1,000 pseudo-absence points
                             values = FALSE, # don't need values
                             na.rm = TRUE,   # don't sample from ocean
                             as.df=T)      # just need coordinates
    background$class <- 0
    presence <- terra::extract(allClimData,dfMonth[,c('lon','lat')],ID=F)
    #####添加标签####
    presence$class <- 1
    presence_background <- rbind(presence,background)
    fwrite(presence_background,paste0('/root/autodl-tmp/GBIF_Data/monthData/',name,'_',m,'.csv'))
    print(paste0(name,'_',m,' has been done!'))
  }
}

writeRaster(ifel(is.na(allClimData[[1]]),-99999,allClimData[[1]]),'/root/autodl-tmp/clim_resampleData/bio_10.tif',overwrite=T)
writeRaster(ifel(is.na(allClimData[[2]]),-99999,allClimData[[2]]),'/root/autodl-tmp/clim_resampleData/bio_12.tif',overwrite=T)
writeRaster(ifel(is.na(allClimData[[3]]),-99999,allClimData[[3]]),'/root/autodl-tmp/clim_resampleData/bio_04.tif',overwrite=T)
writeRaster(ifel(is.na(allClimData[[4]]),-99999,allClimData[[4]]),'/root/autodl-tmp/clim_resampleData/bio_19.tif',overwrite=T)
writeRaster(ifel(is.na(allClimData[[5]]),-99999,allClimData[[5]]),'/root/autodl-tmp/clim_resampleData/ele.tif',overwrite=T)
writeRaster(ifel(is.na(allClimData[[6]]),-99999,allClimData[[6]]),'/root/autodl-tmp/clim_resampleData/waterPer.tif',overwrite=T)



















