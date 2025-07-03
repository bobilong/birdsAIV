library(biomod2)
library(ggplot2)
library(data.table)
library(ggplot2)
library(terra)
library(dplyr)
library(stringr)
library(paletteer)
library(rnaturalearth)
library(tidyr)
library(lubridate)

crs <- '+proj=longlat +datum=WGS84'
`%notin%` <- Negate(`%in%`)

world.map <- rnaturalearth::ne_countries(returnclass = "sf") |> dplyr::filter(continent != "Antarctica")
globalCountry <- vect(world.map) 
globalRaster <- rast(vals=1:259200,nrows=360, ncols=720,xmin=-180, xmax=180,ymin=-90, ymax=90,crs=crs)
coast <- rnaturalearth::ne_coastline(scale = "small", returnclass = "sf")


#####1.SDMs modeling##########
#using thin data and water percentage
climPath <- list.files('/root/climateResample/',full.names = T)
thinData <- list.files('/root/monthThinContainData/',full.names = T)
setwd('path/SDMs/')

for (path in thinData) {
  spMonth <- fread(path)
  spVect <- vect(spMonth[,c('lon','lat')],geom=c('lon','lat'),crs=crs)
  spName <- basename(path) %>% str_sub(.,1,-5)
  spVect$class <- '1'
  
  monthP<- str_split_i(basename(path),'_',2) %>% str_sub(1,-5)
  monthName <- paste0(monthP,'.tif.tif')
  
  preClimate <- grep(monthName,climPath,value=T) %>% rast()
  
  skip_to_next <- FALSE
  tryCatch({
    #1.Construct data, randomly pick 1000 times, repeat three times
    BiomodData <- BIOMOD_FormatingData(resp.var = spVect['class'],
                                       expl.var = preClimate,
                                       resp.name = spName,
                                       PA.strategy = 'random',
                                       PA.nb.rep=3,
                                       PA.nb.absences=1000,
                                       filter.raster = TRUE,
                                       na.rm=T)
    #2.build models
    myBiomodModelOut <- BIOMOD_Modeling(bm.format = BiomodData,
                                        # bm.options = bm.tuning$models.options,
                                        models = c('MAXNET','RF','XGBOOST'),
                                        CV.strategy = 'random',
                                        # CV.nb.rep = 0,
                                        # CV.k=5,
                                        CV.perc = 0.7,
                                        # var.import = 5,
                                        scale.models = T,
                                        prevalence = 0.5,
                                        # nb.cpu = 4,
                                        seed.val = 123,
                                        metric.eval = c('TSS'))
    #3.Integration Models
    myBiomodEM<- BIOMOD_EnsembleModeling(bm.mod = myBiomodModelOut,
                                         models.chosen = 'all',
                                         em.by = 'all',
                                         em.algo = 'EMmean',
                                         metric.select = c('TSS'),
                                         metric.select.thresh = c(0.8),
                                         metric.eval = c('TSS'),
                                         # var.import = 5,
                                         # EMci.alpha = 0.05,
                                         # nb.cpu = 4
    )
    #Single Model outputs
    myBiomodProj <- BIOMOD_Projection(bm.mod = myBiomodModelOut,
                                      proj.name = spName,
                                      new.env = preClimate,
                                      models.chosen = 'all',
                                      metric.binary = 'TSS',
                                      # metric.filter = 'all',
                                      build.clamping.mask = TRUE)
    
    
    #Integration of model outputs
    myBiomodEMProj <- BIOMOD_EnsembleForecasting(bm.em = myBiomodEM, 
                                                 # bm.proj = myBiomodProj,
                                                 new.env = preClimate,
                                                 proj.name=paste0(spName,'_ensem'),
                                                 models.chosen = 'all',
                                                 metric.binary = 'TSS',
                                                 # metric.filter = 'TSS',
                                                 output.format='.tif')
    
    print(paste0(spName,' has been done!'))
  },error = function(e) { 
    print(e)
    skip_to_next <<- TRUE})
}


#####2.SDMs validation##########
#TSS
model_dir <- "/root/autodl-tmp/WAEdata_new_y/SDMs"
# output_dir <- "/root/result"
# dir.create(output_dir, showWarnings = FALSE)
list.files(model_dir)

model_files <- list.files(
  model_dir,
  pattern = "\\.models\\.out$",
  full.names = TRUE,
  recursive = TRUE
)
# remove ensemble.models.out
model_files <- model_files[!grepl("ensemble", model_files)]

calibration_list <- list()
for (file in model_files) {
  model_name <- gsub(".models.out", "", basename(file)) #.ensemble.models.out
  before <- ls()
  load(file)
  new_obj <- setdiff(ls(), before)
  em_obj <- new_obj[sapply(new_obj, function(x) inherits(get(x), "BIOMOD.models.out"))][1]
  eval_data <- tryCatch(get_evaluations(get(em_obj)), error = function(e) NULL)
  # read calibration
  calibration_list[[model_name]] <- data.frame(
    Model = model_name,
    Calibration = eval_data$calibration
  )
  # remove
  rm(list = new_obj)
  gc()
  cat(model_name,'has done\n')
}

all_calib <- do.call(rbind, calibration_list)
fwrite(all_calib, '/root/result/SDMs_TSS.csv')

#plot
ggplot(all_calib, aes(x = "", y = Calibration)) +
  geom_boxplot(fill = "#66c2a5") +
  theme_minimal()
mean_calib <- mean(all_calib$Calibration, na.rm = TRUE)
cat("Mean Calibration：", round(mean_calib, 3), "\n")



#####3.Statistic for Main Figures##########
######(1) Waterbird activity entropy (WAE) calculations###########
speciesPixelNumPath2 <- list.files('/root/singleRast_779',pattern = '.tif',full.names = T)

#Species Richness calculation
spNum <- lapply(speciesPixelNumPath2, function(x){
  r <- rast(x) %>% sum(.,na.rm=T)
  r2 <- ifel(r>0,1,NA)
  names(r2) <- str_sub(basename(x),1,-5)
  return(r2)
})
spNum <- rast(spNum)
spNumTif <- sum(spNum,na.rm=T) %>% mask(globalCountry)
terra::writeRaster(spNumTif,'/root/result/spNumTif.tif')

#Species Richness CV calculation
spMonth <- rast(speciesPixelNumPath2)
sp_perMonth <- tapp(spMonth,names(spMonth),sum,na.rm=T) %>% mask(.,globalCountry)
sp_perMonth <- sp_perMonth[[c('X1','X2','X3','X4','X5','X6','X7','X8','X9','X10','X11','X12')]]

sp_sd<-app(sp_perMonth, sd, na.rm = T)
sp_mean<-app(sp_perMonth, mean, na.rm = T)
sp_meanthan1<-ifel(sp_mean>0.2, sp_mean, NA)
sp_cv <- sp_sd/sp_meanthan1*100
plot(sp_cv)
writeRaster(sp_cv,paste0('/root/result/sp_cv.tif'))

#WAE calculation
speciesPixelNum2 <- rast(speciesPixelNumPath2)
allMonth <- sum(speciesPixelNum2,na.rm=T)%>% mask(globalCountry)

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
writeRaster(AE, paste0('/root/result/WAE.tif'), overwrite = TRUE)


######—— Fig.1###############
#Fig.1 a
#
spNumTif<-rast(paste0('/root/autodl-tmp/WAEdata_new_y/result779/spNumTif.tif'))
pspNum <- ggplot() +
  geom_spatraster(data = spNumTif) +
  geom_spatvector(data=coast,fill=NA)+coord_sf(crs = crs,xlim=c(-160,165),ylim=c(-56,90))+
  theme_void()+
  scale_fill_gradientn(colours = paletteer_c("grDevices::Zissou 1", 30) ,na.value='white')+
  labs(fill="Species Richness")+
  theme(
    plot.title = element_text(hjust=0.5),
    legend.title = element_blank(),
    #legend.title = element_text(hjust=0.5),
    legend.title.align = -10,
    legend.position = 'bottom',#设置图例与主图距离
    legend.direction='horizontal',#图例水平放置vertical，垂直horizontal，
    legend.key.width = unit(1.15,'cm'), #图例宽度
    legend.key.height = unit(0.3,'cm')
  )
pspNum
#
sp_cv <- rast('/root/result/sp_cv.tif')
pSRCV<-ggplot() +
  geom_spatraster(data = sp_cv) +  #sp_cv
  geom_spatvector(data=coast,fill=NA)+coord_sf(crs = crs,xlim=c(-160,165),ylim=c(-56,90))+
  theme_bw()+
  scale_fill_gradientn(colours = paletteer_c("grDevices::Zissou 1", 30) ,na.value='white',values = c(0, 0.3, 1))+
  #labs(title = "CV")+
  theme(
    axis.text = element_text(size=12),
    plot.title = element_text(hjust=0.5),
    # axis.line = element_line(color = 'black'),
    # panel.background = element_rect(fill = "white"),#背景设置
    #legend.direction = "horizontal",
    legend.title = element_blank(),
    #legend.title = element_text(hjust=0.5),
    legend.title.align = -10,
    legend.position = c(0.15, 0.05),#设置图例与主图距离
    legend.direction='horizontal',#图例水平放置vertical，垂直horizontal，
    legend.key.width = unit(1.2,'cm'), #图例宽度
    legend.key.height = unit(0.3,'cm')
    
  )
pSRCV
#
sp_df <- sp_perMonth[[c('X1','X2','X3','X4','X5','X6','X7','X8','X9','X10','X11','X12')]]
names(sp_df) <- c('Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec')
sp_df <- terra::as.data.frame(sp_df, xy = TRUE, na.rm = FALSE) 
head(sp_df)
sp_df_long <- gather(sp_df, key = "Month", value = "Value", -x, -y,na.rm = T)
names(sp_df_long)<-c("lon","lat","Month","Value")
head(sp_df_long)
#sp_df_long$Month <- factor(sp_df_long$Month, levels = c('X1','X2','X3','X4','X5','X6','X7','X8','X9','X10','X11','X12'))
sp_df_long$Month <- factor(sp_df_long$Month, levels = c("Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"))
sp_df_avg <- sp_df_long %>%
  group_by(lat, Month) %>%
  summarise(Avg_Value = mean(Value))
meanValue<-mean(sp_df_long$Value)

pSRCV2<-ggplot(sp_df_avg, aes(x=Avg_Value, y=lat, color = Month, group = Month)) +
  geom_path(linewidth=0.5)+
  xlab('Species Richness')+
  xlim(c(0,90))+
  ylim(c(-56,90))+
  geom_vline(xintercept = meanValue,color='grey',linetype='dashed')+
  theme_bw()+
  coord_fixed(ratio = 1.7)+ #
  # theme_classic()+
  theme(
    legend.position = "none",
    axis.title.y = element_blank(),
    axis.text.y=element_blank(),
    axis.title.x = element_text(size=16),
    axis.text.x=element_text(size=12),
    panel.grid = element_blank(),
    axis.ticks.length = unit(-4, "pt") #
  )

pSRCV+pSRCV2



#Fig.1 b
#
Vborder <- vect('/root/autodl-tmp/Wallace_zoogeographic/newValisBorder.shp') 
vRaster <- rasterize(Vborder,AE,field='name')
valisEntropy <- c(AE,vRaster) %>% terra::as.data.frame() %>%na.omit()
valisEntropy <- valisEntropy[valisEntropy$name != "Antarctica",]   #
Vborder1<-factor(valisEntropy$name,levels = c("Oriental","Ethiopian","Australian","Neotropical","Nearctic","Palearctic")) #

ggplot(valisEntropy, aes(x = Vborder1, y = sum,fill=Vborder1)) + #,
  scale_fill_manual(values = c("#45F3DD","#F6BA26","#E0312D","#3ca1f9","#8204FB","#49FB18"))+
  geom_hline(yintercept = mean(valisEntropy$sum), size = 1, color = "grey", linetype = "dashed")+
  geom_violin(alpha=0.5,position = position_dodge(width = 0.01), scale = 'width') +
  geom_boxplot(alpha=0.8,width=0.45,position=position_dodge(width=0.1),size=0.75,outlier.colour = NA)+
  ylim(c(0,5.5))+
  theme(panel.grid = element_blank(), 
        panel.background = element_rect(fill = 'transparent', color = "black"), 
        legend.position = 'none',
        axis.text = element_text(size = 30),
        text = element_text(size = 34),) +
  labs(y = "WAE", x = " ")
#
AE<-rast(paste0('/root/result/WAE.tif')); plot(AE)
library(paletteer)
pAE <- ggplot() +
  geom_spatraster(data = AE) +
  geom_spatvector(data=coast,fill=NA)+coord_sf(crs = crs,xlim=c(-160,165),ylim=c(-56,90))+
  theme_bw()+
  scale_fill_gradientn(colours = paletteer_c("grDevices::Zissou 1", 30) ,na.value='white',values = c(0, 0.75, 1))+
  #scale_fill_brewer(palette = "RdYlGn")+
  labs(fill="Species Richness")+
  theme(
    axis.text = element_text(size=12),
    plot.title = element_text(hjust=0.5),
    legend.title = element_blank(),
    #legend.title = element_text(hjust=0.5),
    legend.title.align = -10,
    legend.position = c(0.15, 0.05),#
    legend.direction='horizontal',#，
    legend.key.width = unit(1.15,'cm'), #
    legend.key.height = unit(0.3,'cm')
  )
pAE

library(ggpointdensity)
AE_df <- terra::as.data.frame(AE, xy = TRUE, na.rm = T) 
head(AE_df)
names(AE_df)<-c("lon","lat","sum")

pAE2<-ggplot(AE_df,aes(x=sum, y=lat),)+
  geom_pointdensity(adjust = 4, size = 0.5)+                              # 
  scale_color_distiller(palette = "Spectral", direction = -1)+ # 
  xlab('WAE')+
  xlim(c(0,5.5))+
  ylim(c(-56,90))+
  theme_bw()+
  coord_fixed(ratio = 0.105)+ #
  # theme_classic()+
  theme(
    axis.title.x = element_text(size=16),
    axis.text.x=element_text(size=12),
    axis.title.y = element_blank(),
    axis.text.y=element_blank(),
    panel.grid = element_blank(),
    axis.ticks.length = unit(-4, "pt"), #
    legend.position = "none"
  )
pAE2

pAE+pAE2

######(2) Performance of WAE for predicting reported cases of AIV ###########

# Multi-criteria Filtering Process for Restricting Poultry-to-Poultry Transmission Cases (Supplementary Fig. 7)
### result 1 (raw data)
#### all accessible global occurrences of AIVs reports (including those detected in birds, mammals, and humans) from FAO during the period 2000-2022.
outBreakData <- fread('/root/autodl-tmp/result/result779/outBreakData.csv')
outBreakData$label <- str_extract(outBreakData$Serotype,'HPAI|LPAI')
outBreakData$h_label <- str_extract(outBreakData$Serotype,'H[0-9]N[0-9]|H[0-9]')

### Step.1：Subtype-Based Filtering —— AUCresult_1.1.csv
outBreakData_num <- fread('/root/autodl-tmp/result2/result753/outBreakData_num.csv')
outBreakData_num$Serotype <- paste(outBreakData_num$genoType,outBreakData_num$High_low)
con <- outBreakData_num[poultryPropertion>=80,]
df1 <- outBreakData[Serotype%notin%con$Serotype,]
fwrite(df1, '/root/result/AUC_data/AUCresult_1.1.csv')

### Step.2：Spatial -Temporal Sparsification  —— AUCresult_1.2.csv
df1 <- fread('/root/autodl-tmp/WAEdata_new_y/result/nc_review/auc_data/AUCresult_1.1.csv')
df1poul <- subset(df1, dieasSource == 'Poultry')
df1other <- subset(df1, dieasSource != 'Poultry')
# time
df1poul <- df1poul %>% mutate(
    obs_date = as.Date(Observation.date..dd.mm.yyyy., format = "%d/%m/%Y"),
    year = year(obs_date),
    month = month(obs_date)  )
# NA 
df1poul <- df1poul[!is.na(df1poul$year), ]
# 10km
globalRaster10 <- rast(  vals = 1:6480000,  nrows = 1800, ncols = 3600,  xmin = -180, xmax = 180,  ymin = -90, ymax = 90,  crs = "+proj=longlat +datum=WGS84")
df1poul <- df1poul %>%  mutate(rasterID = cellFromXY(globalRaster10, df1poul[, c('Longitude', 'Latitude')]))

final_result <- data.table()
years <- min(df1poul$year):(max(df1poul$year) - 1)
setDT(df1poul)  # 
final_result <- data.table()
for (y in years) {
  # 
  cold_season <- df1poul[
    (year == y & month %in% c(11, 12)) | (year == y + 1 & month %in% 1:3)
  ][order(obs_date)][, season := paste0(y, "_cold")][, .SD[1], by = rasterID]  #
  # 
  warm_season <- df1poul[
    year == y & month %in% 5:9
  ][order(obs_date)][, season := paste0(y, "_warm")][, .SD[1], by = rasterID]  #
  # 
  final_result <- rbind(final_result, cold_season, warm_season)
}

final_result

final_result2 <- final_result[, !c("obs_date", "year", "month", "season"), with = FALSE]
merged_result <- merge(final_result2, df1other, by = intersect(names(final_result2), names(df1other)), all = TRUE)

fwrite(merged_result, '/root/result/AUC_data/AUCresult_1.2.csv')


### Step.3：Spatial Filtering   —— AUCresult_1.3.csv
# crs <- '+proj=longlat +datum=WGS84'
# globalRaster <- rast(vals=1:259200,nrows=360, ncols=720,xmin=-180, xmax=180,ymin=-90, ymax=90,crs=crs)
# 
# df2 <- fread('/root/result/AUC_data/AUCresult_1.2.csv')
# 
# unique(df2$dieasSource)
# otherDf <- df2[dieasSource != 'Poultry']
# WildDf <- df2[dieasSource == 'Wild_bird']
# PoultryDf <- df2[dieasSource == 'Poultry']
# 
# PoultryDf$id <- cellFromXY(globalRaster, PoultryDf[, c('Longitude', 'Latitude')])
# otherDf$id <- cellFromXY(globalRaster, otherDf[, c('Longitude', 'Latitude')])
# 
# otherID <- cellFromXY(globalRaster, otherDf[, c('Longitude', 'Latitude')])
# wildBirdID <- cellFromXY(globalRaster, WildDf[, c('Longitude', 'Latitude')])
# 
# df3poultry <- PoultryDf[id %in% wildBirdID, ]
# df3<- rbind(otherDf, df3poultry)
# 
# length((subset(df3, df3$dieasSource=='Poultry'))$dieasSource)
# unique(df3$Serotype)
# fwrite(df3, '/root/result/AUC_data/AUCresult_1.3.csv')
# 




#ROC calculation
library(pROC)
`%notin%` <- Negate(`%in%`)
crs <- '+proj=longlat +datum=WGS84'
globalRaster <- rast(vals=1:259200,nrows=360, ncols=720,xmin=-180, xmax=180,ymin=-90, ymax=90,crs=crs)
globalSHP <- vect('/root/autodl-tmp/allData/worldBorder/continentNew.shp')
AE <- rast('/root/autodl-tmp/WAEdata_new_y/result779/WAE.tif')

outBreakData <- fread('/root/result/AUC_data/AUCresult_1.2.csv')
#ROC
outBreakData2 <- na.omit(outBreakData[,c('Longitude','Latitude')])
TP_points <- vect(outBreakData2,geom=c('Longitude','Latitude'),crs=crs)
breakRast <- rasterize(TP_points,globalRaster)

sample_TP <- terra::extract(AE,TP_points,min,bind=T) %>% as.data.frame() %>% na.omit()
sample_TP$label <- 'good'

trueNgRast <- mask(AE,breakRast,inverse=T)

set.seed(1234)
sample_TN<- terra::spatSample(trueNgRast,dim(sample_TP)[1],method='random',exhaustive=T,na.rm=T)
sample_TN$label <- 'poor'
sampleData <- rbind(sample_TN,sample_TP[,c('sum','label')])
dfroc1<- roc(sampleData$label, sampleData$sum)



#Accuracy calculation for countries
outBreakData$label <- str_extract(outBreakData$Serotype,'HPAI|LPAI')
outBreakData$h_label <- str_extract(outBreakData$Serotype,'H[0-9]N[0-9]|H[0-9]')

outBreakData3 <- na.omit(outBreakData2[,c('Country','Longitude','Latitude')])
addGeom <- cellFromXY(globalRaster,outBreakData3[,c('Longitude','Latitude')]) %>% 
  cbind(id=.,outBreakData3)
thinData <- unique(data.table(addGeom),by='id') %>% dplyr::select(.,-id)
breakVect <- vect(thinData,geom=c('Longitude','Latitude'),crs=crs)
breakRast <- rasterize(breakVect,globalRaster)
#
AE <- rast('/root/autodl-tmp/result2/result753/WAE.tif')
TP_vect <- vect(outBreakData3,geom=c('Longitude','Latitude'),crs=crs)
sample_TP <- terra::extract(AE,TP_vect,mean,bind=T) %>% as.data.frame() %>% subset(sum>2)
sample_TP$label <- 'good'

trueNgRast <- mask(AE,breakRast,inverse=T)
set.seed(2345)
sample_TN<- terra::spatSample(trueNgRast,dim(sample_TP)[1],method='random',exhaustive=T,na.rm=T,xy=T)
sample_TN$label <- 'poor'

countrRast <- rasterize(globalSHP,AE,'NAME_LONG')
sample_TN_vect <- vect(sample_TN,geom=c('x','y'),crs=crs)
Country <- terra::extract(countrRast,sample_TN_vect)
sample_TN$Country <- Country$NAME_LONG

threshold <- 4.108
result_country <- data.table()
for (countryname in unique(sample_TP$Country)) {
  country <- subset(globalSHP,globalSHP$NAME_LONG==countryname)
  if(length(country)>0){
    df <- subset(thinData,thinData$Country==countryname)
    df_T <- subset(sample_TP,sample_TP$Country==countryname)
    df_N <- subset(sample_TN,sample_TN$Country==countryname)
    countryRast <- crop(AE,country,mask=T)
    allDf <- as.data.frame(countryRast)
    TP <- sum(df_T$sum>=threshold,na.rm=T)
    FN <- dim(df_T)[1]-TP
    TN <- sum(df_N$sum<threshold,na.rm=T)
    FP <- dim(df_N)[1]-TN
    
    result_country <- rbind(result_country, data.frame(country = countryname,df=dim(df)[1],alldf=nrow(allDf),TP=TP,FN=FN,FP=FP,TN=TN))
  }
}
result_country$accuracy <- (result_country$TP+result_country$TN)/(result_country$TP+result_country$TN+result_country$FN+result_country$FP)
result_country$precise <- (result_country$TP)/(result_country$TP+result_country$FP)
result_country$recall <- (result_country$TP)/(result_country$TP+result_country$FN)
result_country$FI <- 2*((result_country$recall*result_country$precise)/(result_country$recall+result_country$precise))

fwrite(result_country,'/root/autodl-tmp/result2/result753/result_country.csv')



#Sensitivity calculation for AIV types
AE <- rast('/root/autodl-tmp/WAEdata_new_y/result779/WAE.tif')
outBreakData <- fread('/root/result/AUC_data/AUCresult_1.2.csv')

outBreakData$label_hl <- str_extract(outBreakData$Serotype,'HPAI|LPAI')
outBreakData$h_label <- str_extract(outBreakData$Serotype,'H[0-9]N[0-9]|H[0-9]')

h5n1_df <- subset(outBreakData,outBreakData$Serotype%in%c('H5N1 HPAI','H5N1 LPAI'))
hx_df<- outBreakData
h5n1_vect <- vect(h5n1_df,geom=c('Longitude','Latitude'),crs=crs)
hx_vect <- vect(hx_df,geom=c('Longitude','Latitude'),crs=crs)

h5n1_value <- terra::extract(AE,h5n1_vect,mean,na.rm=T) %>% 
  as.data.frame() %>% na.omit()
hx_value <- terra::extract(AE,hx_vect,mean,na.rm=T) %>% 
  as.data.frame() %>% na.omit()

H5N1_df <- data.table(label='H5N1', cases=nrow(h5n1_df),recall=nrow(subset(h5n1_value,h5n1_value$sum>4.108))/nrow(h5n1_value))
Hx_df <- data.table(label='all AIV',cases=nrow(hx_df),recall=nrow(subset(hx_value,hx_value$sum>4.108))/nrow(hx_value))

result_st <- rbind(H5N1_df,Hx_df,fill=T)
result_st

h_df <- subset(outBreakData,outBreakData$label_hl%in%c('HPAI'))
l_df <- subset(outBreakData,outBreakData$label_hl%in%c('LPAI'))

h_vect <- vect(h_df,geom=c('Longitude','Latitude'),crs=crs)
l_vect <- vect(l_df,geom=c('Longitude','Latitude'),crs=crs)

h_value <- terra::extract(AE,h_vect,mean,na.rm=T) %>%   as.data.frame() %>% na.omit()
l_value <- terra::extract(AE,l_vect,mean,na.rm=T) %>%   as.data.frame() %>% na.omit()

HPAI_df <- data.table(label='HPAI', cases=nrow(h_df), recall=nrow(subset(h_value,h_value$sum>4.108))/nrow(h_value))

LPAI_df <- data.table(label='LPAI', cases=nrow(l_df), recall=nrow(subset(l_value,l_value$sum>4.108))/nrow(l_value))

result_st4 <- rbind(result_st,HPAI_df,fill=T) %>% rbind(.,LPAI_df,fill=T)
result_st4
fwrite(result_st4, '/root/result/Serotype_edit.csv')



######—— Fig.2###############
#Fig.2 a
#litte
plot(dfroc1,col="red",
     legacy.axes=T,
     print.auc=T,
     expand=c(0,0),
     print.thres=T,
     grid=c(0.2,0.2),grid.col=c("lightgrey","lightgrey"),
     cex.main=1.5,
     cex.sub=1.5,
     cex.axis=1.5,
     cex.lab=1.5)


plotEntropy <-  ifel(WAE>=4.108,1,0)
ggplot() +
  geom_spatraster(data = plotEntropy) +
  scale_fill_gradient(low = "lightgrey",high = "yellow" ,space = "Lab",n.break=2,
                      labels=c('Low risk','High risk'),na.value='white')+
  geom_spatvector(data=coast,fill=NA)+
  #coord_sf(crs = crs,xlim=c(-160,180),ylim=c(-56,90))+
  geom_point(data=outBreakData,aes(Longitude,Latitude, color=newLabel),size=0.8,alpha = 0.8,stroke = 0.3)+   #shape = 21是圆圈，shape = 4是×
  scale_color_manual(values = c('Poultry'='red','Wild_bird'='blue', 'Other'='purple'))+
  scale_shape_binned(c(4,21))+
  theme_bw()+
  #scale_color_manual(values='#fc4e4e', labels='Outbreak of navian influenza') +
  guides(
    fill=guide_legend()  
  )+
  labs(x='',y='')+
  theme(
    text = element_text(size = 15),
    #axis.text = element_text(size = 15),
    plot.title = element_text(hjust=0.5),
    legend.title = element_blank(),
    legend.text = element_text(size = 15),
    legend.position = "none",
    #legend.position = c(0.15, 0.15),#
    legend.direction='horizontal',#
    legend.key.width = unit(0.4,'cm'), #
    legend.key.height = unit(3,'cm')
  )


plotEntropy <-  ifel(WAE>=4.108,1,0)
ones <- sum(values(plotEntropy) == 1, na.rm = TRUE); ones
zeros <- sum(values(plotEntropy) == 0, na.rm = TRUE); zeros

ones/(ones+zeros)
1-ones/(ones+zeros)


#Fig.2 b + c
result_country <- fread('/root/result/result_country.csv')
globalSHP <- vect('/root/autodl-tmp/WAEdata_new_y/result779/otherdata/ne_10m_admin_0_countries.shp')
countryData <- as.data.frame(globalSHP) %>% .[,c('NAME_LONG','POP_EST','POP_RANK','ECONOMY','INCOME_GRP','CONTINENT','GDP_MD')]
result_df4 <- left_join(result_country,countryData,by=c('country'='NAME_LONG'))
result_df4$per <- (result_df4$df/result_df4$alldf)

#GDP
pop_quantiles <- quantile(result_df4$POP_EST, probs=seq(0, 1, by=0.1), na.rm = TRUE)
result_df4$POP_EST_size <- findInterval(result_df4$POP_EST, vec = pop_quantiles)
result_df4$POP_EST <- as.numeric(result_df4$POP_EST)
result_df4 <- mutate(result_df4, POP_EST_label = cut(POP_EST, breaks = c(0, 10000000, 20000000, 50000000, 100000000, 500000000, 1000000000, 1397800000),
                                                     labels = c("<1000", "1000 - 2000", "2000-5000", "5000 - 10000", "10000 - 50000", "50000 - 100000", ">100000")))
result_df4 <- result_df4[order(result_df4$CONTINENT, -result_df4$POP_EST), ]
country5 <- c("China", "United States", "India", "Nigeria")
result_dfc <- result_df4[result_df4$country %in% country5]
result_dfc <- result_dfc %>%
  mutate(country = ifelse(country == "United States", "USA", country))
size_values <- c("<1000" = 3, "1000 - 2000" = 4, "2000-5000" = 6, "5000 - 10000" = 8, "10000 - 50000" = 10, "50000 - 100000" = 12, ">100000" = 18)

my_formula <- y ~ x
p1<-ggplot(result_df4, aes(log(GDP_MD), accuracy, color=CONTINENT)) +
  geom_point(aes(size=POP_EST_label), alpha=0.5,show.legend = F) +
  geom_text(data=result_dfc, aes(label=country, x=log(GDP_MD), y=accuracy),show.legend = F,size=4.5, hjust=0.5, vjust=2) +
  scale_color_manual(values = c("#984EA3", "#E41A1C","#4DAF4A",  "#FF7F00", "#a38900","#377EB8"))+
  #geom_smooth(method = "gam", se = T, fill="lightgrey",color = "grey", alpha=0.3, linetype = "dashed") +  # 
  geom_smooth(aes(log(GDP_MD),accuracy),show.legend = F,fill="lightgrey",color = "grey",method = "lm", formula = y ~ I(x^-1))+
  xlab('GDP (Log)') +
  ylab('Accuracy') +
  ylim(c(0,1))+
  xlim(c(7.5,17))+
  theme_bw() +
  theme(panel.background = element_blank(),
        text = element_text(size=18)) +
  scale_size_manual(values = size_values) +
  guides(size=guide_legend(title=expression(paste("Human population (" ~ 10^4 ~")"))),
         color=guide_legend(title="Continent"))

p1

#Prevalence
result_country <- fread('/root/autodl-tmp/humPoulResult/data/result_country.csv')
globalSHP <- vect('/root/autodl-tmp/worldBorder/ne_10m_admin_0_countries.shp')
countryData <- as.data.frame(globalSHP) %>% .[,c('NAME_LONG','POP_EST','POP_RANK','ECONOMY','INCOME_GRP','CONTINENT','GDP_MD')]
result_df3 <- left_join(result_country,countryData,by=c('country'='NAME_LONG'))
result_df4$per <- (result_df4$df/result_df4$alldf)

# POP_EST
pop_quantiles <- quantile(result_df4$POP_EST, probs=seq(0, 1, by=0.1), na.rm = TRUE)
result_df4$POP_EST_size <- findInterval(result_df4$POP_EST, vec = pop_quantiles)
result_df4$POP_EST <- as.numeric(result_df4$POP_EST)
result_df4 <- mutate(result_df4, POP_EST_label = cut(POP_EST, breaks = c(0, 10000000, 20000000, 50000000, 100000000, 500000000, 1000000000, 1397800000),
                                                     labels = c("<1000", "1000 - 2000", "2000-5000", "5000 - 10000", "10000 - 50000", "50000 - 100000", ">100000")))
result_df4 <- result_df4[order(result_df4$CONTINENT, -result_df4$POP_EST), ]
country5 <- c("China", "United States", "India", "Nigeria")
result_dfc <- result_df4[result_df4$country %in% country5]
result_dfc <- result_dfc %>%
  mutate(country = ifelse(country == "United States", "USA", country))

size_values <- c("<1000" = 3, "1000 - 2000" = 4, "2000-5000" = 6, "5000 - 10000" = 8, "10000 - 50000" = 10, "50000 - 100000" = 12, ">100000" = 18)
my_formula <- y ~ x
library(scales)
my_trans <- trans_new(
  name = "custom",
  transform = function(x) x + 0.1 * (x - 1),
  inverse = function(x) (x - 1) / 1.1 + 1
)
result_df5<-result_df4[result_df4$country!="Cyprus",]
p2<-ggplot(result_df5, aes(per, accuracy, color=CONTINENT)) +
  geom_point(aes(size=POP_EST_label), alpha=0.5) +
  ylim(c(0,1))+
  geom_text(data=result_dfc, aes(label=country, x=per, y=accuracy),size=4.5, hjust=0.5, vjust=2) +
  scale_color_manual(values = c("#984EA3", "#E41A1C","#4DAF4A",  "#FF7F00", "#a38900","#377EB8"))+
  #geom_smooth(method = "gam", se = T, fill="lightgrey",color = "grey", alpha=0.3, linetype = "dashed") +  # 
  geom_smooth(aes(per,accuracy),fill="lightgrey",color = "grey",method = "loess", formula = y ~ I(x^-1))+
  xlab('Prevalence') +
  ylab('Accuracy') +
  theme_bw() +
  theme(panel.background = element_blank(),
        text = element_text(size=18),
        legend.position = 'right') +
  scale_size_manual(values = size_values) +
  guides(size=guide_legend(title=expression(paste("Human population (" ~ 10^4 ~")"))),
         color=guide_legend(title="Continent"))

p1+p2


#Fig.2 d
result_st2 <- fread('/root/result/Serotype_edit.csv')
ggplot(data = result_st2, aes(x = label, y = recall, fill = label)) +
  geom_col(position = 'dodge2') +
  theme_bw() +  # 使用黑白主题作为基础
  labs( y = "Sensitivity") +  
  scale_fill_brewer(palette = "Pastel1") +
  ylim(c(0, 0.9))+
  theme(
    axis.title.x = element_blank(),
    text = element_text(size=18),
    plot.title = element_text(hjust = 0.5),  #
    legend.position = "none"  
  )
result_st3 <- mutate(result_st2, df_label = cut(cases, breaks = c(1500, 2000,3000, 10000, 20000,30000),
                                                labels = c("1500 - 2000", "2000 - 3000","5000 - 6000", "10000 - 20000", "> 25000")))
size_values <- c("1500 - 2000" = 4, "2000 - 3000" = 6,"5000 - 6000" = 8, "10000 - 20000" = 14, "> 25000" =20)
result_st3$Group<-c("Subtype","Subtype","Pathogenicity","Pathogenicity")
result_st3$Group<-factor(result_st3$Group,levels = c("Subtype","Pathogenicity"))
result_st3$label<-factor(result_st3$label,levels = c("H5N1","all AIV","HPAI","LPAI"))

ggplot(data = result_st3, aes(x = label, y = recall)) +
  geom_point(aes(size = df_label, color = Group), alpha = 0.5) +  # 
  geom_hline(aes(yintercept = 0.804), color = "grey", linetype = "dashed", alpha = 1) +
  scale_size_manual(values = size_values) +  # 
  scale_color_manual(values = c("orange", "red")) +  # 
  ylim(c(0.8, 0.9)) +
  labs(size = "Number of Cases", y = "Sensitivity") +  # 
  theme_bw() +
  theme(
    axis.text = element_text(size = 12),
    axis.title.y = element_text(size = 16),
    axis.title.x = element_blank()
  )


######(3) Hotspot identification ###########

#hotspot calculation
crs <- '+proj=longlat +datum=WGS84'
globalRaster <- rast(vals=1:259200,nrows=360, ncols=720,xmin=-180, xmax=180,ymin=-90, ymax=90,crs=crs)

popd2015<-rast("/root/autodl-tmp/WAEdata_new_y/result779/otherdata/pop.tif") %>% resample(globalRaster)
Entropy<-rast("/root/autodl-tmp/WAEdata_new_y/result779/WAE.tif")
poul2015<-rast("/root/autodl-tmp/WAEdata_new_y/result779/otherdata/poultry.tif")%>% resample(globalRaster)
cattle2015<-rast("/root/autodl-tmp/WAEdata_new_y/result779/otherdata/cattle.tif") %>% resample(globalRaster)

quantiles_pop <- global(popd2015,quantile,probs=seq(0, 1, 0.01),na.rm=T)
quantiles_entr<- global(Entropy,quantile,probs=seq(0, 1, 0.01),na.rm=T)
quantiles_poul<- global(poul2015,quantile,probs=seq(0, 1, 0.01),na.rm=T)
quantiles_cattle<- global(cattle2015,quantile,probs=seq(0, 1, 0.01),na.rm=T)

quan_pop <- 60  #
quan_entr<- 4.108  #
quan_poul <- 11000  #
quan_cattle <- 700 #

#reclass
pop_new <- ifel(popd2015>=quan_pop,1,0)
entr_new <- ifel(Entropy>=quan_entr,10,0)
poul_new <- ifel(poul2015>=quan_poul,2,0)
cattle_new <- ifel(cattle2015>=quan_cattle,4,0)

hotentrpoppoulcat <- cover(pop_new + entr_new + poul_new + cattle_new, entr_new)

names(hotentrpoppoulcat) <-'type'
hotAND <- ifel(hotentrpoppoulcat==17,1,0)  ; plot(hotAND); global(hotAND,sum,na.rm=T)
writeRaster(hotentrpoppoulcat, "/root/result/Hot_data/hotentrpoppoulcat.tif")


#hot_statistic
plot(hotentrpoppoulcat)
hotentrpoppoulcat_df <- hotentrpoppoulcat%>% terra::as.data.frame(.,xy=T) %>%na.omit() ; head(hotentrpoppoulcat_df)
hot4_df<-c(hotentrpoppoulcat, popd2015, poul2015 ,cattle2015) %>% terra::as.data.frame(.,xy=T) %>%na.omit()
head(hot4_df)
names(hot4_df)<-c("x","y","type","pop","poul","cattle")

type_counts <- hotentrpoppoulcat_df %>% 
  group_by(type) %>% 
  summarise(count = n())

Num_result <- hot4_df %>%
  group_by(type) %>%
  summarise(across(everything(), ~ sum(.x, na.rm = TRUE)))

Num_result <- left_join(Num_result, type_counts, by = "type")
Num_result

custom_order <- c(12, 14, 11, 16, 15, 13, 17, 10, 0,1,2,3,4,5,6,7) #  
Num_result$type <- factor(Num_result$type, levels = custom_order) #  
Num_result2 <- Num_result %>% arrange(type)
Num_result2

fwrite(Num_result2, '/root/Hot_data/hot_statistic.csv')




#hot_statistic for countries
pop_new <- ifel(popd2015>=quan_pop,1,0)
entr_new <- ifel(Entropy>=quan_entr,10,0)
poul_new <- ifel(poul2015>=quan_poul,2,0)
cattle_new <- ifel(cattle2015>=quan_cattle,4,0)
plot(pop_new+poul_new+cattle_new)

hotentrpoppoulcat <- cover(pop_new + entr_new + poul_new + cattle_new, entr_new)   ;plot(hotentrpoppoulcat)

hotentrpoppoulcat_1 <- ifel(cattle2015>quan_cattle, 1, NA)
countrys <- vect ('/root/autodl-tmp/WAEdata_new_y/result779/otherdata/Con_EU_dis.shp')


regionRast <- rasterize(countrys, hotentrpoppoulcat, field="EU753")
hotentrpoppoulcat_albers <- hotentrpoppoulcat

# 统计不同值的区域数据
hot_11 <- ifel(hotentrpoppoulcat_albers == 11, 1, 0)
names(hot_11) <- 'hot_11'
resultregion_11 <- zonal(hot_11, regionRast, fun = "sum", na.rm=TRUE)

hot_12 <- ifel(hotentrpoppoulcat_albers == 12, 1, 0)
names(hot_12) <- 'hot_12'
resultregion_12 <- zonal(hot_12, regionRast, fun = "sum", na.rm=TRUE)

hot_13 <- ifel(hotentrpoppoulcat_albers == 13, 1, 0)
names(hot_13) <- 'hot_13'
resultregion_13 <- zonal(hot_13, regionRast, fun = "sum", na.rm=TRUE)

hot_14 <- ifel(hotentrpoppoulcat_albers == 14, 1, 0)
names(hot_14) <- 'hot_14'
resultregion_14 <- zonal(hot_14, regionRast, fun = "sum", na.rm=TRUE)

hot_15 <- ifel(hotentrpoppoulcat_albers == 15, 1, 0)
names(hot_15) <- 'hot_15'
resultregion_15 <- zonal(hot_15, regionRast, fun = "sum", na.rm=TRUE)

hot_16 <- ifel(hotentrpoppoulcat_albers == 16, 1, 0)
names(hot_16) <- 'hot_16'
resultregion_16 <- zonal(hot_16, regionRast, fun = "sum", na.rm=TRUE)

hot_17 <- ifel(hotentrpoppoulcat_albers == 17, 1, 0)
names(hot_17) <- 'hot_17'
resultregion_17 <- zonal(hot_17, regionRast, fun = "sum", na.rm=TRUE)

all_results <- Reduce(function(x, y) merge(x, y, by="EU753", all=TRUE),
                      list(resultregion_12, resultregion_14, resultregion_11, resultregion_16,
                           resultregion_15, resultregion_13, resultregion_17))

all_results$all <- all_results$hot_12 + all_results$hot_14 + all_results$hot_11+
  all_results$hot_16 + all_results$hot_15 + all_results$hot_13 + all_results$hot_17
print(all_results)

fwrite(all_results, '/root/result/Hot_data/country_hot.csv')







######—— Fig.3###############
#Fig.3 a
countrys <- vect ('/root/autodl-tmp/WAEdata_new_y/result779/otherdata/Con_EU_dis.shp')
library(rnaturalearth)
coast <- ne_coastline(scale = "small", returnclass = "sf") %>% vect()
crs <- '+proj=longlat +datum=WGS84'
hotentrpoppoulcat_df<-as.data.frame(hotentrpoppoulcat,xy=T) 
names(hotentrpoppoulcat_df)<-c("x","y","hot")
unique(hotentrpoppoulcat_df$hot)

custom_colors <- c(#"0" = "white", "1" = "white","2" = "white","3" = "white","4" = "white",
  #"5" = "white","6" = "white","7" = "white","10" = "white",
  "0" = "lightgrey", "1" = "lightgrey","2" = "lightgrey","3" = "lightgrey","4" = "lightgrey",
  "5" = "lightgrey","6" = "lightgrey","7" = "lightgrey","10" = "#808080",
  "11" = "#FE97A4", "12" = "#71A3F1", "14" = "#EBBF00",  "13" = "#EA68A2", 
  "15" = "#EB7100", "16" = "#9BB31D", "17" = "#E53341")
labels <- c(#"0" = " ", "1" = " ","2" = " ","3" = " ","4" = " ","5" = " ","6" = " ","7" = " ","10" = " ",
  "12" = "Poultry to WAE", 
  "11" = "Human to WAE", 
  "14" = "Cattle to WAE", 
  "16" = "Cattle-Poultry to WAE", 
  "13" = "Human-Poultry to WAE", 
  "15" = "Human-Cattle to WAE", 
  "17" = "Human-Poultry-Cattle to WAE")

hotentrpoppoulcat_df$hot<-factor(hotentrpoppoulcat_df$hot)
ggplot(hotentrpoppoulcat_df) +
  geom_tile(aes(x = x, y = y, fill = hot)) +
  scale_fill_manual(values = custom_colors, 
                    na.value = "white",  # 
                    breaks = c("12", "14", "11", "16", "15", "13", "17"),  # 
                    labels = labels) +  # 
  geom_spatvector(data = coast, fill = NA) +
  geom_spatvector(data = countrys, color= 'black', fill = NA) +
  coord_sf(crs = crs, xlim = c(-160, 165), ylim = c(-56, 90)) +
  labs(x = NULL, y = NULL) +
  theme_bw() +
  theme(
    legend.position = "none",  # 显示图例
    legend.text = element_text(size= 12),
    axis.text = element_text(size = 12),
    panel.grid.major = element_blank()
  )

#Fig.3 b+c
#plot in excel


######(4) Dominant waterbird functional groups ###########
speciesPixelNumPath2 <- list.files('/root/singleRast_779',pattern = '.tif',full.names = T)

# WAE of different waterbird functional groups (Comfirmed Host)
allDf2 <- fread(paste0(basePath,'allDf779_reclass.csv'))
allDf2CH <- allDf2[allDf2$`Host`== "Confirmed Host",]
spName <- basename(speciesPixelNumPath2) %>% str_sub(.,1,-5)

for (fuc in unique(allDf2CH$`Functional Group`)) {
  fucName <- allDf2CH[allDf2CH$`Functional Group`==fuc,]
  paths <- speciesPixelNumPath2[spName%in%fucName$LatName]
  monthNum <- rast(paths) %>% sum(na.rm=T)
  calEntropy <- lapply(paths, function(x){
    r <- rast(x) %>% sum(.,na.rm=T)
    pi <- r/monthNum
    y <- -pi*log(pi)
    names(y) <- str_sub(basename(x),1,-5)
    return(y)
  })
  actEntropy <- rast(calEntropy) %>% sum(na.rm = T)
  writeRaster(actEntropy,paste0('/root/result779/WAE_',fuc,'_CH.tif'))
}
#(Comfirmed Host)
listFuncCH <- list.files('/root/autodl-tmp/WAEdata_new_y/result779/WAE_data/', pattern = 'CH', full.names = T)
AE_FuncCH<-rast(listFuncCH)%>% resample(globalRaster)%>% mask(globalCountry) %>% crop(ext(globalRaster))

pop<-rast("/root/autodl-tmp/WAEdata_new_y/result779/otherdata/pop.tif") %>% resample(globalRaster)%>% mask(globalCountry) %>% crop(ext(globalRaster))
poul<-rast("/root/autodl-tmp/WAEdata_new_y/result779/otherdata/poultry.tif")%>% resample(globalRaster)%>% mask(globalCountry) %>% crop(ext(globalRaster))
cattle<-rast("/root/autodl-tmp/WAEdata_new_y/result779/otherdata/cattle.tif")%>% resample(globalRaster)%>% mask(globalCountry) %>% crop(ext(globalRaster))

pop_log<- rast("/root/autodl-tmp/WAEdata_new_y/result779/otherdata/pop_log.tif")%>% resample(globalRaster)%>% mask(globalCountry) %>% crop(ext(globalRaster))
poul_log<- rast("/root/autodl-tmp/WAEdata_new_y/result779/otherdata/poul_log.tif")%>% resample(globalRaster)%>% mask(globalCountry) %>% crop(ext(globalRaster))
cattle_log<- rast("/root/autodl-tmp/WAEdata_new_y/result779/otherdata/cattle_log.tif")%>% resample(globalRaster)%>% mask(globalCountry) %>% crop(ext(globalRaster))
values(pop_log)[is.infinite(values(pop_log))] <- NA
values(poul_log)[is.infinite(values(poul_log))] <- NA
values(cattle_log)[is.infinite(values(cattle_log))] <- NA

countrys <- vect ('/root/autodl-tmp/WAEdata_new_y/result779/otherdata/Con_EU_dis.shp')
Con_5Raster <- subset(countrys, countrys$EU753 %in% c("China", "India", "EU", "Ethiopia", "United States"))%>% rasterize(.,globalRaster,field='EU753')

Global_CHEntropy_df <-c(AE_FuncCH,pop_log,poul_log,cattle_log,pop,poul,cattle,Con_5Raster) %>% as.data.frame(.,xy=T)

names(Global_CHEntropy_df) <- c('x','y','Large wading birds','Others','Seabirds','Shorebirds','Waterfowl','pop_log','poul_log','cattle_log','pop','poul','cattle',"Country")
Global_CHEntropy_df <- Global_CHEntropy_df%>%
  dplyr::select('x','y','Seabirds','Shorebirds','Waterfowl','Large wading birds','Others','pop_log','poul_log','cattle_log','pop','poul','cattle',"Country")
summary(Global_CHEntropy_df)
fwrite(Global_CHEntropy_df,'/root/result/Global_CHEntropy_df.csv')



## WAE of different waterbird functional groups (All Suspected Host)
allDf2 <- fread(paste0(basePath,'allDf779_reclass.csv'))
spName <- basename(speciesPixelNumPath2) %>% str_sub(.,1,-5)

for (fuc in unique(allDf2$`Functional Group`)) {
  fucName <- allDf2[allDf2$`Functional Group`==fuc,]
  paths <- speciesPixelNumPath2[spName%in%fucName$LatName]
  monthNum <- rast(paths) %>% sum(na.rm=T)
  calEntropy <- lapply(paths, function(x){
    r <- rast(x) %>% sum(.,na.rm=T)
    pi <- r/monthNum
    y <- -pi*log(pi)
    names(y) <- str_sub(basename(x),1,-5)
    return(y)
  })
  actEntropy <- rast(calEntropy) %>% sum(na.rm = T)
  writeRaster(actEntropy,paste0('/root/result779/WAE_',fuc,'.tif'))
}
# (All Suspected Host)
listFunc <- list.files('/root/autodl-tmp/WAEdata_new_y/result779/WAE_data/', pattern = 'ALL', full.names = T)
AE_Func<-rast(listFunc)%>% resample(globalRaster)%>% mask(globalCountry) %>% crop(ext(globalRaster))

pop<-rast("/root/autodl-tmp/WAEdata_new_y/result779/otherdata/pop.tif") %>% resample(globalRaster)%>% mask(globalCountry) %>% crop(ext(globalRaster))
poul<-rast("/root/autodl-tmp/WAEdata_new_y/result779/otherdata/poultry.tif")%>% resample(globalRaster)%>% mask(globalCountry) %>% crop(ext(globalRaster))
cattle<-rast("/root/autodl-tmp/WAEdata_new_y/result779/otherdata/cattle.tif")%>% resample(globalRaster)%>% mask(globalCountry) %>% crop(ext(globalRaster))

pop_log<- rast("/root/autodl-tmp/WAEdata_new_y/result779/otherdata/pop_log.tif")%>% resample(globalRaster)%>% mask(globalCountry) %>% crop(ext(globalRaster))
poul_log<- rast("/root/autodl-tmp/WAEdata_new_y/result779/otherdata/poul_log.tif")%>% resample(globalRaster)%>% mask(globalCountry) %>% crop(ext(globalRaster))
cattle_log<- rast("/root/autodl-tmp/WAEdata_new_y/result779/otherdata/cattle_log.tif")%>% resample(globalRaster)%>% mask(globalCountry) %>% crop(ext(globalRaster))
values(pop_log)[is.infinite(values(pop_log))] <- NA
values(poul_log)[is.infinite(values(poul_log))] <- NA
values(cattle_log)[is.infinite(values(cattle_log))] <- NA

countrys <- vect ('/root/autodl-tmp/WAEdata_new_y/result779/otherdata/Con_EU_dis.shp')
Con_5Raster <- subset(countrys, countrys$EU753 %in% c("China", "India", "EU", "United States"))%>% rasterize(.,globalRaster,field='EU753')

Global_birdEntropy_df <- c(AE_Func,pop_log,poul_log,cattle_log,pop,poul,cattle,Con_5Raster) %>% as.data.frame(.,xy=T)

names(Global_birdEntropy_df) <- c('x','y','Large wading birds','Others','Seabirds','Shorebirds','Waterfowl','pop_log','poul_log','cattle_log','pop','poul','cattle',"Country")
Global_birdEntropy_df <- Global_birdEntropy_df%>%
  dplyr::select('x','y','Seabirds','Shorebirds','Waterfowl','Large wading birds','Others','pop_log','poul_log','cattle_log','pop','poul','cattle',"Country")
summary(Global_birdEntropy_df)
fwrite(Global_birdEntropy_df,'/root/result/Global_birdEntropy_df.csv')



#Person r
Contry5c <- c("China", "India", "EU", "United States")

Global_birdEntropy_df5_long <- gather(Global_birdEntropy_df, key = "Variable", value = "Value",
                                      -`x`,-`y`,-`pop_log`, -`poul_log`, -`cattle_log`,-`pop`, -`poul`, -`cattle`, -`Country`) %>%na.omit()
Global_birdEntropy_df5_long <- subset(Global_birdEntropy_df5_long, Value>0)
summary(Global_birdEntropy_df5_long)#ALL

Global_CHEntropy_df5_long <- gather(Global_CHEntropy_df, key = "Variable", value = "Value",
                                    -`x`,-`y`,-`pop_log`, -`poul_log`, -`cattle_log`,-`pop`, -`poul`, -`cattle`, -`Country`) %>%na.omit()
Global_CHEntropy_df5_long <- subset(Global_CHEntropy_df5_long, Value>0)
summary(Global_CHEntropy_df5_long)#CH

#
All5_poppoulcatresults_df <- data.frame(Variable = character(),  b = numeric(),  SE = numeric(), b_SE = character(), z = numeric(),  P = numeric(),  R2 = numeric(), 
                                        Country = character(), Type = character(), Part = character()  ,stringsAsFactors = FALSE)
Contry5c <- c("China", "India", "EU", "United States")
parts <- c('pop','poul','cattle')
types <- c('All', 'CH')
variables <- c('Seabirds', 'Shorebirds', 'Waterfowl', 'Large wading birds', 'Others')

for(con in Contry5c){#
  
  for(typ in types){
    if(typ == 'All'){#all
      
      for(part in parts){#
        
        for (var in variables) {#
          
          con_birdEntropy_df <- subset(Global_birdEntropy_df5_long,Global_birdEntropy_df5_long$Country==con)
          vardata = subset(con_birdEntropy_df,con_birdEntropy_df$Variable==var)
          
          if(part =='pop'){#
            model <- lm(`pop_log` ~ Value, vardata)
          }
          if(part =='poul'){#
            model <- lm(`poul_log` ~ Value, vardata)
          }
          if(part =='cattle'){#cattle
            model <- lm(`cattle_log` ~ Value, vardata)
            
          }
          model_summary <- summary(model)
          
          b <- model_summary$coefficients[2, 1]
          SE <- model_summary$coefficients[2, 2]
          b_SE <- sprintf("%.2f ± %.2f", b, SE)
          z <- b / SE
          P <- model_summary$coefficients[2, 4]
          R2 <- model_summary$r.squared
         
          All5_poppoulcatresults_df <- rbind(All5_poppoulcatresults_df, data.frame(Variable = var, b = b, SE = SE, b_SE = b_SE, z = z, P = P, R2 = R2, 
                                                                                   Country = con, Type = typ, Part = part))
        }
      }
    }
    
    if(typ == 'CH'){#CH
      
      for(part in parts){#
        
        for (var in variables) {#
          
          con_CHEntropy_df <- subset(Global_CHEntropy_df5_long,Global_CHEntropy_df5_long$Country==con)
          vardata = subset(con_CHEntropy_df,con_CHEntropy_df$Variable==var)
          
          if(part =='pop'){#
            model <- lm(`pop_log` ~ Value, vardata)
          }
          if(part =='poul'){#
            model <- lm(`poul_log` ~ Value, vardata)
          }
          if(part =='cattle'){#cattle
            model <- lm(`cattle_log` ~ Value, vardata)
            
          }
          model_summary <- summary(model)
          
          # 提取所需的参数
          b <- model_summary$coefficients[2, 1]
          SE <- model_summary$coefficients[2, 2]
          b_SE <- sprintf("%.2f ± %.2f", b, SE)
          z <- b / SE
          P <- model_summary$coefficients[2, 4]
          R2 <- model_summary$r.squared
          
          All5_poppoulcatresults_df <- rbind(All5_poppoulcatresults_df, data.frame(Variable = var, b = b, SE = SE, b_SE = b_SE, z = z, P = P, R2 = R2, 
                                                                                   Country = con, Type = typ, Part = part))
        }
        
      }
      
    }
    
  }
}

All5_poppoulcatresults_df
fwrite(All5_poppoulcatresults_df,'/root/result/All5_poppoulcatresults_df.csv')




#Moran's I
FuncTypepath <- list.files('/root/autodl-tmp/WAEdata_new_y/result779/WAE_data/', full.names = T)
countrys <- vect ('/root/autodl-tmp/WAEdata_new_y/result779/otherdata/Con_EU_dis.shp')

crs <- '+proj=longlat +datum=WGS84'
globalRaster <- rast(vals=1:259200,nrows=360, ncols=720,xmin=-180, xmax=180,ymin=-90, ymax=90,crs=crs)
popd2015<-rast("/root/autodl-tmp/WAEdata_new_y/result779/otherdata/pop.tif") %>% resample(globalRaster)
poul2015<-rast("/root/autodl-tmp/WAEdata_new_y/result779/otherdata/poultry.tif")%>% resample(globalRaster)
cattle2015<-rast("/root/autodl-tmp/WAEdata_new_y/result779/otherdata/cattle.tif") %>% resample(globalRaster)

DensityR <- c(popd2015, poul2015, cattle2015)

#
library(spdep)
Moran_Func <- data.table()

for(type in FuncTypepath){
  FuncTypeName <- basename(type) %>% str_sub(1,-5)
  FuncTypeR <- rast(type)
  FuncTypeR2 <- ifel(FuncTypeR>0,FuncTypeR,NA)
  names(FuncTypeR2) <- FuncTypeName
  
  for (DensityName in names(DensityR)) {
    r <- DensityR[DensityName]
    r2 <- c(r,FuncTypeR2)
    #
    countryName <- c("China", "India", "EU", "United States")
    for (n in countryName) {
      countrShp <- subset(countrys,countrys$EU753==n)
      countrR <- mask(r2,countrShp)
      df <- as.data.frame(countrR,xy=T,na.rm=T)
      # 
      coords <- cbind(df$x, df$y)
      # 
      knn <- knearneigh(coords, k =5)
      nb_list <-knn2nb(knn)
      # 
      W_std <-nb2listw(nb_list, style ="W", zero.policy =TRUE)
      res_xy <-  moran_bv(df[,3],df[,4],W_std, nsim=499)
      
      moranSD <- boot::boot.ci(res_xy, conf=c(0.99, 0.95, 0.9), type="basic")
      tmpDf <- data.table(
        countrName=n,
        Func_type=FuncTypeName,
        Density_name=DensityName,
        moranIndex=res_xy$t0,
        sd_1=moranSD$basic[2,4],
        sd_2=moranSD$basic[2,5]
      )
      Moran_Func <- rbind(Moran_Func,tmpDf)
      print(n)
    }
    
  }
}

fwrite(Moran_Func,'/root/result/Moran_Func.csv')


######—— Fig.4###############
#Fig.4 a
Con_popentrpoul_sf <- vect ('/root/autodl-tmp/WAEdata_new_y/result779/worldBorder/continentNew.shp');head(Con_popentrpoul_sf)

Con_Africa<-subset(Con_popentrpoul_sf,Con_popentrpoul_sf$CONTINENT=="Africa")

ggplot() +
  geom_sf(data = world.map, fill = "lightgrey", color = NA) + 
  geom_sf(data = Con_Africa, fill = "#baa1a1", color = NA) +    #
  coord_sf(crs = crs, xlim = c(-65, 60), ylim = c(-50, 50)) +
  theme_bw() + # 使用简洁主题
  theme(legend.position = "none",
        axis.text = element_text(size=12),
        panel.grid.major = element_blank() )


#【US pie】
piedata_US <- lineardata_max[lineardata_max$Country == "United States", ]
piedata_US$Value<-c(1,1,1,1,1,1)

color_US<-c("pop_All"="#B42D34","pop_CH"="#B42D34","poul_All"="#D9705A","poul_CH"="#D9705A","cattle_All"="#D9705A","cattle_CH"="#F2B396")

# 
desired_order <- c("pop_CH","cattle_All","cattle_CH", "poul_All", "poul_CH","pop_All")
piedata_US$Order <- factor(paste(piedata_US$Part, piedata_US$Type, sep = "_"), levels = desired_order)
piedata_US <- piedata_US[order(piedata_US$Order), ]
piedata_US$Color <- factor(piedata_US$Color, levels = unique(piedata_US$Color))


p <- ggplot(piedata_US, aes(x = "", y = Value, fill = Order)) +
  geom_bar(stat = "identity", width = 1) +
  coord_polar("y", start = 0) +
  scale_fill_manual(values = color_US) +
  theme_void() +
  theme(
    panel.background = element_rect(fill = "transparent"), # 
    plot.background = element_rect(fill = "transparent", color = NA), # 
    legend.position = "none"
  )
p
# 
ggsave("/root/autodl-tmp/WAEdata_new_y/result779/piedata_US.png", p, bg = "transparent", width = 10, height = 8, units = "in")


#Fig.4 b
#【US boxplot】
US_birdEntropy_df<-subset(Global_birdEntropy_df,Global_birdEntropy_df$Country=='United States')

US_birdEntropy_df_long <- gather(US_birdEntropy_df, key = "Variable", value = "Value",
                                 -`x`,-`y`,-`pop_log`, -`poul_log`, -`cattle_log`,-`pop`, -`poul`, -`cattle`, -`Country`) %>%na.omit()
US_birdEntropy_df_long <- subset(US_birdEntropy_df_long, Value>0)
summary(US_birdEntropy_df_long)#ALL

x2<-factor(US_birdEntropy_df_long$Variable,levels = c('Others','Large wading birds','Waterfowl','Shorebirds','Seabirds')) #

ggplot(US_birdEntropy_df_long, aes(x = x2, y =Value,fill=x2)) + #,
  scale_fill_manual(values = c("#FF7F00", "#984EA3", "#E41A1C", "#4DAF4A", "#377EB8")) + 
  #scale_fill_manual(values = c("#377EB8", "#E41A1C", "#4DAF4A", "#984EA3", "#FF7F00")) +
  geom_hline(yintercept = mean(US_birdEntropy_df_long$Value), size = 1, color = "grey", linetype = "dashed")+
  geom_violin(alpha=0.5,position = position_dodge(width = 0.01), scale = 'width') +
  geom_boxplot(alpha=0.8,width=0.45,position=position_dodge(width=0.1),size=0.75,outlier.colour = NA)+
  #scale_x_continuous(position = "top")+
  scale_y_continuous(position = "right")+
  labs(y = "", x = " ",title = "")+
  theme_bw()+
  theme(panel.grid = element_blank(), 
        panel.background = element_rect(fill = 'transparent', color = "grey"), 
        legend.position = 'none',
        axis.text.x = element_text(size = 66),
        axis.text.y = element_blank(), 
        text = element_text(size = 26),
  ) +
  coord_flip()#



