library(data.table)
library(ggplot2)
library(terra)
library(dplyr)
library(stringr)
library(paletteer)
library(rnaturalearth)
library(tidyterra)
library(patchwork)
basePath<-"/root/autodl-tmp/root/autodl-tmp/WAEdata_new_y/result/"
crs <- '+proj=longlat +datum=WGS84'
`%notin%` <- Negate(`%in%`)

world.map <- rnaturalearth::ne_countries(returnclass = "sf") |> dplyr::filter(continent != "Antarctica")
globalCountry <- vect(world.map) 
globalRaster <- rast(vals=1:259200,nrows=360, ncols=720,xmin=-180, xmax=180,ymin=-90, ymax=90,crs=crs)
coast <- rnaturalearth::ne_coastline(scale = "small", returnclass = "sf")

allDf <- fread('/root/autodl-tmp/root/autodl-tmp/WAEdata_new_y/result/otherdata/allDf779_reclass.csv')
speciesPixelNumPath <- list.files('/root/autodl-tmp/root/autodl-tmp/WAEdata_new_y/singleRast',pattern = '.tif',full.names = T)

spName <- basename(speciesPixelNumPath) %>% str_sub(.,1,-5)                      
speciesPixelNumPath5 <- speciesPixelNumPath[spName%in%allDf$LatName]




#####Fig1.Species Richness CV & AE###############
#######a-----------
#Species Richness CV calculation

spMonth <- rast(speciesPixelNumPath5)
sp_perMonth <- tapp(spMonth,names(spMonth),sum,na.rm=T) %>% mask(.,globalCountry)
sp_perMonth <- sp_perMonth[[c('X1','X2','X3','X4','X5','X6','X7','X8','X9','X10','X11','X12')]]
#writeRaster(sp_perMonth,'/root/autodl-tmp/root/autodl-tmp/WAEdata_new_y/result/sp_perMonth.tif')

sp_perMonth <- rast('/root/autodl-tmp/root/autodl-tmp/WAEdata_new_y/result/sp_perMonth.tif')
sp_perMonth

sp_sd<-app(sp_perMonth, sd, na.rm = T)
#sp_sd1<-ifel(sp_sd>1, sp_sd, NA)
#plot(sp_sd1)

sp_mean<-app(sp_perMonth, mean, na.rm = T)
sp_meanthan1<-ifel(sp_mean>0.2, sp_mean, NA)
#plot(sp_meanthan1)

sp_cv <- sp_sd/sp_meanthan1*100
#sp_cv <- sp_sd/sp_mean*100
plot(sp_cv)
#writeRaster(sp_cv,paste0('/root/autodl-tmp/root/autodl-tmp/WAEdata_new_y/result/sp_cv.tif'))

#Polt
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
    # panel.background = element_rect(fill = "white"),#
    #legend.direction = "horizontal",
    legend.title = element_blank(),
    #legend.title = element_text(hjust=0.5),
    legend.title.align = -10,
    legend.position = c(0.15, 0.05),#
    legend.direction='horizontal',#
    legend.key.width = unit(1.2,'cm'), #
    legend.key.height = unit(0.3,'cm')
    
  )
pSRCV



#sp_perMonth <- rast(paste0(basePath,'result/sp_perMonth2.tif')) #%>% as.data.table
#sp_df <- sp_perMonth[[c('Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec')]]
sp_df <- sp_perMonth[[c('X1','X2','X3','X4','X5','X6','X7','X8','X9','X10','X11','X12')]]
#names(sp_df) <- c('Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec')

sp_df <- terra::as.data.frame(sp_df, xy = TRUE, na.rm = FALSE) 
head(sp_df)

# 
library(tidyr)
sp_df_long <- gather(sp_df, key = "Month", value = "Value", -x, -y,na.rm = T)
names(sp_df_long)<-c("lon","lat","Month","Value")
head(sp_df_long)
sp_df_long$Month <- factor(sp_df_long$Month, levels = c('X1','X2','X3','X4','X5','X6','X7','X8','X9','X10','X11','X12'))
#sp_df_long$Month <- factor(sp_df_long$Month, levels = c("Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"))
#数据
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
  # labs(y="Latitude(°)",)+
  #geom_ribbon(aes(y=lat, xmin=row_means-row_sds,xmax = row_means+row_sds),fill = "lightgrey", alpha=0.5)+
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

#lit(Species Richness)----
#Species Richness
spNum <- lapply(speciesPixelNumPath5, function(x){
  r <- rast(x) %>% sum(.,na.rm=T)
  r2 <- ifel(r>0,1,NA)
  names(r2) <- str_sub(basename(x),1,-5)
  return(r2)
})
spNum <- rast(spNum)
spNumTif <- sum(spNum,na.rm=T) %>% mask(globalCountry)
terra::writeRaster(spNumTif,'/root/autodl-tmp/root/autodl-tmp/WAEdata_new_y/result/spNumTif.tif')

spNumTif<-rast(paste0('/root/autodl-tmp/root/autodl-tmp/WAEdata_new_y/result/spNumTif.tif'))
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
    legend.position = 'bottom',#
    legend.direction='horizontal',#
    legend.key.width = unit(1.15,'cm'), #
    legend.key.height = unit(0.3,'cm')
  )
pspNum

#
Vborder <- vect('/root/autodl-tmp/allData/Wallace_zoogeographic/newValisBorder.shp') 
vRaster <- rasterize(Vborder,spNumTif,field='name')
valisEntropy <- c(spNumTif,vRaster) %>% terra::as.data.frame() %>%na.omit()


###plot
valisEntropy <- valisEntropy[valisEntropy$name != "Antarctica",]   #
Vborder1<-factor(valisEntropy$name,levels = c("Australian","Ethiopian","Nearctic","Neotropical","Oriental","Palearctic")) #

ggplot(valisEntropy, aes(x = Vborder1, y = sum,fill=Vborder1)) + #,
  scale_fill_manual(values = c("#E0312D","#F6BA26","#8204FB","#3ca1f9","#45F3DD","#49FB18"))+
  geom_hline(yintercept = mean(valisEntropy$sum), size = 1, color = "grey", linetype = "dashed")+
  geom_violin(alpha=0.5,position = position_dodge(width = 0.01), scale = 'width') +
  geom_boxplot(alpha=0.8,width=0.45,position=position_dodge(width=0.1),size=0.75,outlier.colour = NA)+
  theme(panel.grid = element_blank(), 
        panel.background = element_rect(fill = 'transparent', color = "black"), 
        legend.position = 'none',
        axis.text = element_text(size = 30),
        text = element_text(size = 34),) +
  labs(y = "spNumTif", x = " ")

mean(valisEntropy$sum)
Mean_result <- valisEntropy %>%
  group_by(name) %>%
  summarise(across(everything(), mean, na.rm = TRUE))
Mean_result




#b-----------
basePath <- '/root/autodl-tmp/root/autodl-tmp/WAEdata_new_y/result/'
AE<-rast(paste0(basePath,'WAE.tif')); plot(AE)

library(paletteer)
pAE <- ggplot() +
  geom_spatraster(data = AE) +
  geom_spatvector(data=coast,fill=NA)+coord_sf(crs = crs,xlim=c(-160,165),ylim=c(-56,90))+
  theme_bw()+
  #annotation_scale(location='bl') +
  # spatial-aware automagic north arrow
  #annotation_north_arrow(location = "tl", which_north = "true",
  #                      style = north_arrow_fancy_orienteering)+
  scale_fill_gradientn(colours = paletteer_c("grDevices::Zissou 1", 30) ,na.value='white',values = c(0, 0.75, 1))+
  #scale_fill_brewer(palette = "RdYlGn")+
  labs(fill="Species Richness")+
  theme(
    axis.text = element_text(size=12),
    plot.title = element_text(hjust=0.5),
    # axis.line = element_line(color = 'black'),
    # panel.background = element_rect(fill = "white"),#
    #legend.direction = "horizontal",
    legend.title = element_blank(),
    #legend.title = element_text(hjust=0.5),
    legend.title.align = -10,
    legend.position = c(0.15, 0.05),#
    legend.direction='horizontal',#
    legend.key.width = unit(1.15,'cm'), #
    legend.key.height = unit(0.3,'cm')
  )
pAE


library(ggpointdensity)
AE_df <- terra::as.data.frame(AE, xy = TRUE, na.rm = T) 
head(AE_df)
names(AE_df)<-c("lon","lat","sum")

pAE2<-ggplot(AE_df,aes(x=sum, y=lat),)+
  geom_pointdensity(adjust = 4, size = 0.5)+                              # adjust：
  scale_color_distiller(palette = "Spectral", direction = -1)+ # 
  xlab('WAE')+
  xlim(c(0,5.5))+
  ylim(c(-56,90))+
  theme_bw()+
  # labs(y="Latitude(°)",)+
  #geom_ribbon(aes(y=lat, xmin=row_means-row_sds,xmax = row_means+row_sds),fill = "lightgrey", alpha=0.5)+
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


#lit（AE）--------
#
Vborder <- vect('/root/autodl-tmp/Wallace_zoogeographic/newValisBorder.shp') 
vRaster <- rasterize(Vborder,AE,field='name')
valisEntropy <- c(AE,vRaster) %>% terra::as.data.frame() %>%na.omit()


###plot
valisEntropy <- valisEntropy[valisEntropy$name != "Antarctica",]   #
#Vborder1<-factor(valisEntropy$name,levels = c("Australian","Ethiopian","Nearctic","Neotropical","Oriental","Palearctic")) #
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










