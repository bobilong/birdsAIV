

library(terra)
library(data.table)
library(rnaturalearth)
library(cowplot)   
library(ggplot2)   
library(dplyr)
library(tidyterra)
library(ggpmisc)
library(ggpointdensity)
library(patchwork)
library(paletteer)
library(ggpointdensity)



######Fig. 1 Seasonal variation and activity entropy of waterbird species richness###############
basePath<-"/root/autodl-tmp/humPoulResult/data/"

world.map <- rnaturalearth::ne_countries(returnclass = "sf") |> dplyr::filter(continent != "Antarctica")
globalCountry <- vect(world.map) 
globalRaster <- rast(vals=1:259200,nrows=360, ncols=720,xmin=-180, xmax=180,ymin=-90, ymax=90,crs=crs)
coast <- ne_coastline(scale = "small", returnclass = "sf")
crs <- '+proj=longlat +datum=WGS84'

speciesPixelNumPath <- list.files('/root/autodl-tmp/humPoulResult/data/single_model',pattern = '.tif',full.names = T)

#######a-----------
#Species Richness CV

spMonth <- rast(speciesPixelNumPath)
sp_perMonth <- tapp(spMonth,names(spMonth),sum,na.rm=T) %>% mask(.,globalCountry)
#writeRaster(sp_perMonth,paste0(basePath,'WAE_data/sp_perMonth.tif'))

sp_sd<-app(sp_perMonth, sd, na.rm = T)
sp_mean<-app(sp_perMonth, mean, na.rm = T)

sp_cv <- sp_sd/sp_mean*100
plot(sp_cv)
#writeRaster(sp_cv,paste0(basePath,'WAE_data/sp_cv.tif'))

#Polt
pSRCV<-ggplot() +
  geom_spatraster(data = sp_cv) +  #sp_cv
  geom_spatvector(data=coast,fill=NA)+coord_sf(crs = crs,xlim=c(-160,165),ylim=c(-56,90))+
  theme_bw()+
  scale_fill_gradientn(colours = paletteer_c("grDevices::Zissou 1", 30) ,na.value='white',values = c(0, 0.3, 1))+
  #labs(title = "CV")+
  theme(
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


sp_df_avg <- rast(paste0(basePath,'WAE_data/sp_perMonth.tif')) #%>% as.data.table


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
    panel.grid = element_blank(),
    axis.ticks.length = unit(-4, "pt") #
  )

pSRCV+pSRCV2

#little ----
#Species Richness

spNum <- lapply(speciesPixelNumPath, function(x){
  r <- rast(x) %>% sum(.,na.rm=T)
  r2 <- ifel(r>0,1,NA)
  names(r2) <- str_sub(basename(x),1,-5)
  return(r2)
})
spNum <- rast(spNum)
spNumTif <- sum(spNum,na.rm=T) %>% mask(globalCountry)
#terra::writeRaster(spNumTif,paste0(basePath,'WAE_data/spNumTif.tif'))


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
    legend.position = 'none',#
    legend.direction='horizontal',#
    legend.key.width = unit(1.15,'cm'), #
    legend.key.height = unit(0.3,'cm')
  )
pspNum


#b-----------


##WAE calculation-------
#######1.Cumulative species months#########
speciesPixelNum <- rast(speciesPixelNumPath)
allMonth <- sum(speciesPixelNum,na.rm=T)%>% mask(globalCountry)
#writeRaster(allMonth,paste0(basePath,'WAE_data/allMonth.tif'), overwrite=T)

#######2.WAE##########
allMonth <- rast(paste0(basePath,'WAE_data/allMonth.tif'))
calEntropy <- lapply(speciesPixelNumPath, function(x){
  r <- rast(x) %>% sum(.,na.rm=T)
  pi <- r/allMonth
  y <- -pi*log(pi)
  names(y) <- str_sub(basename(x),1,-5)
  return(y)
})
calEntropy <- rast(calEntropy)
WAE<-sum(calEntropy,na.rm = T)%>% mask(globalCountry)
plot(WAE)
#writeRaster(WAE,paste0(basePath,'WAE_data/WAE.tif'), overwrite=T)


#Plot
#WAE<-rast(paste0(basePath,'WAE_data/WAE.tif')); plot(WAE)

pWAE <- ggplot() +
  geom_spatraster(data = WAE) +
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
pWAE



WAE_df <- terra::as.data.frame(WAE, xy = TRUE, na.rm = T) 
head(WAE_df)
names(WAE_df)<-c("lon","lat","sum")

pWAE2<-ggplot(WAE_df,aes(x=sum, y=lat),)+
  geom_pointdensity(adjust = 4, size = 0.5)+                              # adjust
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
    axis.title.y = element_blank(),
    axis.text.y=element_blank(),
    panel.grid = element_blank(),
    axis.ticks.length = unit(-4, "pt"), #
    legend.position = "none"
  )
pWAE2

pWAE+pWAE2


#little violin（WAE）--------
Vborder <- vect('/root/autodl-tmp/Wallace_zoogeographic/newValisBorder.shp') 
vRaster <- rasterize(Vborder,WAE,field='name')
valisEntropy <- c(WAE,vRaster) %>% terra::as.data.frame() %>%na.omit()


###plot
valisEntropy <- valisEntropy[valisEntropy$name != "Antarctica",]   #
Vborder1<-factor(valisEntropy$name,levels = c("Oriental","Ethiopian","Australian","Neotropical","Nearctic","Palearctic")) #

ggplot(valisEntropy, aes(x = Vborder1, y = sum,fill=Vborder1)) + #,
  scale_fill_manual(values = c("#45F3DD","#F6BA26","#E0312D","#3ca1f9","#8204FB","#49FB18"))+
  geom_hline(yintercept = mean(valisEntropy$sum), size = 1, color = "grey", linetype = "dashed")+
  geom_violin(alpha=0.5,position = position_dodge(width = 0.01), scale = 'width') +
  geom_boxplot(alpha=0.8,width=0.45,position=position_dodge(width=0.1),size=0.75,outlier.colour = NA)+
  theme(panel.grid = element_blank(), 
        panel.background = element_rect(fill = 'transparent', color = "black"), 
        legend.position = 'none',
        axis.text = element_text(size = 30),
        text = element_text(size = 34),) +
  labs(y = "WAE", x = " ")






