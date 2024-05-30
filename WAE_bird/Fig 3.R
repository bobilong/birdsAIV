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

library(devtools)
library(tricolore)
library(ggtern)


######Fig. 3 Global hotspots of waterbird carried AIV exposure risk considering waterbird activity entropy, human density, and poultry density.###############
basePath<-"/root/autodl-tmp/humPoulResult/data/"

world.map <- rnaturalearth::ne_countries(returnclass = "sf") |> dplyr::filter(continent != "Antarctica")
globalCountry <- vect(world.map) 
globalRaster <- rast(vals=1:259200,nrows=360, ncols=720,xmin=-180, xmax=180,ymin=-90, ymax=90,crs=crs)
coast <- ne_coastline(scale = "small", returnclass = "sf")
crs <- '+proj=longlat +datum=WGS84'


#a########
popd2015<-rast("/root/autodl-tmp/GWP_v4/popd2015_30.tif") %>% resample(globalRaster)%>% mask(globalCountry)
Entropy<-rast(paste0(basePath,'WAE_data/WWAE.tif'))%>% mask(globalCountry)
poul2015<-rast("/root/autodl-tmp/zyresult/Poultry_duckchic.tif")%>% resample(globalRaster)%>% mask(globalCountry)


quantiles_pop <- global(popd2015,quantile,probs=seq(0, 1, 0.1),na.rm=T)
quantiles_entr<- global(Entropy,quantile,probs=seq(0, 1, 0.1),na.rm=T)
quantiles_poul<- global(poul2015,quantile,probs=seq(0, 1, 0.1),na.rm=T)


spClass <- data.table()
for (i in 1:(length(quantiles_pop)-1)) {
  pop <- data.table(from=quantiles_pop[1,i],to=quantiles_pop[1,i+1],class=i)
  entropy <- data.table(from=quantiles_entr[1,i],to=quantiles_entr[1,i+1],class=i)
  poul <- data.table(from=quantiles_poul[1,i],to=quantiles_poul[1,i+1],class=i)
  
  tmp <- cbind(pop,entropy,poul)
  spClass <- rbind(spClass,tmp)
}

spclmat <- as.matrix(spClass)
class_pop <- classify(popd2015,spclmat[,1:3],include.lowest=TRUE)
class_entropy <- classify(Entropy,spclmat[,4:6],include.lowest=TRUE)
class_poul <- classify(poul2015,spclmat[,7:9],include.lowest=TRUE)

popentrpoul_df <- c(class_pop,class_entropy,class_poul)  %>% as.data.frame(xy=T)
popentrpoul_df2<-na.omit(popentrpoul_df)
#popentrpoul_df2 <- apply(popentrpoul_df, 2, tidyr::replace_na,0) %>% as.data.frame()
names(popentrpoul_df2)<-c("x","y","pop_re","entr_re","poul_re")


#color
col <- Tricolore(popentrpoul_df2, p1 = 'pop_re', p2 = 'poul_re', p3 = 'entr_re', breaks = Inf)  
col$key


#plot
plot_res <- ggplot(popentrpoul_df2) +
  geom_tile(aes(x = x, y = y, fill = rgb)) +
  scale_fill_identity() +
  geom_spatvector(data=coast,fill=NA)+coord_sf(crs = crs,xlim=c(-160,165),ylim=c(-56,90))+
  labs(x = NULL, y = NULL,) +
  #theme_minimal()+
  theme_bw()+
  theme(
    legend.position="none"
  )
print(plot_res)


plot_res+
  annotation_custom(ggplotGrob(col$key),
                    xmin = -170, xmax =-100, ymin = -60, ymax = 10)



#b############

globalSHP <- vect('/root/autodl-tmp/zyresult/Con_popentrpoul_sf_EU.shp')


quantiles_pop <- global(popd2015,quantile,probs=seq(0, 1, 0.01),na.rm=T)
quantiles_entr<- global(Entropy,quantile,probs=seq(0, 1, 0.01),na.rm=T)
quantiles_poul<- global(poul2015,quantile,probs=seq(0, 1, 0.01),na.rm=T)


quan_pop <- 50
quan_entr<- 4.19
quan_poul <- quantiles_poul[1,84]  #83q

pop_new <- ifel(popd2015>quan_pop,1,0)
entr_new <- ifel(Entropy>quan_entr,10,0)
poul_new <- ifel(poul2015>quan_poul,2,0)

hotentrpoppoul <- pop_new+entr_new+poul_new
hotAND <- ifel(hotentrpoppoul==13,1,0)  ; plot(hotAND)
hotOR_pop <- ifel(hotentrpoppoul==11|hotentrpoppoul==13,1,0) ; plot(hotOR_pop)
hotOR_poul <- ifel(hotentrpoppoul==12|hotentrpoppoul==13,1,0)  ; plot(hotOR_poul)


hotentrpoppoul_df<-as.data.frame(hotentrpoppoul,xy=T) 
names(hotentrpoppoul_df)<-c("x","y","hot")
ggplot(hotentrpoppoul_df) +
  geom_tile(aes(x = x, y = y, fill = factor(hot))) +
  scale_fill_manual(values = c("0" = "grey","1" = "grey","2" = "grey", "3" = "grey","10" = "grey","13" = "red","11" = "#f58a2c","12" = "yellow")) +  #b036fd  #ffa500  #3da3ef
  geom_spatvector(data=coast,fill=NA)+coord_sf(crs = crs,xlim=c(-160,165),ylim=c(-56,90))+
  labs(x = NULL, y = NULL,) +
  #theme_minimal()+
  theme_bw()+
  theme(
    legend.position="none"
  )





