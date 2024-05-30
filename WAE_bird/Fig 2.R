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
AE <- rast('/root/autodl-tmp/humPoulResult/data/AE_data/AE.tif') 
sample_TP <- terra::extract(AE,TP_points,mean,bind=T) %>% as.data.frame()
sample_TP$label <- 'good'

trueNgRast <- mask(AE,breakRast,inverse=T)

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
threshold <- dfroc1$auc

WAE<- rast(paste0(basePath,'WAE_data/WAE.tif'))%>% mask(globalCountry)
plotWAE <- ifel(WAE>=threshold,1,0)   #4.19

coast <- rnaturalearth::ne_coastline(scale = "small", returnclass = "sf")
crs <- '+proj=longlat +datum=WGS84'
ggplot() +
  geom_spatraster(data = plotWAE) +
  geom_spatvector(data=coast,fill=NA)+
  #coord_sf(crs = crs,xlim=c(-160,180),ylim=c(-56,90))+
  geom_spatvector(data=outBreak2,size=2, shape = 4,color='#ff4f4c',fill='#f26c6a', alpha = 0.8,stroke = 0.3)+   #shape = 21是圆圈，shape = 4是×
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

#b、c###########
# validation among contries
countrRast <- rasterize(globalSHP,AE,'NAME_LONG')
Country <- terra::extract(countrRast,sample_TN_vect)
sample_TN$Country <- Country$NAME_LONG
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

#plot
attrData <- fread('path')
result_plot <- left_join(result_country,attrData)
#b
pop_quantiles <- quantile(result_df4$POP_EST, probs=seq(0, 1, by=0.1), na.rm = TRUE)

result_plot$POP_EST_size <- findInterval(result_plot$POP_EST, vec = pop_quantiles)
result_plot$POP_EST <- as.numeric(result_plot$POP_EST)
result_plot <- mutate(result_plot, POP_EST_label = cut(POP_EST, breaks = c(0, 10000000, 20000000, 50000000, 100000000, 500000000, 1000000000, 1397800000),
                                                     labels = c("<1000", "1000 - 2000", "2000-5000", "5000 - 10000", "10000 - 50000", "50000 - 100000", ">100000")))


result_plot <- result_plot[order(result_plot$CONTINENT, -result_plot$POP_EST), ]
country5 <- c("China", "United States", "India", "Nigeria")
result_dfc <- result_plot[result_plot$country %in% country5]
result_dfc <- result_dfc %>%
  mutate(country = ifelse(country == "United States", "USA", country))

size_values <- c("<1000" = 3, "1000 - 2000" = 4, "2000-5000" = 6, "5000 - 10000" = 8, "10000 - 50000" = 10, "50000 - 100000" = 12, ">100000" = 18)

ggplot(result_df5, aes(per, accuracy, color=CONTINENT)) +
  geom_point(aes(size=POP_EST_label), alpha=0.5) +
  ylim(c(0,1))+
  geom_text(data=result_dfc, aes(label=country, x=per, y=accuracy),size=4.5, hjust=0.5, vjust=2) +
  scale_color_manual(values = c("#984EA3", "#E41A1C","#4DAF4A",  "#FF7F00", "#a38900","#377EB8"))+
  #geom_smooth(method = "gam", se = T, fill="lightgrey",color = "grey", alpha=0.3, linetype = "dashed") +  # 全部点的趋势线，黑色  #loess  gam
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

#c
pop_quantiles <- quantile(result_plot$POP_EST, probs=seq(0, 1, by=0.1), na.rm = TRUE)
result_plot$POP_EST_size <- findInterval(result_plot$POP_EST, vec = pop_quantiles)
result_plot$POP_EST <- as.numeric(result_plot$POP_EST)
result_plot <- mutate(result_plot, POP_EST_label = cut(POP_EST, breaks = c(0, 10000000, 20000000, 50000000, 100000000, 500000000, 1000000000, 1397800000),
                                                     labels = c("<1000", "1000 - 2000", "2000-5000", "5000 - 10000", "10000 - 50000", "50000 - 100000", ">100000")))

result_plot <- result_plot[order(result_plot$CONTINENT, -result_plot$POP_EST), ]
country5 <- c("China", "United States", "India", "Nigeria")
result_dfc <- result_plot[result_plot$country %in% country5]
result_dfc <- result_dfc %>%
  mutate(country = ifelse(country == "United States", "USA", country))

size_values <- c("<1000" = 3, "1000 - 2000" = 4, "2000-5000" = 6, "5000 - 10000" = 8, "10000 - 50000" = 10, "50000 - 100000" = 12, ">100000" = 18)
p2<-ggplot(result_plot, aes(log(GDP_MD), accuracy, color=CONTINENT)) +
  geom_point(aes(size=POP_EST_label), alpha=0.5,show.legend = F) +
  geom_text(data=result_dfc, aes(label=country, x=log(GDP_MD), y=accuracy),show.legend = F,size=4.5, hjust=0.5, vjust=2) +
  scale_color_manual(values = c("#984EA3", "#E41A1C","#4DAF4A",  "#FF7F00", "#a38900","#377EB8"))+
  #geom_smooth(method = "gam", se = T, fill="lightgrey",color = "grey", alpha=0.3, linetype = "dashed") + 
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



