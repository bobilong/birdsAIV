library(data.table)
library(ggplot2)
library(terra)
library(dplyr)
basePath<-"/root/autodl-tmp/humPoulResult/data/"

world.map <- rnaturalearth::ne_countries(returnclass = "sf") |> dplyr::filter(continent != "Antarctica")
globalCountry <- vect(world.map) 
globalRaster <- rast(vals=1:259200,nrows=360, ncols=720,xmin=-180, xmax=180,ymin=-90, ymax=90,crs=crs)
coast <- rnaturalearth::ne_coastline(scale = "small", returnclass = "sf")
crs <- '+proj=longlat +datum=WGS84'
allDf <- fread(paste0(basePath,'allDf786_reclass.csv'))
speciesPixelNumPath <- list.files('/root/autodl-tmp/humPoulResult/data/single_model',pattern = '.tif',full.names = T)
spName <- basename(speciesPixelNumPath) %>% str_sub(.,1,-5)                      
speciesPixelNumPath2 <- speciesPixelNumPath[spName%in%allDf$LatName]



#####Fig1.Species Richness CV & AE###############
#######a-----------
#Species Richness CV calculation

spMonth <- rast(speciesPixelNumPath2)
sp_perMonth <- tapp(spMonth,names(spMonth),sum,na.rm=T) %>% mask(.,globalCountry)
#writeRaster(sp_perMonth,paste0(basePath,'AE_data/sp_perMonth.tif'))

sp_sd<-app(sp_perMonth, sd, na.rm = T)
sp_mean<-app(sp_perMonth, mean, na.rm = T)
sp_meanthan1<-ifel(sp_mean>0.75, sp_mean, NA)

sp_cv <- sp_sd/sp_meanthan1*100
plot(sp_cv)
#writeRaster(sp_cv,paste0(basePath,'AE_data/sp_cv.tif'))

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



sp_perMonth <- rast(paste0(basePath,'AE_data/sp_perMonth.tif')) #%>% as.data.table
sp_df <- sp_perMonth[[c('X1','X2','X3','X4','X5','X6','X7','X8','X9','X10','X11','X12')]]
names(sp_df) <- c('Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec')

sp_df <- terra::as.data.frame(sp_df, xy = TRUE, na.rm = FALSE) 
head(sp_df)

# 使用 ggplot2 绘制折线图
library(tidyr)
sp_df_long <- gather(sp_df, key = "Month", value = "Value", -x, -y,na.rm = T)
names(sp_df_long)<-c("lon","lat","Month","Value")
head(sp_df_long)
sp_df_long$Month <- factor(sp_df_long$Month, levels = c("Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"))
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
  coord_fixed(ratio = 1.7)+ #调整y轴/x轴的比例
  # theme_classic()+
  theme(
    legend.position = "none",
    axis.title.y = element_blank(),
    axis.text.y=element_blank(),
    axis.title.x = element_text(size=16),
    axis.text.x=element_text(size=12),
    panel.grid = element_blank(),
    axis.ticks.length = unit(-4, "pt") #坐标轴ticks内向
  )

pSRCV+pSRCV2

#lit(Species Richness)----
#Species Richness
spNum <- lapply(speciesPixelNumPath2, function(x){
  r <- rast(x) %>% sum(.,na.rm=T)
  r2 <- ifel(r>0,1,NA)
  names(r2) <- str_sub(basename(x),1,-5)
  return(r2)
})
spNum <- rast(spNum)
spNumTif <- sum(spNum,na.rm=T) %>% mask(globalCountry)
#terra::writeRaster(spNumTif,paste0(basePath,'AE_data/spNumTif.tif'))


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
    legend.position = 'none',#设置图例与主图距离
    legend.direction='horizontal',#图例水平放置vertical，垂直horizontal，
    legend.key.width = unit(1.15,'cm'), #图例宽度
    legend.key.height = unit(0.3,'cm')
  )
pspNum

#
Vborder <- vect('/root/autodl-tmp/Wallace_zoogeographic/newValisBorder.shp') 
vRaster <- rasterize(Vborder,spNumTif,field='name')
valisEntropy <- c(spNumTif,vRaster) %>% terra::as.data.frame() %>%na.omit()


###plot
valisEntropy <- valisEntropy[valisEntropy$name != "Antarctica",]   #去掉不需要的南极界
Vborder1<-factor(valisEntropy$name,levels = c("Australian","Ethiopian","Nearctic","Neotropical","Oriental","Palearctic")) #自定义展示顺序

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
AE<-rast(paste0(basePath,'AE_data/AE.tif')); plot(AE)

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
    # panel.background = element_rect(fill = "white"),#背景设置
    #legend.direction = "horizontal",
    legend.title = element_blank(),
    #legend.title = element_text(hjust=0.5),
    legend.title.align = -10,
    legend.position = c(0.15, 0.05),#设置图例与主图距离
    legend.direction='horizontal',#图例水平放置vertical，垂直horizontal，
    legend.key.width = unit(1.15,'cm'), #图例宽度
    legend.key.height = unit(0.3,'cm')
  )
pAE


library(ggpointdensity)
AE_df <- terra::as.data.frame(AE, xy = TRUE, na.rm = T) 
head(AE_df)
names(AE_df)<-c("lon","lat","sum")

pAE2<-ggplot(AE_df,aes(x=sum, y=lat),)+
  geom_pointdensity(adjust = 4, size = 0.5)+                              # adjust：设置neighbors范围
  scale_color_distiller(palette = "Spectral", direction = -1)+ # 设置连续型颜色
  xlab('WAE')+
  xlim(c(0,5.5))+
  ylim(c(-56,90))+
  theme_bw()+
  # labs(y="Latitude(°)",)+
  #geom_ribbon(aes(y=lat, xmin=row_means-row_sds,xmax = row_means+row_sds),fill = "lightgrey", alpha=0.5)+
  coord_fixed(ratio = 0.105)+ #调整y轴/x轴的比例
  # theme_classic()+
  theme(
    axis.title.x = element_text(size=16),
    axis.text.x=element_text(size=12),
    axis.title.y = element_blank(),
    axis.text.y=element_blank(),
    panel.grid = element_blank(),
    axis.ticks.length = unit(-4, "pt"), #坐标轴ticks内向
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
valisEntropy <- valisEntropy[valisEntropy$name != "Antarctica",]   #去掉不需要的南极界
#Vborder1<-factor(valisEntropy$name,levels = c("Australian","Ethiopian","Nearctic","Neotropical","Oriental","Palearctic")) #自定义展示顺序
Vborder1<-factor(valisEntropy$name,levels = c("Oriental","Ethiopian","Australian","Neotropical","Nearctic","Palearctic")) #自定义展示顺序

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
  labs(y = "AE", x = " ")










#——————————————————————————————————————————————————————————-----------------

#######Fig2. AE validation#########################
#a##########
library(terra)
library(data.table)
library(dplyr)
library(stringr)
library(tidyterra)
# outBreak <- fread('/root/autodl-tmp/YANZHENG/point/allData.csv')
# outBreak$Longitude <- as.numeric(outBreak$Longitude)
# addGeom <- cellFromXY(globalRaster,outBreak[,c('Longitude','Latitude')]) %>% 
#   cbind(id=.,outBreak)
# thinData <- unique(data.table(addGeom),by='id') %>% dplyr::select(.,-id)
# outBreak2 <- vect(thinData,geom=c('Longitude','Latitude'),crs=crs)
# plot(outBreak2)
overEntropy <- rast(paste0(basePath,'AE_data/AE.tif'))%>% mask(globalCountry)
global(overEntropy,quantile,probs=seq(0, 1, 0.05),na.rm=T)


`%notin%` <- Negate(`%in%`)
library(lubridate)
outBreak1 <- fread('/root/autodl-tmp/YANZHENG/point/allData.csv') %>% subset(Diagnosis.status=='Confirmed'&Animal.type%in%c('Domestic','Wild'))

outBreak1$label <- str_extract(outBreak1$Serotype,'HPAI|LPAI')
outBreak1$h_label <- str_extract(outBreak1$Serotype,'H[0-9]N[0-9]|H[0-9]')
outBreak1 <- subset(outBreak1,outBreak1$h_label%notin%c('H9N2','H5N6'))   #这一句是用来运行全部的
#outBreak1 <- subset(outBreak1,outBreak1$h_label=='H5N1'&outBreak1$label=='HPAI')  #这一句是用来仅运行H5N1-HPAI的

outBreak1$Longitude <- as.numeric(outBreak1$Longitude)

addGeom <- cellFromXY(globalRaster,outBreak1[,c('Longitude','Latitude')]) %>% 
  cbind(id=.,outBreak1)
thinData <- unique(data.table(addGeom),by='id') %>% dplyr::select(.,-id)
outBreak2 <- vect(thinData,geom=c('Longitude','Latitude'),crs=crs)


plotEntropy <- ifel(overEntropy>=4.19,1,0)   #4.107
#加载海岸线数据
coast <- rnaturalearth::ne_coastline(scale = "small", returnclass = "sf")
crs <- '+proj=longlat +datum=WGS84'
ggplot() +
  geom_spatraster(data = plotEntropy) +
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
    # panel.background = element_rect(fill = "white"),#背景设置
    legend.title = element_blank(),
    legend.text = element_text(size = 15),
    legend.position = "none",
    #legend.position = c(0.15, 0.15),#设置图例与主图距离
    legend.direction='horizontal',#图例水平放置
    legend.key.width = unit(0.4,'cm'), #图例宽度
    legend.key.height = unit(3,'cm')
  )


#统计high-risk数量：
popd2015<-rast("/root/autodl-tmp/全球人口/GWP_v4/popd2015_30.tif") %>% resample(globalRaster)
poul2015<-rast("/root/autodl-tmp/zyresult/Poultry_duckchic.tif")%>% resample(globalRaster)
plot(plotEntropy);global(plotEntropy, sum, na.rm=T)

plotEntropy_df<-c(plotEntropy,popd2015,poul2015) %>% terra::as.data.frame() %>%na.omit()
head(plotEntropy_df)
names(plotEntropy_df)<-c("plotEntropy","Numpop","Numpoul")

plotEntropy_result <- plotEntropy_df %>%
  group_by(plotEntropy) %>%
  summarise(across(everything(), ~ sum(.x, na.rm = TRUE)))
plotEntropy_result
#fwrite(NumhotOR_result,'/root/autodl-tmp/humPoulResult/data/Hot_data/hotOR_Num_result.csv')

# 统计值为1的像元数量
ones <- sum(values(plotEntropy) == 1, na.rm = TRUE); ones
# 统计值为0的像元数量
zeros <- sum(values(plotEntropy) == 0, na.rm = TRUE); zeros


#a little##########
`%notin%` <- Negate(`%in%`)
library(lubridate)
#outBreak <- fread('/root/autodl-tmp/YANZHENG/point/allData.csv')
outBreak1 <- fread('/root/autodl-tmp/YANZHENG/point/allData.csv') %>% subset(Diagnosis.status=='Confirmed'&Animal.type%in%c('Domestic','Wild'))
# outBreak <- subset(outBreak,outBreak$Diagnosis.status=='Confirmed')
# r <- subset(outBreak,outBreak$Species=='Unspecified Bird'&outBreak$Animal.type=='	
# Domestic')
outBreak1$label <- str_extract(outBreak1$Serotype,'HPAI|LPAI')
outBreak1$h_label <- str_extract(outBreak1$Serotype,'H[0-9]N[0-9]|H[0-9]')
# r <- outBreak[outBreak$Longitude>60&outBreak$Longitude<120&outBreak$Latitude>35&outBreak$Latitude<60,]
#outBreak1 <- subset(outBreak1,outBreak1$h_label%in%c('H9N2','H5N6'))
outBreak1 <- subset(outBreak1,outBreak1$h_label%notin%c('H9N2','H5N6'))
#outBreak1 <- subset(outBreak1,outBreak1$h_label=='H5N1'&outBreak1$label=='HPAI')


outBreak1$Longitude <- as.numeric(outBreak1$Longitude)

# outBreak <- subset(outBreak,outBreak$year%in%2021:2023)
addGeom <- cellFromXY(globalRaster,outBreak1[,c('Longitude','Latitude')]) %>% 
  cbind(id=.,outBreak1)
thinData <- unique(data.table(addGeom),by='id') %>% dplyr::select(.,-id)
outBreak2 <- vect(thinData,geom=c('Longitude','Latitude'),crs=crs)
# r <- terra::extract(v,outBreak2)

overEntropy <- rast(paste0(basePath,'AE_data/AE.tif'))%>% crop(globalCountry,mask=T)
# global(overEntropy,quantile,probs=seq(0, 1, 0.05),na.rm=T)

# overEntropy2 <- ifel(overEntropy>4.152,overEntropy,NA)


df <- terra::extract(overEntropy,outBreak2,mean,na.rm=T,cells=T)%>% as.data.frame()
allDf <- terra::as.data.frame(overEntropy,na.rm=T,cells=T)

library(pROC)

df$label <- 'good'
rocDf <- left_join(allDf,df[,c('cell','label')])
rocDf$label <- ifelse(is.na(rocDf$label),'poor',rocDf$label)
dfroc1<- roc(rocDf$label, rocDf$sum)
plot(dfroc1,col="red",#颜色
     legacy.axes=T,#y轴格式更改
     print.auc=TRUE,#显示AUC面积
     expand=c(0,0),
     print.thres=TRUE,#添加截点和95%CI
     grid=c(0.2,0.2),grid.col=c("lightgrey","lightgrey"),
     #cex.text=1.5,
     cex.main=1,  # 主标题字体放大1.5倍
     cex.sub=1,   # 副标题字体放大1.5倍
     cex.axis=1,  # 坐标轴刻度字体放大1.5倍
     cex.lab=1)#网格线设置



threshold <- 4.107   #4.107

TP <- sum(df$sum>=threshold,na.rm = T)  ;TP
FN <- sum(df$sum<threshold,na.rm = T)   ;FN
FP <- sum(allDf$sum>=threshold,na.rm = T)-TP  ;FP
TN <- nrow(allDf)-TP-FN-FP  ;TN

accuracy <- (TP+TN)/(TP+TN+FN+FP)  ;accuracy
precise <- (TP)/(TP+FP)  ;precise
recall <- (TP)/(TP+FN)  ;recall
FI <- 2*((recall*precise)/(recall+precise))  ;FI




#b、c###########
#新accuracy
result_country <- fread('/root/autodl-tmp/humPoulResult/data/result_country.csv')
healthData <- fread('/root/autodl-tmp/humPoulResult/data/health expenditure.csv',header = T,drop=c(5:44,66:69))
globalSHP <- vect('/root/autodl-tmp/worldBorder/ne_10m_admin_0_countries.shp')
countryData <- as.data.frame(globalSHP) %>% .[,c('NAME_LONG','POP_EST','POP_RANK','ECONOMY','INCOME_GRP','CONTINENT','GDP_MD')]
result_df3 <- left_join(result_country,countryData,by=c('country'='NAME_LONG'))
result_df4 <- left_join(result_df3,healthData,by=c('country'='Country Name')) %>% na.omit()
result_df4$health <- rowMeans(result_df4[,21:41])
result_df4$per <- (result_df4$df/result_df4$alldf)

# result_df3$tt <- result_df3$TP+result_df3$FN
# result_df4 <- subset(result_df3,result_df3$tt>10)

#比较GDP和health哪个好：GDP更好
# p1<-ggplot(data=result_df4)+
#   geom_point(aes(log(GDP_MD),accuracy,color=CONTINENT,size=POP_EST,alpha=0.5))+
#   geom_smooth(aes(log(GDP_MD),accuracy),method = "lm", formula = y ~ I(x^-1))+
#   xlab('GDP_MD (log)')+
#   ylab('Accuracy')+
#   theme_bw()+
#   theme(
#     panel.background = element_blank()
#   )
# 
# p2<-ggplot(data=result_df4)+
#   geom_point(aes(log(health),accuracy,color=CONTINENT,size=POP_EST,alpha=0.5))+
#   geom_smooth(aes(log(health),accuracy),method = "lm", formula = y ~ I(x^-1))+
#   xlab('health (log)')+
#   ylab('Accuracy')+
#   theme_bw()+
#   theme(
#     panel.background = element_blank()
#   )
# 
# library(patchwork)
# p1+p2


#用GDP
# 计算POP_EST的十个分位数
pop_quantiles <- quantile(result_df4$POP_EST, probs=seq(0, 1, by=0.1), na.rm = TRUE)

# 将POP_EST映射到1到10的大小
result_df4$POP_EST_size <- findInterval(result_df4$POP_EST, vec = pop_quantiles)
result_df4$POP_EST <- as.numeric(result_df4$POP_EST)
result_df4 <- mutate(result_df4, POP_EST_label = cut(POP_EST, breaks = c(0, 10000000, 20000000, 50000000, 100000000, 500000000, 1000000000, 1397800000),
                                               labels = c("<1000", "1000 - 2000", "2000-5000", "5000 - 10000", "10000 - 50000", "50000 - 100000", ">100000")))

# 按大洲和人口排序
result_df4 <- result_df4[order(result_df4$CONTINENT, -result_df4$POP_EST), ]

#选择每个大洲人口排名前三的国家
# top_countries_by_continent <- result_df4 %>%
#   group_by(CONTINENT) %>%
#   slice_max(POP_EST, n = 3) %>%
#   ungroup() %>%
#   dplyr::select(country)
# top_countries_by_continent<-c(top_countries_by_continent)
# 
# result_dfc <- result_df4[result_df4$country %in% top_countries_by_continent$country, ]   #[-11]为了去掉Canada
country5 <- c("China", "United States", "India", "Nigeria")
result_dfc <- result_df4[result_df4$country %in% country5]
result_dfc <- result_dfc %>%
  mutate(country = ifelse(country == "United States", "USA", country))
# # 使用 lm() 函数进行线性拟合
# linear_model <- lm(accuracy ~ log(POP_EST), data = result_df4)
# summary(linear_model)
# 绘制图形
#size_labels <- c("<1000", "1000 - 2000", "2000-5000", "5000 - 10000", "10000 - 50000", "50000 - 100000", ">100000")
size_values <- c("<1000" = 3, "1000 - 2000" = 4, "2000-5000" = 6, "5000 - 10000" = 8, "10000 - 50000" = 10, "50000 - 100000" = 12, ">100000" = 18)

my_formula <- y ~ x

#GDP
p2<-ggplot(result_df4, aes(log(GDP_MD), accuracy, color=CONTINENT)) +
  geom_point(aes(size=POP_EST_label), alpha=0.5,show.legend = F) +
  geom_text(data=result_dfc, aes(label=country, x=log(GDP_MD), y=accuracy),show.legend = F,size=4.5, hjust=0.5, vjust=2) +
  scale_color_manual(values = c("#984EA3", "#E41A1C","#4DAF4A",  "#FF7F00", "#a38900","#377EB8"))+
  #geom_smooth(method = "gam", se = T, fill="lightgrey",color = "grey", alpha=0.3, linetype = "dashed") +  # 全部点的趋势线，黑色  #loess  gam
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

p2


#per!!
result_country <- fread('/root/autodl-tmp/humPoulResult/data/result_country.csv')
healthData <- fread('/root/autodl-tmp/humPoulResult/data/health expenditure.csv',header = T,drop=c(5:44,66:69))
globalSHP <- vect('/root/autodl-tmp/worldBorder/ne_10m_admin_0_countries.shp')
countryData <- as.data.frame(globalSHP) %>% .[,c('NAME_LONG','POP_EST','POP_RANK','ECONOMY','INCOME_GRP','CONTINENT','GDP_MD')]
result_df3 <- left_join(result_country,countryData,by=c('country'='NAME_LONG'))
result_df4 <- left_join(result_df3,healthData,by=c('country'='Country Name')) %>% na.omit()
result_df4$per <- (result_df4$df/result_df4$alldf)


# 计算POP_EST的十个分位数
pop_quantiles <- quantile(result_df4$POP_EST, probs=seq(0, 1, by=0.1), na.rm = TRUE)
# 将POP_EST映射到1到10的大小
result_df4$POP_EST_size <- findInterval(result_df4$POP_EST, vec = pop_quantiles)
result_df4$POP_EST <- as.numeric(result_df4$POP_EST)
result_df4 <- mutate(result_df4, POP_EST_label = cut(POP_EST, breaks = c(0, 10000000, 20000000, 50000000, 100000000, 500000000, 1000000000, 1397800000),
                                                     labels = c("<1000", "1000 - 2000", "2000-5000", "5000 - 10000", "10000 - 50000", "50000 - 100000", ">100000")))

# 按大洲和人口排序
result_df4 <- result_df4[order(result_df4$CONTINENT, -result_df4$POP_EST), ]

# # 选择每个大洲人口排名前三的国家
# top_countries_by_continent <- result_df4 %>%
#   group_by(CONTINENT) %>%
#   slice_max(POP_EST, n = 3) %>%
#   ungroup() %>%
#   dplyr::select(country)
# top_countries_by_continent<-c(top_countries_by_continent)
# result_dfc <- result_df4[result_df4$country %in% top_countries_by_continent$country, ]   #[-11]为了去掉Canada
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
p1<-ggplot(result_df5, aes(per, accuracy, color=CONTINENT)) +
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


p1+p2

#病例类型############


birdName <- c("White Stork:ciconia Ciconia(Ciconiidae)", "Saker Falcon (Falco Cherrug)", 
              "Helmeted Guineafowl (Numida Meleagris)", "Spur-Winged Goose", 
              "African Fish Eagle (Haliaeetus Vocifer)", "European Turtle Dove", 
              "Blue Crane (Grus Paradisea)", "Cape Gannet Shy Albatross", 
              "Sandwich Tern", "Hartlaub's Gull (Chroicocephalus Hartlaubii)", 
              "Crowned Cormorant (Microcarbo Coronatus)", "Laughing Dove", 
              "Mediterranean Gull", "Spotted Eagle-Owl (Bubo Africanus)", 
              "African Penguin (Spheniscus Demersus)", "Armenian Gull (Larus Armenicus)", 
              "Eurasian Blackbird (Turdus Merula)", "Silver Teal (Anas Versicolor)", 
              "Kelp Gull", "Arctic Skua (Or Parasitic Jaeger)", 
              "Wild Turkey (Meleagris Gallopavo)", "Black-Headed Heron (Ardea Melanocephala)", 
              "Grey Headed Gull (Chroicocephalus Cirrocephalus)", "Oyster Catcher", 
              "Common Moorhen (Gallinula Chloropus)", "Eurasian Oystercatcher", 
              "Cape Gannet (Morus Capensis)", "Strigidae (Unidentified):Strigidae (Incognita)(Strigidae)", 
              "Columbidae (Unidentified):Columbidae (Incognita)(Columbidae)", 
              "Eurasian Jackdaw (Coloeus Monedula)", "Arctic Tern (Sterna Paradisea)", 
              "Ivory Gull (Pagophila Eburnea)", "Purple Sandpiper (Calidris Maritima)", 
              "Eurasian Tree Sparrow", "Spotted Dove", "Ardeidae (Unspecified)", 
              "Wigeons:mareca Penelope(Anatidae)", "Sulids (Sulidae)", 
              "Slender-Billed Gull", "Glaucous Gull (Larus Hyperboreus)", 
              "Gruidae (Unidentified):Gruidae (Incognita)(Gruidae)", 
              "Northern Fulmar (Fulmarus Glacialis)", "Charadrius Pallidus", 
              "Tawny Owl (Strix Aluco)", "Pink-Footed Goose (Anser Brachyrhynchus)", 
              "Eurasian Collared Dove (Streptopelia Decaocto)", "Green Sandpiper", 
              "Haematopodidae (Unidentified):Haematopodidae (Incognita)(Haematopodidae)", 
              "Greater Rhea (Rhea Americana)", "Common Wood-Pigeon:columba Palumbus(Columbidae)", 
              "Passeridae (Unidentified):Passeridae (Incognita)(Passeridae)", 
              "Red Knot", "Semipalmated Sandpiper", "Thrush (Turdidae)", 
              "Sanderling", "Numididae (Unidentified)", "Great Blue Heron", 
              "Western Grebe (Aechmophorus Occidentalis):Podicipedidae-Procellariiformes", 
              "Common Snipe", "Black-Necked Grebe (Podiceps Nigricollis):Podicipedidae-Procellariiformes", 
              "Red-Billed Gull (Larus Scopulinus):Laridae-Charadriiformes", 
              "Red Junglefowl", "Little Owl (Athene Noctua):Strigidae-Suliformes", 
              "Griffon Vulture (Gyps Fulvus)", "Cackling Goose (Branta Hutchinsii)", 
              "Pied Avocet", "Black-Winged Stilt", "Ferruginous Duck (Aythya Nyroca)", 
              "Ruff", "Little Ringed Plover", "Kentish Plover", 
              "Common Ringed Plover", "Little Stint", "Wader", 
              "Gull-Billed Tern (Gelochelidon Nilotica)", "Black Guillemot (Cepphus Grylle):Alcidae-Charadriiformes", 
              "Bearded Vulture (Gypaetus Barbatus)", "Black Guillemot (Cepphus Grylle)", 
              "Atlantic Puffin (Fratercula Arctica)", "Northern Lapwing", 
              "Common Loon (Gavia Immer)", "Northern Crested Caracara (Caracara Cheriway)", 
              "Variable Hawk (Geranoaetus Polyosoma)", "Northern Bald Ibis (Geronticus Eremita)", 
              "Chilean Flamingo (Phoenicopterus Chilensis)", "Royal Tern", 
              "Arenaria Interpres", "Brown Skua (Stercorarius Antarcticus)", 
              "Cabot's Tern (Thalasseus Acuflavidus)", "South American Tern", 
              "Brown Booby (Sula Leucogaster)", "Neotropic Cormorant (Nannopterum Brasilianum)", 
              "Gray Hawk (Buteo Plagiatus)", "Magnificent Frigatebird (Fregata Magnificens)", 
              "Manx Shearwater (Puffinus Puffinus)", "Magellanic Penguin", 
              "Southern Lapwing", "White-Faced Ibis (Plegadis Chihi)",
              "Black-Necked Swan (Cygnus Melancoryphus)", "Iceland Gull (Larus Glaucoides)",
              "American Crow (Corvus Brachyrhynchos)", "Corvus Brachyrhynchos: (Corvidae) American Crow",
              "Turkey Vulture (Cathartes Aura)", "Adélie Penguin (Pygoscelis Adeliae)",
              "Gentoo Penguin (Pygoscelis Papua)", "Southern Fulmar (Fulmarus Glacialoides)",
              "Black-Browed Albatross (Thalassarche Melanophris)", "Penguin",
              "Red-Breasted Merganser (Mergus Serrator)", "Bald Eagle (Haliaeetus Leucocephalus)",
              "Red Breasted Goose", "Red-Tailed Hawk", "Andean Goose (Chloephaga Melanoptera)",
              "Ring-Billed Gull (Larus Delawarensis)", "Double-Crested Cormorant (Nannopterum Auritum)",
              "Great Horned Owl (Bubo Virginianus)", "Imperial Shag", "Grey Gull(Leucophaeus Modestus)",
              "Guanay Cormorant (Leucocarbo Bougainvilliorum)", "Peruvian Pelican (Pelecanus Thagus)",
              "Humboldt Penguin (Spheniscus Humboldti)", "Peruvian Booby (Sula Variegata)",
              "American Oystercatcher", "Black Skimmer", "Red-Gartered Coot (Fulica Armillata)",
              "Guanay Cormorant Or Guanay Shag (Leucocarbo Bougainvilliorum)", "Snow Goose (Anser Caerulescens)",
              "Rough-Legged Hawk", "Surf Scoter (Melanitta Perspicillata):Anatidae-Anseriformes",
              "Wood Duck (Aix Sponsa) (Anatidae)", "Red Crested Pochard (Netta Rufina)",
              "American Wigeon", "Brown Pelican (Pelecanus Occidentalis)", "American Green-Winged Teal",
              "Trumpeter Swan (Cygnus Buccinator)", "Blue-Footed Booby (Sula Nebouxii)",
              "Green-Winged Teal/Eurasian Teal", "Thick-Billed Murre (Uria Lomvia):Alcidae-Charadriiformes",
              "Cooper's Hawk (Accipiter Cooperii)", "Blue Winged Teal (Anas Discors):Anatidae-Anseriformes",
              "Swainson's Hawk (Buteo Swainsoni)", "Black-Billed Magpie (Pica Hudsonia):Corvidae-Passeriformes",
              "American White Pelican (Pelecanus Erythrorhynchos):Pelecanidae-Phaethontiformes",
              "Great-Tailed Grackle (Quiscalus Mexicanus):Icteridae-Pelecaniformes", "Common Grackle",
              "Barn Swallow", "Franklin's Gull", "Greater Yellowlegs", "Eared Grebe",
              "Barred Owl (Strix Varia)", "Ross's Goose (Anser Rossii)", "Cinnamon Teal",
              "Red-Footed Booby (Sula Sula)", "Lesser Snow Goose", "Ring Necked Duck",
              "Greater Frigate Bird", "Red-Gartered Coot",
              "Helmeted Guineafowl (Numida Meleagris)", "Wild Turkey (Meleagris Gallopavo)", "Red Junglefowl")


poultryName <- c()


mammalName <- c("Stone Marten", "Common Raccoon (Procyon Lotor)", 
                "South American Coati (Nasua Nasua):Procyonidae-Carnivora", 
                "Gray Seal (Halichoerus Grypus):Phocidae-Carnivora", 
                "Polecat", "Ferret", "South America Fur Seal (Arctophoca Australis)", 
                "Southern Elephant Seal (M. Leonina)", "South American Sea Lion (Otaria Flavescens)",
                "Coati (Gen. Nasua)", "Skunk (Mephitis Mephitis)",
                "Marine Otter (Lontra Felina)", "Sea Otter (Enhydra Lutris)",
                "Burmeister's Porpoise (Phocoena Spinipinnis)", "Sea Lion",
                "Southern River Otter (Lontra Provocax)", "Amur Leopard (Panthera Pardus Orientalis)",
                "Wild Fox", "Bottlenose Dolphin (Tursiops Truncatus)",
                "Virginia Opossum (Didelphis Virginiana)", "Cattle",
                "Coyote (Canis Latrans)", "Bobcat (Lynx Rufus)",
                "Fisher (Pekania Pennanti)", "American Black Bear (Ursus Americanus)",
                "Cougar (Puma Concolor)", "Kodiak Grizzly Bear (Ursus Arctos Horribilis)",
                "Wild Cat", "Rodents", "Polar Bear (Ursus Maritimus)",
                "Chilean Dolphin (Cephalorhynchus Eutropia)")


EnvironmentalName <- c() # 没有环境名称







#五个国家多少病例数？#######
`%notin%` <- Negate(`%in%`)
library(lubridate)
outBreak1 <- fread('/root/autodl-tmp/YANZHENG/point/allData.csv') %>% subset(Diagnosis.status=='Confirmed'&Animal.type%in%c('Domestic','Wild'))

#outBreak1$label <- str_extract(outBreak1$Serotype,'HPAI|LPAI')
#outBreak1$h_label <- str_extract(outBreak1$Serotype,'H[0-9]N[0-9]|H[0-9]')
#outBreak1 <- subset(outBreak1,outBreak1$h_label%notin%c('H9N2','H5N6'))   #这一句是用来运行全部的
#outBreak1 <- subset(outBreak1,outBreak1$h_label=='H5N1'&outBreak1$label=='HPAI')  #这一句是用来仅运行H5N1-HPAI的

outBreak1$Longitude <- as.numeric(outBreak1$Longitude)

addGeom <- cellFromXY(globalRaster,outBreak1[,c('Longitude','Latitude')]) %>% 
  cbind(id=.,outBreak1)
thinData <- unique(data.table(addGeom),by='id') %>% dplyr::select(.,-id)
outBreak2 <- vect(thinData,geom=c('Longitude','Latitude'),crs=crs)

plot(outBreak2)


globalSHP <- vect('/root/autodl-tmp/zyresult/Con_popentrpoul_sf_EU.shp')
globalSHP2 <- terra::aggregate(globalSHP,'name_ec')
Contry5c <-c('China','India','European Country','Nigeria','United States')
country <- subset(globalSHP2, globalSHP2$name_ec %in% Contry5c ) 

# 
# points_in_country <- extract(country, outBreak2)
# num_points_in_country <- sum(!is.na(points_in_country[,1])) ; num_points_in_country
# total_points <- nrow(outBreak2) ; total_points
# 
# proportion <- num_points_in_country / total_points  ; proportion

# 转换为sf对象
outBreak2 <- as(sf::st_as_sf(outBreak2), "sf")
country <- as(sf::st_as_sf(country), "sf")

outBreak2 <- st_transform(outBreak2, st_crs(country))
points_in_country <- st_intersects(outBreak2, country, sparse = FALSE)

num_points_in_country <- sum(points_in_country) ; num_points_in_country
total_points <- nrow(outBreak2) ; total_points
proportion <- num_points_in_country / total_points ; proportion



#——————————————————————————————————————————————————————————-----------------

#######Fig3. Hotspots#########################
#a########
library(devtools)
library(tricolore)

popd2015<-rast("/root/autodl-tmp/全球人口/GWP_v4/popd2015_30.tif") %>% resample(globalRaster)%>% mask(globalCountry)
Entropy<-rast(paste0(basePath,'AE_data/AE.tif'))%>% mask(globalCountry)
poul2015<-rast("/root/autodl-tmp/zyresult/Poultry_duckchic.tif")%>% resample(globalRaster)%>% mask(globalCountry)

# 计算分位数
quantiles_pop <- global(popd2015,quantile,probs=seq(0, 1, 0.1),na.rm=T)
quantiles_entr<- global(Entropy,quantile,probs=seq(0, 1, 0.1),na.rm=T)
quantiles_poul<- global(poul2015,quantile,probs=seq(0, 1, 0.1),na.rm=T)

#重分类
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



# 查看结果
head(popentrpoul_df2)

# 使用 Tricolore 创建配色方案
library(tricolore)  #加载后ggplot会出错

# col <- Tricolore(popentrpoul_df2, p1 = 'pop_re', p2 = 'poul_re', p3 = 'entr_re',
#                  hue =  0.7, # 品红、黄色、青色的大致色调值
#                  chroma = 0.9,            # 颜色的饱和度
#                  lightness = 0.7)         # 颜色的亮度
col <- Tricolore(popentrpoul_df2, p1 = 'pop_re', p2 = 'poul_re', p3 = 'entr_re', breaks = Inf)  
col$key
popentrpoul_df2$rgb <- col$rgb


#出图
plot_res <- ggplot(popentrpoul_df2) +
  geom_tile(aes(x = x, y = y, fill = rgb)) +
  scale_fill_identity() +
  geom_spatvector(data=coast,fill=NA)+coord_sf(crs = crs,xlim=c(-160,165),ylim=c(-56,90))+
  labs(x = NULL, y = NULL,) +
  #theme_minimal()+
  theme_bw()+
  theme(
    legend.position="none",
        axis.text = element_text(size=12),
    panel.grid.major = element_blank(),  # 删除主要的灰色虚线
  )
print(plot_res)

library(ggtern)
plot_res+
  annotation_custom(ggplotGrob(col$key),
                    xmin = -170, xmax =-100, ymin = -60, ymax = 10)



#b############

popd2015<-rast("/root/autodl-tmp/全球人口/GWP_v4/popd2015_30.tif") %>% resample(globalRaster)
Entropy<-rast(paste0(basePath,'AE_data/AE.tif'))%>% mask(globalCountry)
poul2015<-rast("/root/autodl-tmp/zyresult/Poultry_duckchic.tif")%>% resample(globalRaster)


#####国家矢量
globalSHP <- vect('/root/autodl-tmp/zyresult/Con_popentrpoul_sf_EU.shp')
globalSHP2 <- terra::aggregate(globalSHP,'name_ec')


###取阈值
quantiles_pop <- global(popd2015,quantile,probs=seq(0, 1, 0.01),na.rm=T)
quantiles_entr<- global(Entropy,quantile,probs=seq(0, 1, 0.01),na.rm=T)
quantiles_poul<- global(poul2015,quantile,probs=seq(0, 1, 0.01),na.rm=T)

allData <- data.table()
for (i in 51:100) {           #6:10
  #取阈值
  quan_pop <- quantiles_pop[1,i]
  quan_entr<- quantiles_entr[1,i]
  quan_poul <- quantiles_poul[1,i]
  #重分类
  pop_new <- ifel(popd2015>quan_pop,1,0)
  entr_new <- ifel(Entropy>quan_entr,1,0)
  poul_new <- ifel(poul2015>quan_poul,1,0)
  
  hotArea <- pop_new+entr_new+poul_new
  hotArea2 <- ifel(hotArea==3,1,NA)
  #区域统计
  hotResult <- zonal(hotArea,z=globalSHP2,'sum',na.rm=T)
  hotResult2 <- data.table(countryName=globalSHP2$name_ec,hotSum=hotResult$popd2015_30,threshold=i-1)   #threshold=(i-1)*10
  allData <- rbind(allData,hotResult2)
  # #计算面积
  # crs(hotArea2) <- "+proj=aea +lat_1=first_standard_parallel +lat_2=second_standard_parallel +lat_0=latitude_of_origin +lon_0=central_meridian +x_0=false_easting +y_0=false_northing +datum=WGS84"   #Albers Conic Equal Area投影
  # crs(globalSHP2) <- crs(hotArea2)
  # hot_zonal <- zonal(hotArea2 == 1, globalSHP2, fun='sum', na.rm=TRUE)
  # area_per_cell <- (111.32 * cos(lat_from_y(hotArea2_raster))) * 111.32 * (res(hotArea2_raster)[1] * res(hotArea2_raster)[2])
  # 
}


#######(-)重新画一下entr约登指数但50人分位数-----------
quan_pop <- 50
quan_entr<- 4.107 #4.107
quan_poul <- quantiles_poul[1,84]
#重分类
pop_new <- ifel(popd2015>quan_pop,1,0)
entr_new <- ifel(Entropy>quan_entr,10,0)
poul_new <- ifel(poul2015>quan_poul,2,0)

hotentrpoppoul <- pop_new+entr_new+poul_new   ;plot(hotentrpoppoul)
hotAND <- ifel(hotentrpoppoul==13,1,0)  ; plot(hotAND); global(hotAND,sum,na.rm=T)
hotOR_pop <- ifel(hotentrpoppoul==11|hotentrpoppoul==13,1,0) ; plot(hotOR_pop)
hotOR_poul <- ifel(hotentrpoppoul==12|hotentrpoppoul==13,1,0)  ; plot(hotOR_poul)
Nonehot<- ifel(hotentrpoppoul==0,1,0)  ; plot(Nonehot); global(Nonehot,sum,na.rm=T)

# writeRaster(hotAND,'/root/autodl-tmp/humPoulResult/data/Hot_data/hotAND.tif')
# writeRaster(hotOR_pop,'/root/autodl-tmp/humPoulResult/data/Hot_data/hotOR_pop.tif')
# writeRaster(hotOR_poul,'/root/autodl-tmp/humPoulResult/data/Hot_data/hotOR_poul.tif')
# writeRaster(Nonehot,'/root/autodl-tmp/humPoulResult/data/Hot_data/Nonehot.tif')
#writeRaster(hotentrpoppoul,'/root/autodl-tmp/humPoulResult/data/Hot_data/hotentrpoppoul.tif')
library(rnaturalearth)
coast <- ne_coastline(scale = "small", returnclass = "sf") %>% vect()
crs <- '+proj=longlat +datum=WGS84'
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
    legend.position="none",
    axis.text = element_text(size=12),
    panel.grid.major = element_blank(),  # 删除主要的灰色虚线
    #panel.grid.minor = element_blank(),   # 删除次要的灰色虚线
    #panel.background = element_rect(fill="white")#背景设置
  )


#主图3红橙黄统计   求出暴露的人口总数和家禽总数--------------
popd2015<-rast("/root/autodl-tmp/全球人口/GWP_v4/popd2015_30.tif") %>% resample(globalRaster)
#Entropy<-rast(paste0(basePath,'AE_data/AE.tif'))%>% mask(globalCountry)
poul2015<-rast("/root/autodl-tmp/zyresult/Poultry_duckchic.tif")%>% resample(globalRaster)

#像元数
hotAND <- ifel(hotentrpoppoul==13,1,0)  ; plot(hotAND); global(hotAND,sum,na.rm=T)
hotAND_pop <- ifel(hotentrpoppoul==11,1,0) ; plot(hotAND_pop); global(hotAND_pop,sum,na.rm=T)
hotAND_poul <- ifel(hotentrpoppoul==12,1,0)  ; plot(hotAND_poul); global(hotAND_poul,sum,na.rm=T)
Nonehot<- ifel(hotentrpoppoul==0,1,0)  ; plot(Nonehot); global(Nonehot,sum,na.rm=T)

Numpop_hotAND<-hotAND*popd2015; plot(Numpop_hotAND); global(Numpop_hotAND,sum,na.rm=T)
Numpoul_hotAND<-hotAND*poul2015; plot(Numpoul_hotAND); global(Numpoul_hotAND,sum,na.rm=T)

Numpop_hotAND_pop<-hotAND_pop*popd2015; plot(Numpop_hotAND_pop); global(Numpop_hotAND_pop,sum,na.rm=T)
Numpoul_hotAND_pop<-hotAND_pop*poul2015; plot(Numpoul_hotAND_pop); global(Numpoul_hotAND_pop,sum,na.rm=T)

Numpop_hotAND_poul<-hotAND_poul*popd2015; plot(Numpop_hotAND_poul); global(Numpop_hotAND_poul,sum,na.rm=T)
Numpoul_hotAND_poul<-hotAND_poul*poul2015; plot(Numpoul_hotAND_poul); global(Numpoul_hotAND_poul,sum,na.rm=T)

#统计总数
plot(hotentrpoppoul)

hotentrpoppoul_df<-c(hotentrpoppoul,popd2015,poul2015) %>% terra::as.data.frame() %>%na.omit()
head(hotentrpoppoul_df)
names(hotentrpoppoul_df)<-c("type","pop","poul")

Num_result <- hotentrpoppoul_df %>%
  group_by(type) %>%
  summarise(across(everything(), ~ sum(.x, na.rm = TRUE)))
Num_result




#各大洲求出暴露的人口总数和家禽总数--------------
popd2015<-rast("/root/autodl-tmp/全球人口/GWP_v4/popd2015_30.tif") %>% resample(globalRaster)
Entropy<-rast(paste0(basePath,'AE_data/AE.tif'))%>% mask(globalCountry)
poul2015<-rast("/root/autodl-tmp/zyresult/Poultry_duckchic.tif")%>% resample(globalRaster)

Continent_vect<- vect ("/root/autodl-tmp/worldBorder/continentNew.shp") |> dplyr::filter(CONTINENT != "Antarctica")
ContinentRaster <- rasterize(Continent_vect,Entropy,field='CONTINENT')

Numpop_hotAND<-hotAND*popd2015; plot(Numpop_hotAND)
Numpoul_hotAND<-hotAND*poul2015; plot(Numpoul_hotAND)

AllhotAND_df <-c(Entropy,Numpop_hotAND,Numpoul_hotAND,ContinentRaster) %>% terra::as.data.frame() %>%na.omit()
head(AllhotAND_df)
names(AllhotAND_df)<-c("entr","Numpop_hotAND","Numpoul_hotAND","continent")

NumhotAND_result <- AllhotAND_df %>%
  group_by(continent) %>%
  summarise(across(everything(), ~ sum(.x, na.rm = TRUE)))
NumhotAND_result
#fwrite(Num_result,'/root/autodl-tmp/humPoulResult/data/Hot_data/hot_Num_result.csv')



#饼图
NumhotAND_result <- NumhotAND_result %>%
  mutate(
    Percent_Numpop = round((Numpop_hotAND / sum(Numpop_hotAND)) * 100, 2),
    Percent_Numpoul = round((Numpoul_hotAND / sum(Numpoul_hotAND)) * 100, 2)
  )
NumhotAND_result$continent<-factor(NumhotAND_result$continent,levels = c("Asia","Africa","Europe","North America","South America","Oceania"))
colors1 <- c( "#df7415", "#f57f17", "#f58a2c", "#f69641","#f7a156", "#f8ad6b")
colors2 <- c( "#a38900", "#ba9c00", "#d1b000", "#e8c300","#ffd600", "#ffe45c")

# Numpop_hotAND 的饼图
pNumpop_hotAND<-ggplot(NumhotAND_result, aes(x="", y=Numpop_hotAND, fill=continent)) +
  geom_bar(width = 1, stat = "identity") +
  coord_polar("y", start=0) +
  scale_fill_manual(values=colors1) +
  #geom_text(aes(label = paste0(round(Percent_Numpop, 1), "%")), position = position_stack(vjust = 0.5)) +
  theme_void() +
  labs(fill="Continent")
pNumpop_hotAND
# Numpoul_hotAND 的饼图
pNumpoul_hotAND<-ggplot(NumhotAND_result, aes(x="", y=Numpoul_hotAND, fill=continent)) +
  geom_bar(width = 1, stat = "identity") +
  coord_polar("y", start=0) +
  scale_fill_manual(values=colors2) +
  #geom_text(aes(label = paste0(round(Percent_Numpoul, 1), "%")), position = position_stack(vjust = 0.5)) +
  theme_void() +
  labs(fill="Continent")
pNumpoul_hotAND


pNumpop_hotAND+pNumpoul_hotAND


#各国家求出暴露的人口总数和家禽总数------------------
popd2015<-rast("/root/autodl-tmp/全球人口/GWP_v4/popd2015_30.tif") %>% resample(globalRaster)
Entropy<-rast(paste0(basePath,'AE_data/AE.tif'))%>% mask(globalCountry)
poul2015<-rast("/root/autodl-tmp/zyresult/Poultry_duckchic.tif")%>% resample(globalRaster)

globalSHP <- vect('/root/autodl-tmp/zyresult/Con_popentrpoul_sf_EU.shp')
globalSHP2 <- terra::aggregate(globalSHP,'name_ec')
# countries <- c("China", "India", "European Country", "Nigeria", "United States",
#                "Thailand", "Indonesia", "Ethiopia", "Russian Federation", "Myanmar",
#                "Japan", "Pakistan", "Vietnam", "Mexico", "Brazil",
#                "Uganda", "Ukraine", "Burkina Faso", "United Kingdom", "Philippines")
# selectedCountries <- globalSHP2[globalSHP2$name_ec %in% countries,]
# countryRaster <- rasterize(selectedCountries, Entropy, field='name_ec')
countryRaster <- rasterize(globalSHP2, Entropy, field='name_ec')


Numpop_hotAND<-hotAND*popd2015; plot(Numpop_hotAND)
Numpoul_hotAND<-hotAND*poul2015; plot(Numpoul_hotAND)

AllhotANDcountry_df <-c(Entropy,Numpop_hotAND,Numpoul_hotAND,countryRaster) %>% terra::as.data.frame() %>%na.omit()
head(AllhotANDcountry_df)
names(AllhotANDcountry_df)<-c("entr","Numpop_hotAND","Numpoul_hotAND","country")

NumhotANDcountry_result <- AllhotANDcountry_df %>%
  group_by(country) %>%
  summarise(across(everything(), ~ sum(.x, na.rm = TRUE)))
NumhotANDcountry_result
#fwrite(NumhotANDcountry_result,'/root/autodl-tmp/humPoulResult/data/Hot_data/NumhotANDcountry_result.csv')








#——————————————————————————————————————————————————————————-----------------
basePath<-"/root/autodl-tmp/humPoulResult/data/"

world.map <- rnaturalearth::ne_countries(returnclass = "sf") |> dplyr::filter(continent != "Antarctica")
globalCountry <- vect(world.map) 
globalRaster <- rast(vals=1:259200,nrows=360, ncols=720,xmin=-180, xmax=180,ymin=-90, ymax=90,crs=crs)
coast <- ne_coastline(scale = "small", returnclass = "sf")
crs <- '+proj=longlat +datum=WGS84'
allDf <- fread(paste0(basePath,'allDf786_reclass.csv'))
speciesPixelNumPath <- list.files('/root/autodl-tmp/humPoulResult/data/single_model',pattern = '.tif',full.names = T)
spName <- basename(speciesPixelNumPath) %>% str_sub(.,1,-5)                      
speciesPixelNumPath2 <- speciesPixelNumPath[spName%in%allDf$LatName]

#Fig4 Mantel test------------
#计算分功能群活动熵###############



# #分功能群SH计算活动熵
# allDf2 <- fread('/root/autodl-tmp/zyresult/allDf786_reclass.csv')
# allDf2SH <- allDf2[allDf2$`Host`== "Suspected Host",]
# spName <- basename(speciesPixelNumPath2) %>% str_sub(.,1,-5)
# for (fam in unique(allDf2SH$`Family`)) {
#   famName <- allDf2SH[allDf2SH$`Family`==fam,]
#   paths <- speciesPixelNumPath2[spName%in%famName$LatName]
#   monthNum <- rast(paths) %>% sum(na.rm=T)
#   calEntropy <- lapply(paths, function(x){
#     r <- rast(x) %>% sum(.,na.rm=T)
#     pi <- r/monthNum
#     y <- -pi*log(pi)
#     names(y) <- str_sub(basename(x),1,-5)
#     return(y)
#   })
#   actEntropy <- rast(calEntropy) %>% sum(na.rm = T)
#   #writeRaster(actEntropy,paste0('/root/autodl-tmp/humPoulResult/data/AE_data/AE_',fuc,'_SH.tif'))
# }


#分功能群All-CH计算活动熵
allDf2 <- fread('/root/autodl-tmp/zyresult/allDf786_reclass.csv')
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
  #writeRaster(actEntropy,paste0('/root/autodl-tmp/humPoulResult/data/AE_data/AE_',fuc,'_CH.tif'))
}

# #分功能群SH计算活动熵
# allDf2 <- fread('/root/autodl-tmp/zyresult/allDf786_reclass.csv')
# allDf2SH <- allDf2[allDf2$`Host`== "Suspected Host",]
# spName <- basename(speciesPixelNumPath2) %>% str_sub(.,1,-5)
# for (fuc in unique(allDf2SH$`Functional Group`)) {
#   fucName <- allDf2SH[allDf2SH$`Functional Group`==fuc,]
#   paths <- speciesPixelNumPath2[spName%in%fucName$LatName]
#   monthNum <- rast(paths) %>% sum(na.rm=T)
#   calEntropy <- lapply(paths, function(x){
#     r <- rast(x) %>% sum(.,na.rm=T)
#     pi <- r/monthNum
#     y <- -pi*log(pi)
#     names(y) <- str_sub(basename(x),1,-5)
#     return(y)
#   })
#   actEntropy <- rast(calEntropy) %>% sum(na.rm = T)
#   writeRaster(actEntropy,paste0('/root/autodl-tmp/humPoulResult/data/AE_data/AE_',fuc,'_SH.tif'))
# }





#c########
###合并环境数据和活动熵数据
#1.环境数据
AE<-rast('/root/autodl-tmp/humPoulResult/data/Mantel_data/corData/AE.tif') %>% resample(globalRaster)%>% mask(globalCountry) %>% crop(ext(globalRaster))
spNum<-rast('/root/autodl-tmp/humPoulResult/data/Mantel_data/corData/spNumTif.tif') %>% resample(globalRaster)%>% mask(globalCountry)%>% crop(ext(globalRaster))
spNum_CV<-rast('/root/autodl-tmp/humPoulResult/data/Mantel_data/corData/sp_cv.tif') %>% resample(globalRaster)%>% mask(globalCountry)%>% crop(ext(globalRaster))
Allmonth<-rast('/root/autodl-tmp/humPoulResult/data/Mantel_data/corData/allMonth.tif') %>% resample(globalRaster)%>% mask(globalCountry)%>% crop(ext(globalRaster))
pop_log<-rast("/root/autodl-tmp/humPoulResult/data/Mantel_data/corData/pop_log.tif") %>% resample(globalRaster)%>% mask(globalCountry) %>% crop(ext(globalRaster))
poul_log<-rast("/root/autodl-tmp/humPoulResult/data/Mantel_data/corData/poul_log.tif") %>% resample(globalRaster)%>% mask(globalCountry) %>% crop(ext(globalRaster))
GPP<-rast('/root/autodl-tmp/humPoulResult/data/Mantel_data/corData/GPP.tif') %>% resample(globalRaster)%>% mask(globalCountry)%>% crop(ext(globalRaster))
GPP_CV<-rast('/root/autodl-tmp/humPoulResult/data/Mantel_data/corData/GPP_CV.tif') %>% resample(globalRaster)%>% mask(globalCountry)%>% crop(ext(globalRaster))
below0Days<-rast('/root/autodl-tmp/humPoulResult/data/Mantel_data/corData/below0Days.tif') %>% resample(globalRaster)%>% mask(globalCountry) %>% crop(ext(globalRaster))
bio1temp<-rast('/root/autodl-tmp/humPoulResult/data/Mantel_data/corData/bio1temp.tif') %>% resample(globalRaster)%>% mask(globalCountry) %>% crop(ext(globalRaster))
bio12prec<-rast('/root/autodl-tmp/humPoulResult/data/Mantel_data/corData/bio12prec.tif') %>% resample(globalRaster)%>% mask(globalCountry) %>% crop(ext(globalRaster))
# otherClimateData<-rast('/root/autodl-tmp/humPoulResult/data/Mantel_data/corData/otherClimateData.tif') %>% resample(globalRaster)%>% mask(globalCountry) %>% crop(ext(globalRaster))
# otherClimateData<-otherClimateData[[c("NDVI", "NDWI","LAI","Npp","nightLight", "roadDensity")]]

newclimateBirdData<-rast('/root/autodl-tmp/humPoulResult/data/Mantel_data/corData/climateBirdData.tif')%>% resample(globalRaster)%>% mask(globalCountry) %>% crop(ext(globalRaster))
otherClimateData<-newclimateBirdData[[c("NDVI", "NDWI","LAI","Npp","nightLight", "roadDensity")]]

Global_corData <- c(AE, spNum, Allmonth,spNum_CV,
                    bio12prec,bio1temp, below0Days,
                    GPP, GPP_CV ,otherClimateData, pop_log, poul_log)

globalSHP <- vect('/root/autodl-tmp/zyresult/Con_popentrpoul_sf_EU.shp')
globalSHP2 <- terra::aggregate(globalSHP,'name_ec')
Contry5c <-c('China','India','European Country','Nigeria','United States')
countryRaster <- subset(globalSHP2, globalSHP2$name_ec %in% Contry5c ) %>% rasterize(.,globalRaster,field='name_ec')

Global_corData_df <- c(Global_corData,countryRaster) %>% as.data.frame(xy=T)
names(Global_corData_df)<-c('Lontitude', 'Latitude' ,'AE','Species Richness', "Cumulative Species-months", 'Species Richness CV',
                            'Precipitation', 'Temperature', 'Frost Days',
                            'GPP', 'GPP CV',
                            "NDVI", "NDWI","LAI","Npp", 
                            'Light Density', 'Road Density', 
                            'Population(Log)', 'Poultry(Log)',
                            'Country')

Global_corData_df <- Global_corData_df %>%
  dplyr::select(#'AE','Species Richness', "Cumulative Species-months", 'Species Richness CV',
         'Lontitude', 'Latitude' ,'Precipitation', 'Temperature', 'Frost Days',
         'GPP', 'GPP CV',"Npp", "NDVI", "NDWI","LAI",
         'Light Density', 'Road Density', 'Population(Log)', 'Poultry(Log)',
         'Country')
summary(Global_corData_df)


#2.功能群活动熵数据
listFunc<-list.files('/root/autodl-tmp/humPoulResult/data/Mantel_data/birdFuncAE',full.names = T)
AE_Func<-rast(listFunc)%>% resample(globalRaster)%>% mask(globalCountry) %>% crop(ext(globalRaster))

globalSHP <- vect('/root/autodl-tmp/zyresult/Con_popentrpoul_sf_EU.shp')
globalSHP2 <- terra::aggregate(globalSHP,'name_ec')
Contry5c <-c('China','India','European Country','Nigeria','United States')
countryRaster <- subset(globalSHP2, globalSHP2$name_ec %in% Contry5c ) %>% rasterize(.,globalRaster,field='name_ec')
Global_birdEntropy_df <- c(AE_Func,countryRaster) %>% as.data.frame()

names(Global_birdEntropy_df) <- c('wadingbirds','unclassified','seabirds','shorebirds','waterfowls',"Country")
Global_birdEntropy_df <- Global_birdEntropy_df%>%
  dplyr::select('seabirds','shorebirds','waterfowls','wadingbirds','unclassified',"Country")
summary(Global_birdEntropy_df)


#计算分国家的mantel结果------------
#Global_birdEntropy_df5<-subset(Global_birdEntropy_df, Global_birdEntropy_df$Country %in%`Contry5c`)
Contry5c <-c('China','India','European Country','Nigeria','United States')
All5corData_result <- data.table()
All5mantel_result <- data.table()
for (con in Contry5c) {
  con_corData_df <- subset(Global_corData_df,Global_corData_df$Country==con) %>% dplyr::select(-ncol(.))
  con_birdEntropy_df <- subset(Global_birdEntropy_df,Global_birdEntropy_df$Country==con,) %>% dplyr::select(-ncol(.))
  #去掉有NA值的行(否则相关性分析会有NA)
  na_rows_cor <- which(rowSums(is.na(con_corData_df)) > 0)
  na_rows_bird <- which(rowSums(is.na(con_birdEntropy_df)) > 0)
  na_rows <- unique(c(na_rows_cor, na_rows_bird))
  con_corData_df_clean <- con_corData_df[-na_rows, ]
  con_birdEntropy_df_clean <- con_birdEntropy_df[-na_rows, ]
  #输出环境数据
  con_birdEntropy_df_clean2<-mutate(con_corData_df_clean, Country=con)
  All5corData_result<-rbind(All5corData_result,con_birdEntropy_df_clean2)
  #求mantel
  mantel <- mantel_test(con_birdEntropy_df_clean,    ## 影响因子数据 con_birdEntropy_df_clean
                        con_corData_df_clean,   ## 分类数据 con_corData_df_clean
                        spec_select = list(
                          seabirds = 1,
                          shorebirds = 2,
                          wadingbirds = 3,
                          waterfowls = 4,
                          unclassified = 5
                        )) %>%
    mutate(rd = cut(r, breaks = c(-Inf, 0, 0.4, 0.5, 0.6, 0.7, 0.8, Inf),
                    labels = c("< 0", "< 0.4", "0.4 - 0.5","0.5 - 0.6", "0.6 - 0.7", "0.7 - 0.8", ">= 0.8")),
           pd = cut(p, breaks = c(-Inf, 0.01, 0.05, Inf),
                    labels = c("< 0.01", "0.01 - 0.05", ">= 0.05")),
           Country=con)
  All5mantel_result <- rbind(All5mantel_result,mantel)
}

fwrite(All5corData_result,'/root/autodl-tmp/humPoulResult/data/Mantel_data/Mantel_result/All5corData_result0526.csv')
fwrite(All5mantel_result,'/root/autodl-tmp/humPoulResult/data/Mantel_data/Mantel_result/All5mantel_result0526.csv')


#画图China----------
library(vegan)
library(linkET)

China_corData_df_clean<-fread('/root/autodl-tmp/humPoulResult/data/Mantel_data/Mantel_result/All5corData_result0526.csv')%>%
  subset(.,.$Country=='China') %>% dplyr::select(-ncol(.))
China_mantel <- fread('/root/autodl-tmp/humPoulResult/data/Mantel_data/Mantel_result/All5mantel_result0526.csv')%>%
  subset(.,.$Country=='China') %>% dplyr::select(-ncol(.))

head(China_mantel)
China_mantel2<-subset(China_mantel, r>0.4)


# China_corData_df_clean2<-China_corData_df_clean[,c('AE','Species Richness', "Cumulative Species-months", 'Species Richness CV',
#                                                    'Lontitude', 'Latitude' , 'Precipitation', 'Temperature', 'Frost Days',
#                                                    'GPP', 'GPP CV', 'NPP', 
#                                                    'NDVI', 'NDWI', 'LAI', 
#                                                    'Light Density', 'Road Density', 
#                                                    'Human density(Log)', 'Poultry density(Log)')]
mdata <- correlate(China_corData_df_clean)
mdata_df <- as.data.frame(mdata$r)
#fwrite(mdata_df,'/root/autodl-tmp/zyresult/Mantel_Data/China_Correlate.csv')

qcorrplot(mdata,
          type = "lower", # 热图展示下半部分
          diag = F,
)+
  geom_square()+
  scale_fill_gradientn(colours = RColorBrewer::brewer.pal(11, "RdBu"))

China_mantel2 <- China_mantel2 %>%
  mutate(spec = case_when(
    spec == "seabirds" ~ "Seabirds",
    spec == "shorebirds" ~ "Shorebirds",
    spec == "wadingbirds" ~ "Large wading birds",
    spec == "waterfowls" ~ "Waterfowls",
    spec == "unclassified" ~ "Others",
    # 添加更多的替换规则
    TRUE ~ spec  # 默认情况下保持原值
  ))


qcorrplot(China_corData_df_clean, type = "lower", diag = FALSE, is_corr = TRUE) +  #China_corData_df_clean
  geom_square() +   ## 相关性热图的形状
  ## 
  geom_couple(aes(colour = rd, linetype = pd), size=1,
              data = China_mantel, label.size = 6,
              curvature = nice_curvature()) +
  ## 颜色参数调整
  scale_fill_gradientn(colours = RColorBrewer::brewer.pal(11, "RdBu")) +
  scale_linetype_manual(values = c("< 0.01" = "solid", 
                                   "0.01 - 0.05" = "dashed", 
                                   ">= 0.05" = "blank")) + # "blank"表示不显示线
  scale_size_manual(values = c(1, 2)) +
  #scale_colour_manual(values = color_pal(2)) +
  scale_colour_manual(values = c("#e8f0ff","#bad4ff","#5db2f6","#3871d1","#2357a9"))+
  guides(
    # size = guide_legend(title = "Mantel's p",
    #                          override.aes = list(colour = "grey35"),
    #                          order = 2),
    colour = guide_legend(title = "Mantel's r",
                          override.aes = list(size = 3), 
                          order = 1),
    linetype = guide_legend(title = "Mantel's p",
                            override.aes = list(size = 3), 
                            order = 3), # 添加linetype的设置
    fill = guide_colorbar(title = "Pearson's r", order = 3))+
  theme(text = element_text(size = 15),
        axis.text = element_text(size = 15),)









#a#############

All5corData_result<- fread('/root/autodl-tmp/humPoulResult/data/Mantel_data/Mantel_result/All5corData_result0526.csv')


# 加载必要的库
library(dplyr)
library(tidyr)
library(purrr)
library(corrr)

# # 将数据按照Country分组并嵌套，然后计算每个分组的相关性矩阵
# correlation_results <- All5corData_result %>%
#   group_by(Country) %>%
#   nest() %>%
#   mutate(correlation = map(data, ~correlate(select_if(.x, is.numeric)))) %>%
#   select(Country, correlation) %>%
#   unnest(correlation)

# 将数据按照Country分组并嵌套，然后计算每个分组的相关性矩阵
correlation_results <- All5corData_result %>%
  group_by(Country) %>%
  nest() %>%
  mutate(correlation = map(data, ~ {
    # 排除标准差为零的变量
    numeric_data <- select_if(.x, is.numeric)
    numeric_data <- numeric_data[, sapply(numeric_data, function(x) sd(x) != 0)]
    # 计算相关性矩阵
    correlate(numeric_data)
  })) %>%
  select(Country, correlation) %>%
  unnest(correlation)


# 首先将相关性矩阵转换为长格式
long_correlation_results <- correlation_results %>%
  select(-term) %>%  # 排除term列
  mutate(id = row_number()) %>%
  pivot_longer(cols = -c(Country, id), names_to = "variable", values_to = "correlation") %>%
  filter(variable %in% c("Lontitude", "Latitude", "Precipitation", "Temperature", "Frost Days", 
                         "GPP", "CV of GPP", "Npp", "NDVI", "NDWI", "LAI",
                         "Light Density", "Road Density", "Population(Log)", "Poultry(Log)"))
long_correlation_results <- correlation_results %>%
  select(-term) %>%  # 排除term列
  pivot_longer(cols = -Country, names_to = "variable", values_to = "correlation")
# 分组绘制柱形图
# 第一组：Longitude, Latitude, Precipitation, Temperature, Frost Days
ggplot(long_correlation_results %>% filter(variable %in% c("Lontitude", "Latitude", "Precipitation", "Temperature", "Frost Days")), 
       aes(x = variable, y = correlation, fill = Country)) +
  geom_bar(stat = "identity", position = "dodge") +
  facet_wrap(~ Country) +
  labs(title = "Group 1: Environmental Factors", x = "Variable", y = "Correlation") +
  theme_minimal()

# 第二组：GPP, GPPCV, NPP, NDVI, NDWI, LAI
ggplot(long_correlation_results %>% filter(variable %in% c("GPP", "GPPCV", "NPP", "NDVI", "NDWI", "LAI")), 
       aes(x = variable, y = correlation, fill = Country)) +
  geom_bar(stat = "identity", position = "dodge") +
  facet_wrap(~ Country) +
  labs(title = "Group 2: Vegetation Indices and Productivity", x = "Variable", y = "Correlation") +
  theme_minimal()

# 第三组：Light Density, Road Density, Population, Poultry
ggplot(long_correlation_results %>% filter(variable %in% c("Light Density", "Road Density", "Population", "Poultry")), 
       aes(x = variable, y = correlation, fill = Country)) +
  geom_bar(stat = "identity", position = "dodge") +
  facet_wrap(~ Country) +
  labs(title = "Group 3: Human Activity and Livestock", x = "Variable", y = "Correlation") +
  theme_minimal()


#####地理与气候因子的数据
climatedata <- data.frame(
  Country = c("China", "India", "EU", "Nigeria", "US"),
  Longitude = c(0.39523392, -0.138386625, 0.162993606, -0.223024804, 0.475806627),
  Latitude = c(-0.409592792, -0.390704017, 0.217154251, 0.355986578, -0.307689008),
  Precipitation = c(0.74364176, 0.15271758, -0.13481555, -0.225767993, 0.712091197),
  Temperature = c(0.744273909, 0.661957413, 0.042850816, 0.346339967, 0.417501088),
  Subzero_Temperature_Days = c(-0.640081483, -0.565841963, -0.139838478, NA, -0.464274536)
)

# 将数据从宽格式转换为长格式
climatedata_long <- climatedata %>% 
  gather(key = "Attribute", value = "Value", -Country)
# 对每个国家的数据进行排序
climatedata_long2 <- climatedata_long %>% 
  group_by(Attribute) %>% 
  mutate(Rank = rank(-Value)) %>% 
  ungroup() %>% 
  arrange(Attribute, Rank)%>%
  group_by(Country) %>% 
  mutate(Rank2 = rank(-Value)) %>% 
  ungroup() %>% 
  arrange(Country, Rank2)

#climatedata_long2$Attribute<-factor(climatedata_long2$Attribute,levels = c("Longitude","Latitude","Precipitation","Temperature","Subzero_Temperature_Days"))

climatedata_long2 <- climatedata_long2 %>%
  mutate(Attribute = case_when(
    Attribute == "Subzero_Temperature_Days" ~ "Frost Days",
    TRUE ~ Attribute  # 默认情况下保持原值
  ))

climatedata_long2$Attribute<-factor(climatedata_long2$Attribute,levels = c("Longitude","Latitude","Precipitation","Temperature","Frost Days"))

#####植被指标的数据

vegetationdata <- data.frame(
  Country = c("China", "India", "European Country", "Nigeria", "United States"),
  GPP = c(0.657613199, 0.158686035, -0.023214674, -0.405162254, 0.69214518),
  CV_of_GPP = c(-0.65866009, -0.509137502, 0.051341938, 0.31214534, -0.428621043),
  NPP = c(0.628789993, -0.044672024, 0.069611827, -0.430419677, 0.714024786),
  NDVI = c(0.646875931, 0.220660987, 0.131865748, -0.370403682, 0.705716504),
  NDWI = c(-0.444388109, -0.36178946, 0.083415481, 0.158774818, -0.611855474),
  LAI = c(0.440519166, -0.071703767, -0.058383767, -0.423725974, 0.5237572)
)

# 将数据从宽格式转换为长格式
vegetationdata_long <- vegetationdata %>% 
  gather(key = "Attribute", value = "Value", -Country)
# 对每个国家的数据进行排序
vegetationdata_long2 <- vegetationdata_long %>% 
  group_by(Attribute) %>% 
  mutate(Rank = rank(-Value)) %>% 
  ungroup() %>% 
  arrange(Attribute, Rank)%>%
  group_by(Country) %>% 
  mutate(Rank2 = rank(-Value)) %>% 
  ungroup() %>% 
  arrange(Country, Rank2)


vegetationdata_long2 <- vegetationdata_long2 %>%
  mutate(Attribute = case_when(
    Attribute == "CV_of_GPP" ~ "GPP CV",
    TRUE ~ Attribute  # 默认情况下保持原值
  ))


vegetationdata_long2$Attribute<-factor(vegetationdata_long2$Attribute,levels = c("GPP","GPP CV","NPP","NDVI","NDWI","LAI"))





#####人类活动
humandata <- data.frame(
  Country = c("China", "United States", "India", "European Country", "Nigeria"),
  Light_Density = c(0.41433085, 0.494451945, 0.326066598, 0.09459512, 0.055669773),
  Road_Density = c(0.404520561, 0.377699523, 0.212541503, 0.268782901, 0.171802425),
  Population_Log = c(0.675935177, 0.637477794, 0.654811034, 0.283260059, 0.083748191),
  Poultry_Log = c(0.72515455, 0.561622506, 0.35928713, 0.122675018, 0.10380714)
)

# 将数据从宽格式转换为长格式
humandata_long <- humandata %>% 
  gather(key = "Attribute", value = "Value", -Country)
# 对每个国家的数据进行排序
humandata_long2 <- humandata_long %>% 
  group_by(Attribute) %>% 
  mutate(Rank = rank(-Value)) %>% 
  ungroup() %>% 
  arrange(Attribute, Rank)%>%
  group_by(Country) %>% 
  mutate(Rank2 = rank(-Value)) %>% 
  ungroup() %>% 
  arrange(Country, Rank2)

#humandata_long2$Attribute<- factor(humandata_long2$Attribute, levels = c("Light_Density","Road_Density","Population_Log","Poultry_Log"))
humandata_long2 <- humandata_long2 %>%
  mutate(Attribute = case_when(
    Attribute == "Light_Density" ~ "Light Density",
    Attribute == "Road_Density" ~ "Road Density",
    Attribute == "Population_Log" ~ "Human density(Log)",
    Attribute == "Poultry_Log" ~ "Poultry density(Log)",
    TRUE ~ Attribute  # 默认情况下保持原值
  ))

humandata_long2$Attribute<- factor(humandata_long2$Attribute, levels = c("Light Density","Road Density","Human density(Log)","Poultry density(Log)"))

########## _绘制直方图
line_data1 <- data.frame(Attribute = rep(unique(climatedata_long2$Attribute), each = 1),
                         y = c(-0.221247174689965, -0.109375908212166, 0.268893708, 0.263637402, -0.341974174))
p1<-ggplot(climatedata_long2, aes(x = Attribute, y = Value, fill = reorder(Country, -Value))) +
  geom_bar(stat = "identity", position = "dodge") +
  #facet_grid(~ Attribute, scales = "free") +
  geom_segment(aes(x = 0.55, xend = 1.45, y = -0.221247174689965), linetype = "dashed", color = "#ffbcbc",size=0.3) +   #Lon
  geom_segment(aes(x = 1.55, xend = 2.45, y = -0.109375908212166), linetype = "dashed", color = "#ffbcbc",size=0.3) +   #Lat
  geom_segment(aes(x = 2.55, xend = 3.45, y = 0.268893708), linetype = "dashed", color = "#ffbcbc",size=0.3) +   #Per
  geom_segment(aes(x = 3.55, xend = 4.45, y = 0.263637402), linetype = "dashed", color = "#ffbcbc",size=0.3) +   #Temp
  geom_segment(aes(x = 4.55, xend = 5.45, y = -0.341974174), linetype = "dashed", color = "#ffbcbc",size=0.3) +   #Temp<0days
  # geom_segment(data = line_data1,
  #                aes(x = 0.55, xend = 1.45, y = y, yend = y),
  #                linetype = "dashed", color = "#ffbcbc", size = 0.3) +
  scale_fill_manual(values = c("India"="#377EB8", "China"="#E41A1C", "EU"="#4DAF4A", "Nigeria"="#984EA3", "US"="#f7a156")) +  #FF7F00
  theme_bw() +
  #scale_y_continuous(limits = c(-0.5,1))+
  labs(title = "Geographic and Climate Factors", x = NULL,y = "Person's r \n ( ~ AE )", fill = NULL) +
  theme(axis.text.x = element_text(size=12),   #angle = 15, hjust = 0.9,
        text = element_text(size=14),
        legend.position = "bottom",
        plot.title = element_text(hjust = 0.5,size=14))

p2<-ggplot(vegetationdata_long2, aes(x = Attribute, y = Value, fill = reorder(Country, Rank))) +
  geom_bar(stat = "identity", position = "dodge") +
  geom_segment(aes(x = 0.55, xend = 1.45, y = 0.352070516), linetype = "dashed", color = "#ffbcbc",size=0.3) +   #GPP
  geom_segment(aes(x = 1.55, xend = 2.45, y = -0.287011164), linetype = "dashed", color = "#ffbcbc",size=0.3) +   #GPPCV
  geom_segment(aes(x = 2.55, xend = 3.45, y = 0.292114987), linetype = "dashed", color = "#ffbcbc",size=0.3) +   #NPP
  geom_segment(aes(x = 3.55, xend = 4.45, y = 0.435593707), linetype = "dashed", color = "#ffbcbc",size=0.3) +   #NDVI
  geom_segment(aes(x = 4.55, xend = 5.45, y = -0.367204322), linetype = "dashed", color = "#ffbcbc",size=0.3) +   #NDWI
  geom_segment(aes(x = 5.55, xend = 6.45, y = 0.130518149), linetype = "dashed", color = "#ffbcbc",size=0.3) +   #LAI
  #facet_wrap(~ Country, scales = "free") +
  scale_fill_manual(values = c("India"="#377EB8", "China"="#E41A1C", "European Country"="#4DAF4A", "Nigeria"="#984EA3", "United States"="#f7a156")) +  #FF7F00
  #scale_fill_manual(values = c("#377EB8", "#E41A1C", "#4DAF4A", "#984EA3", "#FF7F00")) +  #FF7F00
  #scale_fill_manual(values = c("India"="#4C8BC0", "China"="#C75C64", "European Country"="#32A852", "Nigeria"="#984EA3", "United States"="#F0B57D")) +  #F0B57D
  theme_bw() +
  #scale_y_continuous(limits = c(0,0.8))+
  labs(title = "Vegetation Indicators",x = NULL, y = NULL, fill = NULL) +
  theme(axis.text.x = element_text(size=12),  #angle = 15, hjust = 0.9,
        text = element_text(size=14),
        legend.position = "bottom",
        plot.title = element_text(hjust = 0.5,size=14))

p3<-ggplot(humandata_long2, aes(x = Attribute, y = Value, fill = reorder(Country, Rank))) +
  geom_bar(stat = "identity", position = "dodge") +
  #facet_wrap(~ Country, scales = "free") +
  geom_segment(aes(x = 0.55, xend = 1.45, y = 0.259050496), linetype = "dashed", color = "#ffbcbc",size=0.3) +   #Light
  geom_segment(aes(x = 1.55, xend = 2.45, y = 0.241559629), linetype = "dashed", color = "#ffbcbc",size=0.3) +   #road
  geom_segment(aes(x = 2.55, xend = 3.45, y = 0.354607465), linetype = "dashed", color = "#ffbcbc",size=0.3) +   #Pop
  geom_segment(aes(x = 3.55, xend = 4.45, y = 0.267712171), linetype = "dashed", color = "#ffbcbc",size=0.3) +   #poul
  scale_fill_manual(values = c("India"="#377EB8", "China"="#E41A1C", "European Country"="#4DAF4A", "Nigeria"="#984EA3", "United States"="#f7a156")) +  #FF7F00
  #scale_fill_manual(values = c("#377EB8", "#E41A1C", "#4DAF4A", "#984EA3", "#FF7F00")) +  #FF7F00
  #scale_fill_manual(values = c("India"="#4C8BC0", "China"="#C75C64", "European Country"="#32A852", "Nigeria"="#984EA3", "United States"="#F0B57D")) +  #F0B57D
  theme_bw() +
  scale_y_continuous(limits = c(0,0.8))+
  labs(title = "Human Activity Indicators",x = NULL, y = NULL, fill = NULL) +
  theme(axis.text.x = element_text(size=12),   #angle = 15, hjust = 0.9,
        text = element_text(size=14),
        legend.position = "bottom",
        plot.title = element_text(hjust = 0.5,size=14))



library(patchwork)
p1 + p2 + p3 + plot_layout(widths = c(2/3, 4/5, 8/15))



#b#######################

popd2015<-rast("/root/autodl-tmp/全球人口/GWP_v4/popd2015_30.tif") %>% resample(globalRaster)%>% mask(globalCountry)
Entropy<-rast(paste0(basePath,'AE_data/AE.tif'))%>% mask(globalCountry)
poul2015<-rast("/root/autodl-tmp/zyresult/Poultry_duckchic.tif")%>% resample(globalRaster)%>% mask(globalCountry)


globalSHP <- vect('/root/autodl-tmp/zyresult/Con_popentrpoul_sf_EU.shp')
globalSHP2 <- terra::aggregate(globalSHP,'name_ec')
Con_5 <- subset(globalSHP2, globalSHP2$name_ec %in% c("China", "India", "European Country", "Nigeria", "United States"))
cRaster <- rasterize(Con_5,Entropy,field='name_ec')
All5_df <- c(Entropy,popd2015,poul2015,cRaster) %>% terra::as.data.frame() %>%na.omit()
names(All5_df)<-c("entr","pop","poul","country")




#####人口线性
# 过滤掉pop值小于0.0001的数据行
All5_df1 <- All5_df[!(All5_df$pop < 0.0001),]

p1<-ggplot(data = All5_df1, aes(x = entr, y = log(pop), color = country)) +
  geom_point(size=1) +  #,shape=21
  scale_color_manual(values = c("India"="#377EB8", "China"="#E41A1C", "European Country"="#4DAF4A", "Nigeria"="#984EA3", "United States"="#f7a156")) +  #FF7F00
  #scale_color_manual(values = c("#377EB8", "#E41A1C", "#4DAF4A", "#984EA3", "#FF7F00")) +
  geom_smooth(method = "lm",linewidth=0.5, se = T, lty=5,aes(color = country)) +
  #geom_smooth(method = "lm", se = T, color = "black", linetype = "dashed") +  # 全部点的趋势线，黑色
  #geom_smooth(method = "lm",linewidth=0.5, se = T, lty=5,color = "black") +  #分界
  stat_poly_eq(aes(label = after_stat(paste(eq.label, rr.label,p.value.label, sep = "~~~"))),  #,p.value.label  eq.label,
               formula = y ~ x,
               coef.digits = 3, rr.digits = 2, parse = TRUE, size =3) +
  #facet_wrap(~ country) +  #如果要分界
  labs(x = NULL, y = "Human density (Log)") +   #Waterbird Activity Entropy
  ylim(-10,10)+
  theme_bw()+ 
  theme(legend.position = "none",
        legend.title = element_blank(),
        axis.text = element_text(size = 14),
        plot.title = element_text(size = 14, hjust = 0.5),
        axis.title = element_text(size = 16))
# stat_poly_eq(aes(label = after_stat(paste(eq.label, rr.label, sep = "~~~"))),
#              formula = y ~ x,
#              coef.digits = 3, rr.digits = 2, parse = TRUE, size = 3,
#              label.x = "left", label.y = "top",
#              color = "black")
print(p1)

#不log
ggplot(data = All5_df1, aes(x = sum, y = pop, color = country)) +
  geom_point(size=1) +  #,shape=21
  scale_color_manual(values = c("#377EB8", "#E41A1C", "#4DAF4A", "#984EA3", "#FF7F00")) +
  geom_smooth(method = "lm",linewidth=0.5, se = T, lty=5,aes(color = country)) +
  #geom_smooth(method = "lm",linewidth=0.5, se = T, lty=5,color = "black") +  #分界
  stat_poly_eq(aes(label = after_stat(paste(eq.label, rr.label, sep = "~~~"))),
               formula = y ~ x, 
               coef.digits = 3, rr.digits = 2, parse = TRUE, size = 3) +
  #facet_wrap(~ country) +  #如果要分界
  labs(x = "Entropy", y = "Population") +
  theme_minimal() +
  theme(legend.position = "bottom",
        legend.title = element_blank(),
        plot.title = element_text(size = 12, hjust = 0.5),
        axis.title = element_text(size = 14))



#统计
india_count <- nrow(All5_df1[All5_df1$country == 'India', ]); india_count  # 1195
China_count <- nrow(All5_df1[All5_df1$country == 'China', ]); China_count  # 3798
EU_count <- nrow(All5_df1[All5_df1$country == 'European Country', ]); EU_count  # 2273
Nigeria_count <- nrow(All5_df1[All5_df1$country == 'Nigeria', ]); Nigeria_count  # 330
US_count <- nrow(All5_df1[All5_df1$country == 'United States', ]); US_count  # 3770

# 拟合线性模型
model <- lm(log(pop) ~ entr, data = All5_df1[All5_df1$country == 'India', ])
model_summary <- summary(model);model_summary  # 获取模型摘要
b <- model_summary$coefficients[2, 1] ; b # 提取回归系数b和标准误差SE
SE <- model_summary$coefficients[2, 2];SE
z <- b / SE   ;z
p <- model_summary$coefficients[2, 4];  p
# cor(x, y, method = "pearson")
# stat_text <- paste("b =", round(b, 2), "\nSE =", round(SE, 2), "\nz =", round(z, 2), "\np =", round(p, 2));stat_text


model <- lm(log(pop) ~ entr, data = All5_df1[All5_df1$country == 'China', ])
model_summary <- summary(model);model_summary  # 获取模型摘要
b <- model_summary$coefficients[2, 1] ; b # 提取回归系数b和标准误差SE
SE <- model_summary$coefficients[2, 2];SE
z <- b / SE   ;z
p <- model_summary$coefficients[2, 4];  p
stat_text <- paste("b =", round(b, 2), "\nSE =", round(SE, 2), "\nz =", round(z, 2), "\np =", round(p, 2));stat_text

model <- lm(log(pop) ~ entr, data = All5_df1[All5_df1$country == 'European Country', ])
model_summary <- summary(model);model_summary  # 获取模型摘要
b <- model_summary$coefficients[2, 1] ; b # 提取回归系数b和标准误差SE
SE <- model_summary$coefficients[2, 2];SE
z <- b / SE   ;z
p <- model_summary$coefficients[2, 4];  p
stat_text <- paste("b =", round(b, 2), "\nSE =", round(SE, 2), "\nz =", round(z, 2), "\np =", round(p, 2));stat_text

model <- lm(log(pop) ~ entr, data = All5_df1[All5_df1$country == 'Nigeria', ])
model_summary <- summary(model);model_summary  # 获取模型摘要
b <- model_summary$coefficients[2, 1] ; b # 提取回归系数b和标准误差SE
SE <- model_summary$coefficients[2, 2];SE
z <- b / SE   ;z
p <- model_summary$coefficients[2, 4];  p
stat_text <- paste("b =", round(b, 2), "\nSE =", round(SE, 2), "\nz =", round(z, 2), "\np =", round(p, 2));stat_text

model <- lm(log(pop) ~ entr, data = All5_df1[All5_df1$country == 'United States', ])
model_summary <- summary(model);model_summary  # 获取模型摘要
b <- model_summary$coefficients[2, 1] ; b # 提取回归系数b和标准误差SE
SE <- model_summary$coefficients[2, 2];SE
z <- b / SE   ;z
p <- model_summary$coefficients[2, 4];  p
stat_text <- paste("b =", round(b, 2), "\nSE =", round(SE, 2), "\nz =", round(z, 2), "\np =", round(p, 2));stat_text


#####家禽线性
# 过滤掉poul值小于0.0001
All5_df2 <- All5_df[!(All5_df$poul < 0.0001) ,]

p2<-ggplot(data = All5_df2, aes(x = entr, y = log(poul), color = country)) +
  geom_point(size=1) +  #,shape=21
  scale_color_manual(values = c("India"="#377EB8", "China"="#E41A1C", "European Country"="#4DAF4A", "Nigeria"="#984EA3", "United States"="#f7a156")) +  #FF7F00
  #scale_color_manual(values = c("#377EB8", "#E41A1C", "#4DAF4A", "#984EA3", "#FF7F00")) +
  geom_smooth(method = "lm",linewidth=0.5, se = T, lty=5,aes(color = country)) +
  #geom_smooth(method = "lm", se = T, color = "black", linetype = "dashed") +  # 全部点的趋势线，黑色
  #geom_smooth(method = "lm",linewidth=0.5, se = T, lty=5,color = "black") +  #分界
  stat_poly_eq(aes(label = after_stat(paste(eq.label, rr.label,p.value.label, sep = "~~~"))),
               formula = y ~ x,
               coef.digits = 3, rr.digits = 2, parse = TRUE, size = 3) +
  #facet_wrap(~ country) +  #如果要分界
  labs(x = "AE", y = "Poultry density (Log)") +
  #ylim(-7,15)+
  theme_bw() +
  theme(legend.position = "none",
        legend.title = element_blank(),
        axis.text = element_text(size = 12),
        plot.title = element_text(size = 12, hjust = 0.5),
        axis.title = element_text(size = 14))
# stat_poly_eq(aes(label = after_stat(paste(eq.label, rr.label, sep = "~~~"))),
#              formula = y ~ x,
#              coef.digits = 3, rr.digits = 2, parse = TRUE, size = 3,
#              label.x = "left", label.y = "top",
#              color = "black")
print(p2)

library(patchwork)
p1 + p2 + plot_layout(ncol = 1)

#不log
ggplot(data = All5_df2, aes(x = entr, y = poul, color = country)) +
  geom_point(size=1) +  #,shape=21
  scale_color_manual(values = c("#377EB8", "#E41A1C", "#4DAF4A", "#984EA3", "#FF7F00")) +
  geom_smooth(method = "lm",linewidth=0.5, se = T, lty=5,aes(color = country)) +
  #geom_smooth(method = "lm",linewidth=0.5, se = T, lty=5,color = "black") +  #分界
  stat_poly_eq(aes(label = after_stat(paste(eq.label, rr.label, sep = "~~~"))),
               formula = y ~ x, 
               coef.digits = 3, rr.digits = 2, parse = TRUE, size = 3) +
  #facet_wrap(~ country) +  #如果要分界
  labs(x = "Entropy", y = "Poultry") +
  theme_minimal() +
  theme(legend.position = "bottom",
        legend.title = element_blank(),
        plot.title = element_text(size = 14, hjust = 0.5),
        axis.title = element_text(size = 12))


#统计
india_count <- nrow(All5_df2[All5_df2$country == 'India', ]); india_count  # 1114
China_count <- nrow(All5_df2[All5_df2$country == 'China', ]); China_count  # 3817
EU_count <- nrow(All5_df2[All5_df2$country == 'European Country', ]); EU_count  # 2163
Nigeria_count <- nrow(All5_df2[All5_df2$country == 'Nigeria', ]); Nigeria_count  # 300
US_count <- nrow(All5_df2[All5_df2$country == 'United States', ]); US_count  # 4353

# 拟合线性模型
model <- lm(log(poul) ~ entr, data = All5_df2[All5_df2$country == 'India', ])
model_summary <- summary(model);model_summary  # 获取模型摘要
b <- model_summary$coefficients[2, 1] ; b # 提取回归系数b和标准误差SE
SE <- model_summary$coefficients[2, 2];SE
z <- b / SE   ;z
p <- model_summary$coefficients[2, 4];  p
stat_text <- paste("b =", round(b, 2), "\nSE =", round(SE, 2), "\nz =", round(z, 2), "\np =", round(p, 2));stat_text


model <- lm(log(poul) ~ entr, data = All5_df2[All5_df2$country == 'China', ])
model_summary <- summary(model);model_summary  # 获取模型摘要
b <- model_summary$coefficients[2, 1] ; b # 提取回归系数b和标准误差SE
SE <- model_summary$coefficients[2, 2];SE
z <- b / SE   ;z
p <- model_summary$coefficients[2, 4];  p
stat_text <- paste("b =", round(b, 2), "\nSE =", round(SE, 2), "\nz =", round(z, 2), "\np =", round(p, 2));stat_text

model <- lm(log(poul) ~ entr, data = All5_df2[All5_df2$country == 'European Country', ])
model_summary <- summary(model);model_summary  # 获取模型摘要
b <- model_summary$coefficients[2, 1] ; b # 提取回归系数b和标准误差SE
SE <- model_summary$coefficients[2, 2];SE
z <- b / SE   ;z
p <- model_summary$coefficients[2, 4];  p
stat_text <- paste("b =", round(b, 2), "\nSE =", round(SE, 2), "\nz =", round(z, 2), "\np =", round(p, 2));stat_text

model <- lm(log(poul) ~ entr, data = All5_df2[All5_df2$country == 'Nigeria', ])
model_summary <- summary(model);model_summary  # 获取模型摘要
b <- model_summary$coefficients[2, 1] ; b # 提取回归系数b和标准误差SE
SE <- model_summary$coefficients[2, 2];SE
z <- b / SE   ;z
p <- model_summary$coefficients[2, 4];  p
stat_text <- paste("b =", round(b, 2), "\nSE =", round(SE, 2), "\nz =", round(z, 2), "\np =", round(p, 2));stat_text

model <- lm(log(poul) ~ entr, data = All5_df2[All5_df2$country == 'United States', ])
model_summary <- summary(model);model_summary  # 获取模型摘要
b <- model_summary$coefficients[2, 1] ; b # 提取回归系数b和标准误差SE
SE <- model_summary$coefficients[2, 2];SE
z <- b / SE   ;z
p <- model_summary$coefficients[2, 4];  p
stat_text <- paste("b =", round(b, 2), "\nSE =", round(SE, 2), "\nz =", round(z, 2), "\np =", round(p, 2));stat_text




#——————————————————————————————————————————————————————————-----------------
basePath<-"/root/autodl-tmp/humPoulResult/data/"

world.map <- rnaturalearth::ne_countries(returnclass = "sf") |> dplyr::filter(continent != "Antarctica")
globalCountry <- vect(world.map) 
globalRaster <- rast(vals=1:259200,nrows=360, ncols=720,xmin=-180, xmax=180,ymin=-90, ymax=90,crs=crs)
coast <- ne_coastline(scale = "small", returnclass = "sf")
crs <- '+proj=longlat +datum=WGS84'
allDf <- fread(paste0(basePath,'allDf786_reclass.csv'))
speciesPixelNumPath <- list.files('/root/autodl-tmp/humPoulResult/data/single_model',pattern = '.tif',full.names = T)
spName <- basename(speciesPixelNumPath) %>% str_sub(.,1,-5)                      
speciesPixelNumPath2 <- speciesPixelNumPath[spName%in%allDf$LatName]




#Fig.5 ############
#a ： Functional Group###########
#2.总功能群活动熵数据
listFunc<-list.files('/root/autodl-tmp/humPoulResult/data/Mantel_data/birdFuncAE',full.names = T)
AE_Func<-rast(listFunc)%>% resample(globalRaster)%>% mask(globalCountry) %>% crop(ext(globalRaster))
pop_log<-rast("/root/autodl-tmp/humPoulResult/data/Mantel_data/corData/pop_log.tif") %>% resample(globalRaster)%>% mask(globalCountry) %>% crop(ext(globalRaster))
poul_log<-rast("/root/autodl-tmp/humPoulResult/data/Mantel_data/corData/poul_log.tif") %>% resample(globalRaster)%>% mask(globalCountry) %>% crop(ext(globalRaster))

pop<-rast("/root/autodl-tmp/全球人口/GWP_v4/popd2015_30.tif") %>% resample(globalRaster)%>% mask(globalCountry) %>% crop(ext(globalRaster))
poul<-rast("/root/autodl-tmp/zyresult/Poultry_duckchic.tif") %>% resample(globalRaster)%>% mask(globalCountry) %>% crop(ext(globalRaster))

globalSHP <- vect('/root/autodl-tmp/zyresult/Con_popentrpoul_sf_EU.shp')
globalSHP2 <- terra::aggregate(globalSHP,'name_ec')
Con_5Raster <- subset(globalSHP2, globalSHP2$name_ec %in% c("China", "India", "European Country", "Nigeria", "United States"))%>% rasterize(.,globalRaster,field='name_ec')

Global_birdEntropy_df <- c(AE_Func,pop_log,poul_log,pop,poul,Con_5Raster) %>% as.data.frame()

names(Global_birdEntropy_df) <- c('Large wading birds','Others','Seabirds','Shorebirds','Waterfowls','Population(Log)','Poltry(Log)','pop','poul',"Country")
Global_birdEntropy_df <- Global_birdEntropy_df%>%
  dplyr::select('Seabirds','Shorebirds','Waterfowls','Large wading birds','Others','Population(Log)','Poltry(Log)','pop','poul',"Country")
summary(Global_birdEntropy_df)
#fwrite(Global_birdEntropy_df,'/root/autodl-tmp/zyresult/Global_birdEntropy_df.csv')


#宿主功能群
listFuncCH<-list.files('/root/autodl-tmp/humPoulResult/data/AE_data',pattern = 'CH.tif',full.names = T)
AE_FuncCH<-rast(listFuncCH)%>% resample(globalRaster)%>% mask(globalCountry) %>% crop(ext(globalRaster))
pop_log<-rast("/root/autodl-tmp/humPoulResult/data/Mantel_data/corData/pop_log.tif") %>% resample(globalRaster)%>% mask(globalCountry) %>% crop(ext(globalRaster))
poul_log<-rast("/root/autodl-tmp/humPoulResult/data/Mantel_data/corData/poul_log.tif") %>% resample(globalRaster)%>% mask(globalCountry) %>% crop(ext(globalRaster))

pop<-rast("/root/autodl-tmp/全球人口/GWP_v4/popd2015_30.tif") %>% resample(globalRaster)%>% mask(globalCountry) %>% crop(ext(globalRaster))
poul<-rast("/root/autodl-tmp/zyresult/Poultry_duckchic.tif") %>% resample(globalRaster)%>% mask(globalCountry) %>% crop(ext(globalRaster))

globalSHP <- vect('/root/autodl-tmp/zyresult/Con_popentrpoul_sf_EU.shp')
globalSHP2 <- terra::aggregate(globalSHP,'name_ec')
Con_5Raster <- subset(globalSHP2, globalSHP2$name_ec %in% c("China", "India", "European Country", "Nigeria", "United States"))%>% rasterize(.,globalRaster,field='name_ec')

Global_CHEntropy_df <- c(AE_FuncCH,pop_log,poul_log,pop,poul,Con_5Raster) %>% as.data.frame()

names(Global_CHEntropy_df) <- c('Large wading birds','Others','Seabirds','Shorebirds','Waterfowls','Population(Log)','Poltry(Log)','pop','poul',"Country")
Global_CHEntropy_df <- Global_CHEntropy_df%>%
  dplyr::select('Seabirds','Shorebirds','Waterfowls','Large wading birds','Others','Population(Log)','Poltry(Log)','pop','poul',"Country")
summary(Global_CHEntropy_df)

#fwrite(Global_CHEntropy_df,'/root/autodl-tmp/zyresult/Global_CHEntropy_df.csv')




# Global_CHEntropy_df3 <- data.table::melt(data.table(Global_birdEntropy_df),id=6:10,measure=1:5)
# result_chpop <- Global_CHEntropy_df3[value>0,.(slope=summary(lm(`Population(Log)`~value))$coef[2],
#                                                r2=summary(lm(`Population(Log)`~value))$r.squared,
#                                                se=summary(lm(`Population(Log)`~value))$coef[4],
#                                                p=summary(lm(`Population(Log)`~value))$coef[8],
#                                                label='ch'),.(Country,variable)]
# 
# 

######循环：分别与人口\家禽进行线性拟合------------
Contry5c <-c('China','India','European Country','Nigeria','United States')

Global_birdEntropy_df5_long <- gather(Global_birdEntropy_df, key = "Variable", value = "Value",
                    -`Population(Log)`, -`Poltry(Log)`, -`pop`, -`poul`, -`Country`) %>%na.omit()
Global_birdEntropy_df5_long <- subset(Global_birdEntropy_df5_long, Value>0)
summary(Global_birdEntropy_df5_long)#ALL

Global_CHEntropy_df5_long <- gather(Global_CHEntropy_df, key = "Variable", value = "Value",
                                      -`Population(Log)`, -`Poltry(Log)`, -`pop`, -`poul`, -`Country`) %>%na.omit()
Global_CHEntropy_df5_long <- subset(Global_CHEntropy_df5_long, Value>0)
summary(Global_CHEntropy_df5_long)#CH



# 初始化一个空的数据框来存储结果
All5_poppoulresults_df <- data.frame(Variable = character(),  b = numeric(),  SE = numeric(), b_SE = character(), z = numeric(),  P = numeric(),  R2 = numeric(), 
                                     Country = character(), Type = character(), Part = character()  ,stringsAsFactors = FALSE)

Contry5c <-c('China','India','European Country','Nigeria','United States')
parts <- c('pop','poul')
types <- c('All', 'CH')
variables <- c('Seabirds', 'Shorebirds', 'Waterfowls', 'Large wading birds', 'Others')

# 循环每个变量，进行回归分析，并将结果添加到数据框
for(con in Contry5c){#国家
  
      for(typ in types){
        if(typ == 'All'){#全部水鸟
          
          for(part in parts){#人口还是家禽
            
            for (var in variables) {#哪个功能群
              
              con_birdEntropy_df <- subset(Global_birdEntropy_df5_long,Global_birdEntropy_df5_long$Country==con)
              vardata = subset(con_birdEntropy_df,con_birdEntropy_df$Variable==var)
              
              if(part =='pop'){#人口
                model <- lm(`Population(Log)` ~ Value, vardata)
                }
              else{#家禽
                model <- lm(`Poltry(Log)` ~ Value, vardata)
                
              }
              model_summary <- summary(model)
              
              # 提取所需的参数
              b <- model_summary$coefficients[2, 1]
              SE <- model_summary$coefficients[2, 2]
              b_SE <- sprintf("%.2f ± %.2f", b, SE)
              z <- b / SE
              P <- model_summary$coefficients[2, 4]
              R2 <- model_summary$r.squared
              
              # 将结果添加到数据框
              All5_poppoulresults_df <- rbind(All5_poppoulresults_df, data.frame(Variable = var, b = b, SE = SE, b_SE = b_SE, z = z, P = P, R2 = R2, 
                                                                                 Country = con, Type = typ, Part = part))
            }
          }
        }
        
        if(typ == 'CH'){#CH水鸟
          
          for(part in parts){#人口还是家禽
            
            for (var in variables) {#哪个功能群
              
            con_CHEntropy_df <- subset(Global_CHEntropy_df5_long,Global_CHEntropy_df5_long$Country==con)
            vardata = subset(con_CHEntropy_df,con_CHEntropy_df$Variable==var)
            
            if(part =='pop'){#人口
              model <- lm(`Population(Log)` ~ Value, vardata)
            }
            else{#家禽
              model <- lm(`Poltry(Log)` ~ Value, vardata)
              
            }
            model_summary <- summary(model)
            
            # 提取所需的参数
            b <- model_summary$coefficients[2, 1]
            SE <- model_summary$coefficients[2, 2]
            b_SE <- sprintf("%.2f ± %.2f", b, SE)
            z <- b / SE
            P <- model_summary$coefficients[2, 4]
            R2 <- model_summary$r.squared
            
            # 将结果添加到数据框
            All5_poppoulresults_df <- rbind(All5_poppoulresults_df, data.frame(Variable = var, b = b, SE = SE, b_SE = b_SE, z = z, P = P, R2 = R2, 
                                                                               Country = con, Type = typ, Part = part))
            }
            
          }
          
        }
        
      }
}


All5_poppoulresults_df
#fwrite(All5_poppoulresults_df,'/root/autodl-tmp/humPoulResult/data/Mantel_data/All5_poppoulresults_df.csv')





##########绘制中国的情况
head(Global_birdEntropy_df5_long)
head(Global_CHEntropy_df5_long)#CH


Con_China_df_long<-subset(Global_birdEntropy_df5_long,Global_birdEntropy_df5_long$Country=='China')
Con_China_df_long$Variable<-factor(Con_China_df_long$Variable, levels = c('Seabirds', 'Shorebirds', 'Waterfowls', 'Large wading birds', 'Others'))

popall<-ggplot(data = Con_China_df_long, aes(x = Value, y = log(pop), color = Variable)) +
  geom_point(size=2) +  #,shape=21
  scale_color_manual(values = c("#377EB8", "#E41A1C", "#4DAF4A", "#984EA3", "#FF7F00")) +
  geom_smooth(method = "lm",linewidth=0.5, se = T, lty=5,aes(color = Variable)) +  #  , color="black"
  #geom_smooth(method = "lm", se = T, color = "black", linetype = "dashed") +  # 全部点的趋势线，黑色
  #geom_smooth(method = "lm",linewidth=0.5, se = T, lty=5,color = "black") +  #分界
  stat_poly_eq(aes(label = after_stat(paste(eq.label, rr.label,p.value.label, sep = "~~~"))),
               formula = y ~ x,  label.x = Inf, # 右边界
               #label.y=-Inf,
               coef.digits = 3, rr.digits = 2, parse = TRUE, size = 3) +
  #facet_wrap(~ variable) +  #如果要分界
  labs(x = "AE", y = "Population (Log)", title = "All Waterbirds") +
  #ylim(-7,10)+
  theme_bw() +
  theme(legend.position = "right",
        legend.title = element_blank(),
        plot.title = element_text(size = 24, hjust = 0.5),
        axis.title = element_text(size = 20),
        text = element_text(size = 18)
        )
popall


Con_China_CHdf_long<-subset(Global_CHEntropy_df5_long,Global_CHEntropy_df5_long$Country=='China')
Con_China_CHdf_long$Variable<-factor(Con_China_CHdf_long$Variable, levels = c('Seabirds', 'Shorebirds', 'Waterfowls', 'Large wading birds', 'Others'))

popch<-ggplot(data = Con_China_CHdf_long, aes(x = Value, y = `Population(Log)`, color = Variable)) +
  geom_point(size=2) +  #,shape=21
  scale_color_manual(values = c("#377EB8", "#E41A1C", "#4DAF4A", "#984EA3", "#FF7F00")) +
  geom_smooth(method = "lm",linewidth=0.5, se = T, lty=5,aes(color = Variable)) +  #,aes(color = variable)   , color="black"
  #geom_smooth(method = "lm", se = T, color = "black", linetype = "dashed") +  # 全部点的趋势线，黑色
  #geom_smooth(method = "lm",linewidth=0.5, se = T, lty=5,color = "black") +  #分界
  stat_poly_eq(aes(label = after_stat(paste(eq.label, rr.label,p.value.label, sep = "~~~"))),
               formula = y ~ x,  label.x = Inf, # 右边界
               #label.y=-Inf,
               coef.digits = 3, rr.digits = 2, parse = TRUE, size = 3) +
  #facet_wrap(~ variable) +  #如果要分界
  labs(x = "AE", y = "Population (Log)", title = "Confirmed Host") +
  #ylim(-7,10)+
  theme_bw() +
  theme(legend.position = "right",
        legend.title = element_blank(),
        plot.title = element_text(size = 24, hjust = 0.5),
        axis.title = element_text(size = 20),
        text = element_text(size = 18))
popch
library(patchwork)
popall+popch



######饼图，最重要的------------
All5_poppoulresults_df<-fread('/root/autodl-tmp/humPoulResult/data/Mantel_data/All5_poppoulresults_df.csv')
lineardatalong <- All5_poppoulresults_df %>%data.table()
lineardatalong$r <- sqrt(lineardatalong$R2)
lineardatalong<-lineardatalong %>%
  mutate(Pd = cut(P, breaks = c(-Inf, 0.01, 0.05, Inf),
                  labels = c("< 0.01", "0.01 - 0.05", "> 0.05")))%>%
  # mutate(R2d = cut(R2, breaks = c(-Inf, 0.3, 0.4,0.5,0.6, Inf),
  #                  labels = c("< 0.3", "0.3 - 0.4", "0.4 - 0.5", "0.5 - 0.6", "> 0.6"))
  mutate(rd = cut(r, breaks = c(-Inf, 0.4,0.5,0.6,0.7, Inf),
                          labels = c("< 0.4", "0.4 - 0.5", "0.5 - 0.6", "0.6 - 0.7", "> 0.7")))
head(lineardatalong)

# 定义一个函数来筛选每个组合的最大r值对应的行
get_max_r2_rows <- function(dt, type) {
  max_r2 <- dt[Type == type, .(Max_r = max(r)), by = .(Country, Part)]    # 计算每个Country和Part组合的最大R2值
  setkey(dt, Country, Part, r)                # 将最大r值合并回原始数据集
  setkey(max_r2, Country, Part, Max_r)
  max_rows <- max_r2[dt, nomatch=0]      # 筛选出最大R2值对应的行
  #max_rows[, Max_r := NULL]    # Max_r
  return(max_rows)
}

# 分别为ALL和CH类型筛选数据
max_r2_rows_ALL <- get_max_r2_rows(lineardatalong, "All") 
max_r2_rows_CH <- get_max_r2_rows(lineardatalong, "CH")
lineardata_max <- rbind(max_r2_rows_ALL, max_r2_rows_CH)
print(lineardata_max)

# # 使用ifelse语句来创建颜色列
# lineardata_max$Color <- ifelse(lineardata_max$Type == "All" & lineardata_max$rd == "< 0.4", "#F5EAD0",
#                                ifelse(lineardata_max$Type == "All" & lineardata_max$rd == "0.4 - 0.5", "#E3CB8F",
#                                       ifelse(lineardata_max$Type == "All" & lineardata_max$rd == "0.5 - 0.6", "#BF8D44",
#                                              ifelse(lineardata_max$Type == "All" & lineardata_max$rd == "0.6 - 0.7", "#8E581C",
#                                                     ifelse(lineardata_max$Type == "All" & lineardata_max$rd == "> 0.7", "#59310D",
#                                                            ifelse(lineardata_max$Type == "CH" & lineardata_max$rd == "< 0.4", "#FAE3D7",
#                                                                   ifelse(lineardata_max$Type == "CH" & lineardata_max$rd == "0.4 - 0.5", "#F2B396",
#                                                                          ifelse(lineardata_max$Type == "CH" & lineardata_max$rd == "0.5 - 0.6", "#D9705A",
#                                                                                 ifelse(lineardata_max$Type == "CH" & lineardata_max$rd == "0.6 - 0.7", "#B42D34",
#                                                                                        ifelse(lineardata_max$Type == "CH" & lineardata_max$rd == "> 0.7", "#6E0D20",
#                                                                                               NA))))))))))
# 使用ifelse语句来创建颜色列
lineardata_max$Color <- ifelse(lineardata_max$Type == "All" & lineardata_max$rd == "< 0.4", "#FAE3D7",
                               ifelse(lineardata_max$Type == "All" & lineardata_max$rd == "0.4 - 0.5", "#F2B396",
                                      ifelse(lineardata_max$Type == "All" & lineardata_max$rd == "0.5 - 0.6", "#D9705A",
                                             ifelse(lineardata_max$Type == "All" & lineardata_max$rd == "0.6 - 0.7", "#B42D34",
                                                    ifelse(lineardata_max$Type == "All" & lineardata_max$rd == "> 0.7", "#6E0D20",
                                                           ifelse(lineardata_max$Type == "CH" & lineardata_max$rd == "< 0.4", "#FAE3D7",
                                                                  ifelse(lineardata_max$Type == "CH" & lineardata_max$rd == "0.4 - 0.5", "#F2B396",
                                                                         ifelse(lineardata_max$Type == "CH" & lineardata_max$rd == "0.5 - 0.6", "#D9705A",
                                                                                ifelse(lineardata_max$Type == "CH" & lineardata_max$rd == "0.6 - 0.7", "#B42D34",
                                                                                       ifelse(lineardata_max$Type == "CH" & lineardata_max$rd == "> 0.7", "#6E0D20",
                                                                                              NA))))))))))


lineardata_max
#fwrite(lineardata_max,"/root/autodl-tmp/humPoulResult/data/Mantel_data/lineardata_max.csv")
######——制作图例---------
colors_All <- c("< 0.4"="#FAE3D7", "0.4 - 0.5" = "#F2B396", "0.5 - 0.6" = "#D9705A", "0.6 - 0.7" = "#B42D34", "> 0.7" = "#6E0D20")
#colors_CH <- c("< 0.3"="#E7F2F1", "0.3 - 0.4" = "#AFDCD5", "0.4 - 0.5" = "#6DB2AD", "0.5 - 0.6" = "#347C78", "> 0.6" = "#174A41")
colors_CH <- c("< 0.4"="#FAE3D7", "0.4 - 0.5" = "#F2B396", "0.5 - 0.6" = "#D9705A", "0.6 - 0.7" = "#B42D34", "> 0.7" = "#6E0D20")

data_All <- lineardata_max[lineardata_max$Type == "All", ]
data_CH <- lineardata_max[lineardata_max$Type == "CH", ]
p1<-ggplot() +
  geom_bar(data=data_All, aes(x=Country, y=Max_r, fill=rd), stat="identity", position=position_dodge())+
  scale_fill_manual(name="r (All)", values=colors_All, breaks=names(colors_All), labels=names(colors_All))+
  theme(text = element_text(size=20),
        legend.position = "bottom",
        legend.background = element_rect(fill = "transparent")) # 透明的图例背景
#legend.box.background = element_rect(fill = "transparent")) # 透明的图例面板)
p2<-ggplot() +
  geom_bar(data=data_CH, aes(x=Country, y=Max_r, fill=rd), stat="identity", position=position_dodge())+
  scale_fill_manual(name="r (CH)", values=colors_CH, breaks=names(colors_CH), labels=names(colors_CH))+
  theme(text = element_text(size=20),
        legend.position = "bottom",
        legend.background = element_rect(fill = "transparent")) # 透明的图例背景
#legend.box.background = element_rect(fill = "transparent")) # 透明的图例面板)

p1+p2
#ggsave("/root/autodl-tmp/zyresult/Entropy/poppulresult/piedata_Lengend.png", p1+p2, bg = "transparent", width = 18, height = 8, units = "in")



colors_gradient <- c("#FAE3D7", "#F2B396", "#D9705A", "#B42D34", "#6E0D20")
library(scales)
# 创建一个连续的图例
ggplot(data_All, aes(x=Country, y=Max_r, fill=Max_r)) +
  geom_bar(stat="identity", position=position_dodge()) +
  scale_fill_gradientn(colors=colors_gradient, 
                       values=rescale(c(0.4, 0.5, 0.6, 0.7, 0.8)),
                       name="r (All)",
                       limits=c(0.4, 0.8),
                       breaks=c(0.4, 0.5, 0.6, 0.7, 0.8),
                       labels=c("< 0.4", "0.5", "0.6", "0.7", "> 0.8")) +
  theme(text = element_text(size=10),
        legend.position = "bottom",
        legend.background = element_rect(fill="transparent"),
        legend.title = element_blank())


######——制作底图---------
#世界地图
world.map <- rnaturalearth::ne_countries(returnclass = "sf") %>% filter(continent != "Antarctica")

Con_popentrpoul_sf <- st_read("/root/autodl-tmp/zyresult/Con_popentrpoul_sf_EU.shp");head(Con_popentrpoul_sf)
Con_China<-subset(Con_popentrpoul_sf,name_ec=="China")
Con_India<-subset(Con_popentrpoul_sf,name_ec=="India") 
Con_EU<-subset(Con_popentrpoul_sf,name_ec=="European Country")
Con_Nigeria<-subset(Con_popentrpoul_sf,name_ec=="Nigeria")
Con_US<-subset(Con_popentrpoul_sf,name_ec=="United States") 

ggplot() +
  geom_sf(data = world.map, fill = "lightgrey", color = NA) + 
  geom_sf(data = Con_China, fill = "#baa1a1", color = NA) +    #lightcoral
  geom_sf(data = Con_India, fill = "#baa1a1", color = NA) + 
  geom_sf(data = Con_EU, fill = "#baa1a1", color = NA) + 
  geom_sf(data = Con_Nigeria, fill = "#baa1a1", color = NA) + 
  geom_sf(data = Con_US, fill = "#baa1a1", color = NA) + 
  coord_sf(crs = "+proj=longlat +datum=WGS84", xlim = c(-160, 165), ylim = c(-56, 90)) +
  theme_bw() + # 使用简洁主题
  theme(legend.position = "none",
        axis.text = element_text(size=12),
        panel.grid.major = element_blank() ) # 删除主要的灰色虚线)



######【China饼图】------
piedata_China <- lineardata_max[lineardata_max$Country == "China", ]
piedata_China$Value<-c(1,1,1,1)
color_China<-c("pop_All"="#6E0D20","pop_CH"="#B42D34","poul_CH"="#B42D34","poul_All"="#6E0D20")

# 顺序
desired_order <- c("pop_CH", "pop_All", "poul_All", "poul_CH")
piedata_China$Order <- factor(paste(piedata_China$Part, piedata_China$Type, sep = "_"), levels = desired_order)
piedata_China <- piedata_China[order(piedata_China$Order), ]
piedata_China$Color <- factor(piedata_China$Color, levels = unique(piedata_China$Color))
#c("#8E581C","#D9705A","#B42D34","#59310D")

p <- ggplot(piedata_China, aes(x = "", y = Value, fill = Order)) +
  geom_bar(stat = "identity", width = 1) +
  coord_polar("y", start = 0) +
  scale_fill_manual(values = color_China) +
  theme_void() +
  theme(
    panel.background = element_rect(fill = "transparent"), # 透明的面板背景
    plot.background = element_rect(fill = "transparent", color = NA), # 透明的绘图背景
    legend.position = "none"
  )
p
# 保存为 PNG 图像，背景透明
ggsave("/root/autodl-tmp/zyresult/Entropy/poppulresult/piedata_China.png", p, bg = "transparent", width = 10, height = 8, units = "in")


######【India 饼图】------
piedata_India <- lineardata_max[lineardata_max$Country == "India", ]
piedata_India$Value<-c(1,1,1,1)

# 顺序
desired_order <- c("pop_CH", "pop_All", "poul_All", "poul_CH")
piedata_India$Order <- factor(paste(piedata_India$Part, piedata_India$Type, sep = "_"), levels = desired_order)
piedata_India <- piedata_India[order(piedata_India$Order), ]
piedata_India$Color <- factor(piedata_India$Color, levels = unique(piedata_India$Color))


p <- ggplot(piedata_India, aes(x = "", y = Value, fill = Color)) +
  geom_bar(stat = "identity", width = 1) +
  coord_polar("y", start = 0) +
  scale_fill_manual(values = levels(piedata_India$Color)) +
  theme_void() +
  theme(
    panel.background = element_rect(fill = "transparent"), # 透明的面板背景
    plot.background = element_rect(fill = "transparent", color = NA), # 透明的绘图背景
    legend.position = "none"
  )
p
# 保存为 PNG 图像，背景透明
ggsave("/root/autodl-tmp/zyresult/Entropy/poppulresult/piedata_India.png", p, bg = "transparent", width = 10, height = 8, units = "in")


######【EU 饼图】------
piedata_EU <- lineardata_max[lineardata_max$Country == "European Country", ]
piedata_EU$Value<-c(1,1,1,1)
color_EU<-c("pop_All"="#D9705A","pop_CH"="#F2B396","poul_CH"="#F2B396","poul_All"="#F2B396")  #8E581C
# 顺序
desired_order <- c("pop_CH", "pop_All", "poul_All", "poul_CH")
piedata_EU$Order <- factor(paste(piedata_EU$Part, piedata_EU$Type, sep = "_"), levels = desired_order)
piedata_EU <- piedata_EU[order(piedata_EU$Order), ]
piedata_EU$Color <- factor(piedata_EU$Color, levels = unique(piedata_EU$Color))


p <- ggplot(piedata_EU, aes(x = "", y = Value, fill = Order)) +
  geom_bar(stat = "identity", width = 1) +
  coord_polar("y", start = 0) +
  scale_fill_manual(values = color_EU) +
  theme_void() +
  theme(
    panel.background = element_rect(fill = "transparent"), # 透明的面板背景
    plot.background = element_rect(fill = "transparent", color = NA), # 透明的绘图背景
    legend.position = "none"
  )
p
# 保存为 PNG 图像，背景透明
ggsave("/root/autodl-tmp/zyresult/Entropy/poppulresult/piedata_EU.png", p, bg = "transparent", width = 10, height = 8, units = "in")

######【Nigeria 饼图】------
piedata_Nigeria <- lineardata_max[lineardata_max$Country == "Nigeria", ]
piedata_Nigeria$Value<-c(1,1,1,1)

color_Nigeria<-c("pop_All"="#FAE3D7","pop_CH"="#FAE3D7","poul_CH"="#FAE3D7","poul_All"="#FAE3D7")

# 顺序
desired_order <- c("pop_CH", "pop_All", "poul_All", "poul_CH")
piedata_Nigeria$Order <- factor(paste(piedata_Nigeria$Part, piedata_Nigeria$Type, sep = "_"), levels = desired_order)
piedata_Nigeria <- piedata_Nigeria[order(piedata_Nigeria$Order), ]
piedata_Nigeria$Color <- factor(piedata_Nigeria$Color, levels = unique(piedata_Nigeria$Color))


p <- ggplot(piedata_Nigeria, aes(x = "", y = Value, fill = Order)) +
  geom_bar(stat = "identity", width = 1) +
  coord_polar("y", start = 0) +
  scale_fill_manual(values = color_Nigeria) +
  theme_void() +
  theme(
    panel.background = element_rect(fill = "transparent"), # 透明的面板背景
    plot.background = element_rect(fill = "transparent", color = NA), # 透明的绘图背景
    legend.position = "none"
  )
p
# 保存为 PNG 图像，背景透明
ggsave("/root/autodl-tmp/zyresult/Entropy/poppulresult/piedata_Nigeria.png", p, bg = "transparent", width = 10, height = 8, units = "in")


######【US 饼图】------
piedata_US <- lineardata_max[lineardata_max$Country == "United States", ]
piedata_US$Value<-c(1,1,1,1)

color_US<-c("pop_All"="#B42D34","pop_CH"="#D9705A","poul_CH"="#D9705A","poul_All"="#D9705A")

# 顺序
desired_order <- c("pop_CH", "pop_All", "poul_All", "poul_CH")
piedata_US$Order <- factor(paste(piedata_US$Part, piedata_US$Type, sep = "_"), levels = desired_order)
piedata_US <- piedata_US[order(piedata_US$Order), ]
piedata_US$Color <- factor(piedata_US$Color, levels = unique(piedata_US$Color))


p <- ggplot(piedata_US, aes(x = "", y = Value, fill = Order)) +
  geom_bar(stat = "identity", width = 1) +
  coord_polar("y", start = 0) +
  scale_fill_manual(values = color_US) +
  theme_void() +
  theme(
    panel.background = element_rect(fill = "transparent"), # 透明的面板背景
    plot.background = element_rect(fill = "transparent", color = NA), # 透明的绘图背景
    legend.position = "none"
  )
p
# 保存为 PNG 图像，背景透明
ggsave("/root/autodl-tmp/zyresult/Entropy/poppulresult/piedata_US.png", p, bg = "transparent", width = 10, height = 8, units = "in")










#b ： 弦图############

######1.全部 -----------
Global_birdEntropy_dfna<-Global_birdEntropy_df %>% na.omit()
summary(Global_birdEntropy_dfna)#ALL
# 比较五列的值，并创建新列Guide
Global_birdEntropy_dfna <- Global_birdEntropy_dfna %>%
  mutate(Guide = pmax(Seabirds, Shorebirds, Waterfowls, `Large wading birds`, Others, na.rm = TRUE)) %>%
  mutate(Guide = case_when(
    Guide == Seabirds ~ "Seabirds",
    Guide == Shorebirds ~ "Shorebirds",
    Guide == Waterfowls ~ "Waterfowls",
    Guide == `Large wading birds` ~ "Large wading birds",
    Guide == Others ~ "Others"
  ))

# 统计每个国家的Guide类别的像元数量
country_guide_count <- Global_birdEntropy_dfna %>%
  count(Country, Guide)
country_guide_count <- country_guide_count %>%
  group_by(Guide) %>%
  mutate(percentage = n / sum(n))
country_guide_count <- country_guide_count %>%
  group_by(Country) %>%
  mutate(percentage_c = n / sum(n))

# 使用弦图可视化
library(circlize)

# 准备数据
chord_data <- country_guide_count %>%
  filter(n > 0) %>%
  select(Guide, Country, n)%>%
  mutate(Country = case_when(
    Country == "European Country" ~ "EU",
    Country == "United States" ~ "US",
    TRUE ~ Country ))


# 生成弦图
#chordDiagram(chord_data, transparency = 0.5)


#设定颜色
grid.col = NULL
grid.col[c("China", "India", "EU", "Nigeria","US")] = c("#E3CB8F", "#E3CB8F","#E3CB8F",  "#E3CB8F", "#E3CB8F")#c("#FF6699","#33CCFF","#99FFC9","#E600E6","#FFD700")
grid.col[c("Seabirds", "Shorebirds", "Waterfowls", "Large wading birds","Others")] = c("#377EB8","#E41A1C","#4DAF4A","#984EA3","#f7a156") #c("#E41A1C","#4DAF4A","#377EB8","#984EA3","#f7a156") 
#grid.col[rownames(chord_data)] = c("#E41A1C","#377EB8","#4DAF4A","#984EA3","#f7a156")  #c("#FF6699","#33CCFF","#99FFC9","#E600E6","#FFD700")


#第一种
par(cex = 1.5) # 这里的1.5是字体大小的倍数
chordDiagram(chord_data,
             #annotationTrack = "grid", # 不创建表示标签的轨迹
             big.gap = 10,
             diffHeight = 0.06, # 外圈与连线间隔高度
             grid.col = grid.col, # 线条颜色
             link.lwd = 0.02, # 线条宽度
             transparency = 0.5) # 连接颜色透明度

circos.clear()

#第二种
par(cex = 2) # 这里的1.5是字体大小的倍数
chordDiagram(chord_data, big.gap = 10, grid.col = grid.col, annotationTrack = "grid", preAllocateTracks = list(track.height = 0.1))

# 在每个扇区添加横向文本
circos.track(track.index = 1, panel.fun = function(x, y) {
  circos.text(CELL_META$xcenter, CELL_META$ylim[1] + 0.1, 
              CELL_META$sector.index, facing = "clockwise", 
              niceFacing = TRUE, adj = c(0, 0.5))
}, bg.border = NA)

# 清除图层，以便后续操作
circos.clear()









######2.宿主----------
Global_CHEntropy_dfna<-Global_CHEntropy_df%>% na.omit()
summary(Global_CHEntropy_dfna)#CH

# 比较五列的值，并创建新列Guide
Global_CHEntropy_dfna <- Global_CHEntropy_dfna %>%
  mutate(Guide = pmax(Seabirds, Shorebirds, Waterfowls, `Large wading birds`, Others, na.rm = TRUE)) %>%
  mutate(Guide = case_when(
    Guide == Seabirds ~ "Seabirds",
    Guide == Shorebirds ~ "Shorebirds",
    Guide == Waterfowls ~ "Waterfowls",
    Guide == `Large wading birds` ~ "Large wading birds",
    Guide == Others ~ "Others"
  ))

# 统计每个国家的Guide类别的像元数量
country_guide_count <- Global_CHEntropy_dfna %>%
  count(Country, Guide)
country_guide_count <- country_guide_count %>%
  group_by(Guide) %>%
  mutate(percentage = n / sum(n))
country_guide_count <- country_guide_count %>%
  group_by(Country) %>%
  mutate(percentage_c = n / sum(n))
# 使用弦图可视化
# 使用弦图可视化
library(circlize)

# 准备数据
chord_data <- country_guide_count %>%
  filter(n > 0) %>%
  select(Guide, Country, n)%>%
  mutate(Country = case_when(
    Country == "European Country" ~ "EC",
    Country == "United States" ~ "US",
    TRUE ~ Country ))
    

# 生成弦图
#chordDiagram(chord_data, transparency = 0.5)


#设定颜色
grid.col = NULL
grid.col[c("Seabirds", "Shorebirds", "Waterfowls", "Large wading birds","Others")] = c("#377EB8","#E41A1C","#4DAF4A","#984EA3","#f7a156") #c("#E41A1C","#4DAF4A","#377EB8","#984EA3","#f7a156") 
#grid.col[rownames(chord_data)] = c("#E41A1C","#377EB8","#4DAF4A","#984EA3","#f7a156")  #c("#FF6699","#33CCFF","#99FFC9","#E600E6","#FFD700")
grid.col[c("China", "India", "EC", "Nigeria","US")] = c("#D9705A", "#D9705A","#D9705A",  "#D9705A", "#D9705A")  #c("#FF6699","#33CCFF","#99FFC9","#E600E6","#FFD700")


#第一种
par(cex = 1.5) # 这里的1.5是字体大小的倍数
chordDiagram(chord_data,
             #annotationTrack = "grid", # 不创建表示标签的轨迹
             big.gap = 10,
             diffHeight = 0.06, # 外圈与连线间隔高度
             grid.col = grid.col, # 线条颜色
             link.lwd = 0.02, # 线条宽度
             transparency = 0.5) # 连接颜色透明度

circos.clear()

#第二种
par(cex = 2) # 这里的1.5是字体大小的倍数
chordDiagram(chord_data, big.gap = 10, grid.col = grid.col, annotationTrack = "grid", preAllocateTracks = list(track.height = 0.1))

# 在每个扇区添加横向文本
circos.track(track.index = 1, panel.fun = function(x, y) {
  circos.text(CELL_META$xcenter, CELL_META$ylim[1] + 0.1, 
              CELL_META$sector.index, facing = "clockwise", 
              niceFacing = TRUE, adj = c(0, 0.5))
}, bg.border = NA)

# 清除图层，以便后续操作
circos.clear()









#####计算lineardatalong#####
#####1.China####
######【人口】散点图并添加线性回归拟合线---------
#log不分面
popAll<-ggplot(data = Con_China_df7_long2, aes(x = value, y = log(pop), color = variable)) +
  geom_point(size=1) +  #,shape=21
  scale_color_manual(values = c("#377EB8", "#E41A1C", "#4DAF4A", "#984EA3", "#FF7F00")) +
  geom_smooth(method = "lm",linewidth=0.5, se = T, lty=5,aes(color = variable)) +  #,aes(color = variable)   , color="black"
  #geom_smooth(method = "lm", se = T, color = "black", linetype = "dashed") +  # 全部点的趋势线，黑色
  #geom_smooth(method = "lm",linewidth=0.5, se = T, lty=5,color = "black") +  #分界
  # stat_poly_eq(aes(label = after_stat(paste(eq.label, rr.label,p.value.label, sep = "~~~"))),
  #              formula = y ~ x,  label.x = Inf, # 右边界
  #              label.y=-Inf,
  #              coef.digits = 3, rr.digits = 2, parse = TRUE, size = 3) +
  #facet_wrap(~ variable) +  #如果要分界
  labs(x = "AE", y = "Population (Log)", title = "All Waterbirds") +
  #ylim(-7,10)+
  theme_bw() +
  theme(legend.position = "none",
        legend.title = element_blank(),
        plot.title = element_text(size = 24, hjust = 0.5),
        axis.title = element_text(size = 20),
        text = element_text(size = 18))
popAll
#log分面
ggplot(data = Con_China_df7_long2, aes(x = value, y = log(pop), color = variable)) +
  geom_point(size=1) +  #,shape=21
  scale_color_manual(values = c("#377EB8", "#E41A1C", "#4DAF4A", "#984EA3", "#FF7F00")) +
  geom_smooth(method = "lm",linewidth=0.5, se = T, lty=5, color="black") +  #,aes(color = variable)   , color="black"
  #geom_smooth(method = "lm", se = T, color = "black", linetype = "dashed") +  # 全部点的趋势线，黑色
  #geom_smooth(method = "lm",linewidth=0.5, se = T, lty=5,color = "black") +  #分界
  stat_poly_eq(aes(label = after_stat(paste(eq.label, rr.label, sep = "~~~"))),  #,p.value.label
               formula = y ~ x,  label.x = Inf, # 右边界
               label.y=-Inf,
               coef.digits = 3, rr.digits = 2, parse = TRUE, size = 3) +
  facet_wrap(~ variable) +  #如果要分界
  labs(x = "Waterbird Activity Entropy", y = "Population (Log)", title = "All") +
  #ylim(-7,10)+
  theme_bw() +
  theme(legend.position = "bottom",
        legend.title = element_blank(),
        plot.title = element_text(size = 24, hjust = 0.5),
        axis.title = element_text(size = 20),
        text = element_text(size = 18))

#计算参数
head(Con_China_df7_long2)
# 初始化一个空的数据框来存储结果
popresults_df <- data.frame(  Variable = character(),  b = numeric(),  SE = numeric(),  z = numeric(),  P = numeric(),  R2 = numeric(),  stringsAsFactors = FALSE)
variables <- c('Seabirds', 'Shorebirds', 'Waterfowls', 'Large wading bird', 'Others')

# 循环每个变量，进行回归分析，并将结果添加到数据框
for (var in variables) {
  model <- lm(log(pop) ~ value, data = Con_China_df7_long2[Con_China_df7_long2$variable == var, ])
  model_summary <- summary(model)
  
  # 提取所需的参数
  b <- model_summary$coefficients[2, 1]
  SE <- model_summary$coefficients[2, 2]
  z <- b / SE
  P <- model_summary$coefficients[2, 4]
  R2 <- model_summary$r.squared
  
  # 将结果添加到数据框
  popresults_df <- rbind(popresults_df, data.frame(Variable = var, b_pop = b, SE_pop = SE, z_pop = z, P_pop = P, R2_pop = R2))
}
head(popresults_df)



###### 【家禽】散点图并添加线性回归拟合线-----
pouAll<-ggplot(data = Con_China_df7_long2, aes(x = value, y = log(poul), color = variable)) +
  geom_point(size=1) +  #,shape=21
  scale_color_manual(values = c("#377EB8", "#E41A1C", "#4DAF4A", "#984EA3", "#FF7F00")) +
  geom_smooth(method = "lm",linewidth=0.5, se = T, lty=5,aes(color = variable)) +  #,aes(color = variable)   , color="black"
  #geom_smooth(method = "lm", se = T, color = "black", linetype = "dashed") +  # 全部点的趋势线，黑色
  #geom_smooth(method = "lm",linewidth=0.5, se = T, lty=5,color = "black") +  #分界
  # stat_poly_eq(aes(label = after_stat(paste(eq.label, rr.label,p.value.label, sep = "~~~"))),
  #              formula = y ~ x,  label.x = Inf, # 右边界
  #              label.y=-Inf,
  #              coef.digits = 3, rr.digits = 2, parse = TRUE, size = 3) +
  # facet_wrap(~ variable) +  #如果要分界
  labs(x = "AE", y = "Poultry (Log)", title = "All Waterbirds") +
  #ylim(-7,10)+
  theme_bw() +
  theme(legend.position = "none",
        legend.title = element_blank(),
        plot.title = element_text(size = 24, hjust = 0.5),
        axis.title = element_text(size = 20),
        text = element_text(size = 18))
pouAll

ggplot(data = Con_China_df7_long2, aes(x = value, y = log(poul), color = variable)) +
  geom_point(size=1) +  #,shape=21
  scale_color_manual(values = c("#377EB8", "#E41A1C", "#4DAF4A", "#984EA3", "#FF7F00")) +
  geom_smooth(method = "lm",linewidth=0.5, se = T, lty=5, color="black") +  #,aes(color = variable)   , color="black"
  #geom_smooth(method = "lm", se = T, color = "black", linetype = "dashed") +  # 全部点的趋势线，黑色
  #geom_smooth(method = "lm",linewidth=0.5, se = T, lty=5,color = "black") +  #分界
  stat_poly_eq(aes(label = after_stat(paste(eq.label, rr.label, sep = "~~~"))), #,p.value.label
               formula = y ~ x,  label.x = Inf, # 右边界
               label.y=-Inf,
               coef.digits = 3, rr.digits = 2, parse = TRUE, size = 3) +
  facet_wrap(~ variable) +  #如果要分界
  labs(x = "Waterbird Activity Entropy", y = "Poultry (Log)", title = "All") +
  #ylim(-7,10)+
  theme_bw() +
  theme(legend.position = "bottom",
        legend.title = element_blank(),
        plot.title = element_text(size = 24, hjust = 0.5),
        axis.title = element_text(size = 20),
        text = element_text(size = 18))


#计算参数
head(Con_China_df7_long2)
# 初始化一个空的数据框来存储结果
poulresults_df <- data.frame(  Variable = character(),  b = numeric(),  SE = numeric(),  z = numeric(),  P = numeric(),  R2 = numeric(),  stringsAsFactors = FALSE)
variables <- c('Seabirds', 'Shorebirds', 'Waterfowls', 'Large wading bird', 'Others')

# 循环每个变量，进行回归分析，并将结果添加到数据框
for (var in variables) {
  model <- lm(log(poul) ~ value, data = Con_China_df7_long2[Con_China_df7_long2$variable == var, ])
  model_summary <- summary(model)
  # 提取所需的参数
  b <- model_summary$coefficients[2, 1]
  SE <- model_summary$coefficients[2, 2]
  z <- b / SE
  P <- model_summary$coefficients[2, 4]
  R2 <- model_summary$r.squared
  poulresults_df <- rbind(poulresults_df, data.frame(Variable = var, b_poul = b, SE_poul = SE, z_poul = z, P_poul = P, R2_poul = R2))
}
head(poulresults_df)

#合并人口和家禽的线性拟合结果
China_poppoul_r<-left_join(popresults_df, poulresults_df, by="Variable")
#fwrite(China_poppoul_r,"/root/autodl-tmp/zyresult/Entropy/poppulresult/China_poppoul_r.csv")







