
basePath<-"/root/autodl-tmp/humPoulResult/data/"
crs <- '+proj=longlat +datum=WGS84'
world.map <- rnaturalearth::ne_countries(returnclass = "sf") |>filter(continent != "Antarctica")
globalCountry <- vect(world.map) 
globalRaster <- rast(vals=1:259200,nrows=360, ncols=720,xmin=-180, xmax=180,ymin=-90, ymax=90,crs=crs)
coast <- rnaturalearth::ne_coastline(scale = "small", returnclass = "sf")

allDf <- fread(paste0(basePath,'allDf786_reclass.csv'))
speciesPixelNumPath <- list.files('/root/autodl-tmp/humPoulResult/data/single_model',pattern = '.tif',full.names = T)
spName <- basename(speciesPixelNumPath) %>% str_sub(.,1,-5)                      
speciesPixelNumPath2 <- speciesPixelNumPath[spName%in%allDf$LatName]

`%notin%` <- Negate(`%in%`)

#####Fig1.CV of Species Richness & AE###############
#######a-----------
#CV of Species Richness calculation

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
    panel.grid = element_blank(),
    axis.ticks.length = unit(-4, "pt") #坐标轴ticks内向
  )

pSRCV+pSRCV2

#lit----
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
  xlab('AE')+
  xlim(c(0,5.5))+
  ylim(c(-56,90))+
  theme_bw()+
  # labs(y="Latitude(°)",)+
  #geom_ribbon(aes(y=lat, xmin=row_means-row_sds,xmax = row_means+row_sds),fill = "lightgrey", alpha=0.5)+
  coord_fixed(ratio = 0.105)+ #调整y轴/x轴的比例
  # theme_classic()+
  theme(
    axis.title.y = element_blank(),
    axis.text.y=element_blank(),
    panel.grid = element_blank(),
    axis.ticks.length = unit(-4, "pt"), #坐标轴ticks内向
    legend.position = "none"
  )
pAE2

pAE+pAE2


#lit--------
 #
Vborder <- vect('/root/autodl-tmp/Wallace_zoogeographic/newValisBorder.shp') 
vRaster <- rasterize(Vborder,AE,field='name')
valisEntropy <- c(AE,vRaster) %>% terra::as.data.frame() %>%na.omit()


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
  labs(y = "AE", x = " ")











#######Fig2. AE validation#########################
#a##########
library(terra)
library(data.table)
library(dplyr)
outBreak <- fread('/root/autodl-tmp/YANZHENG/point/allData.csv')
outBreak$Longitude <- as.numeric(outBreak$Longitude)
addGeom <- cellFromXY(globalRaster,outBreak[,c('Longitude','Latitude')]) %>% 
  cbind(id=.,outBreak)
thinData <- unique(data.table(addGeom),by='id') %>% dplyr::select(.,-id)
outBreak2 <- vect(thinData,geom=c('Longitude','Latitude'),crs=crs)
plot(outBreak2)
overEntropy <- rast(paste0(basePath,'AE_data/AE.tif'))%>% mask(globalCountry)
global(overEntropy,quantile,probs=seq(0, 1, 0.05),na.rm=T)

plotEntropy <- ifel(overEntropy>=4.107,1,0)
#加载海岸线数据
coast <- rnaturalearth::ne_coastline(scale = "small", returnclass = "sf")
crs <- '+proj=longlat +datum=WGS84'
ggplot() +
  geom_spatraster(data = plotEntropy) +
  geom_spatvector(data=coast,fill=NA)+
  coord_sf(crs = crs,xlim=c(-160,180),ylim=c(-56,90))+
  geom_spatvector(data=outBreak2,size=2, shape = 4, aes(color='Outbreak of\navian influenza'))+
  theme_bw()+
  scale_color_manual(values='#fc4e4e', labels='Outbreak of\navian influenza') +
  scale_fill_gradient(low = "#132B43",high = "yellow" ,space = "Lab",n.break=2,
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

#b##########
library(lubridate)
outBreak <- fread('/root/autodl-tmp/YANZHENG/point/allData.csv') %>% subset(Diagnosis.status=='Confirmed'&Animal.type%in%c('Domestic','Wild'))
outBreak1 <- outBreak

# r <- subset(outBreak,outBreak$Species=='Unspecified Bird'&outBreak$Animal.type=='	
# Domestic')
outBreak1$label <- str_extract(outBreak1$Serotype,'HPAI|LPAI')
outBreak1$h_label <- str_extract(outBreak1$Serotype,'H[0-9]N[0-9]|H[0-9]')
# r <- outBreak[outBreak$Longitude>60&outBreak$Longitude<120&outBreak$Latitude>35&outBreak$Latitude<60,]
outBreak1 <- subset(outBreak1,outBreak1$h_label%notin%c('H9N2','H5N6'))

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
     cex.main=1.5,  # 主标题字体放大1.5倍
     cex.sub=1.5,   # 副标题字体放大1.5倍
     cex.axis=1.5,  # 坐标轴刻度字体放大1.5倍
     cex.lab=1.5)#网格线设置



threshold <- 4.107

TP <- sum(df$sum>=threshold,na.rm = T)  ;TP
FN <- sum(df$sum<threshold,na.rm = T)   ;FN
FP <- sum(allDf$sum>=threshold,na.rm = T)-TP  ;FP
TN <- nrow(allDf)-TP-FN-FP  ;TN

accuracy <- (TP+TN)/(TP+TN+FN+FP)  ;accuracy
precise <- (TP)/(TP+FP)  ;precise
recall <- (TP)/(TP+FN)  ;recall
FI <- 2*((recall*precise)/(recall+precise))  ;FI


#c###########
#health expenditure
healthData <- fread('/root/autodl-tmp/humPoulResult/data/health expenditure.csv',header = T,drop=c(5:44,66:69))
globalSHP <- vect('/root/autodl-tmp/worldBorder/ne_10m_admin_0_countries.shp')
result_df2<-fread("/root/autodl-tmp/result/validationCountry.csv")
countryData <- as.data.frame(globalSHP) %>% .[,c('NAME_LONG','POP_EST','POP_RANK','ECONOMY','INCOME_GRP','CONTINENT','GDP_MD')]

result_country_tmp <- fread('/root/autodl-tmp/humPoulResult/data/result_country_tmp.csv')
result_df3 <- left_join(result_country_tmp,countryData,by=c('country'='NAME_LONG'))
result_df4 <- left_join(result_df3,healthData,by=c('country'='Country Name'))
result_df4$health <- rowMeans(result_df4[,21:41])
# result_df3$tt <- result_df3$TP+result_df3$FN
# result_df4 <- subset(result_df3,result_df3$tt>10)
# result_df4$health <- ifelse(is.na(result_df4$health),0,result_df4$health)
# result_df4$health_log <- log(result_df4$health)

a <- subset(result_df4,log(result_df4$health)<4.5&result_df4$accuracy>0.75)
b <- select(a,c(1:11,41:42))
result_df4$per <- result_df4$df/result_df4$alldf

ggplot(data=subset(result_df4,result_df4$alldf>10&result_df4$country!='Guatemala'))+
  geom_point(aes(log(health),accuracy,color=CONTINENT,size=POP_EST,alpha=0.5))+
  geom_smooth(aes(log(health),accuracy),method = "lm",formula = )+
  xlab('health')+
  ylab('Accuracy')+
  theme_bw()+
  theme(
    panel.background = element_blank()
  )
# 计算POP_EST的十个分位数
pop_quantiles <- quantile(result_df3$POP_EST, probs=seq(0, 1, by=0.1), na.rm = TRUE)

# 将POP_EST映射到1到10的大小
result_df3$POP_EST_size <- findInterval(result_df3$POP_EST, vec = pop_quantiles)
result_df3$POP_EST <- as.numeric(result_df3$POP_EST)
result_df3 <- mutate(result_df3, POP_EST_label = cut(POP_EST, breaks = c(35000, 250000, 4970000, 7800000, 10500000, 14700000, 22800000, 34800000, 49700000, 103500000, 1397800000),
                                               labels = c("3.5 - 25", "25 - 497", "497 - 780", "780 - 1050", "1050 - 1470", "1470-2280", "2280 - 3480", "3480 - 4970", "4970 - 10350", "10350 - 139780")))
# 创建一个向量，包含每个分位数区间的标签
# size_labels <- sapply(2:11, function(x) {
#   paste(format(round(pop_quantiles[x-1]/10000)), "-", format(round(pop_quantiles[x]/10000)))
# })

# 按大洲和人口排序
result_df3 <- result_df3[order(result_df3$CONTINENT, -result_df3$POP_EST), ]

# 选择每个大洲人口排名前三的国家
top_countries_by_continent <- result_df3 %>%
  group_by(CONTINENT) %>%
  slice_max(POP_EST, n = 3) %>%
  ungroup() %>%
  dplyr::select(country)
top_countries_by_continent<-c(top_countries_by_continent)

result_df4 <- result_df3[result_df3$country %in% top_countries_by_continent$country[-11], ]   #[-11]为了去掉Canada
# 使用 lm() 函数进行线性拟合
linear_model <- lm(accuracy ~ log(GDP_MD), data = result_df3)
summary(linear_model)
# 绘制图形
size_labels <- c("3.5 - 25", "25 - 497", "497 - 780", "780 - 1050", "1050 - 1470", "1470-2280", "2280 - 3480", "3480 - 4970", "4970 - 10350", "10350 - 139780")
my_formula <- y ~ x
ggplot(result_df3, aes(log(GDP_MD), accuracy, color=CONTINENT)) +
  geom_point(aes(size=POP_EST_size), alpha=0.5) +
  geom_text(data=result_df4, aes(label=country, x=log(GDP_MD), y=accuracy),size=4.5, hjust=0.5, vjust=2) +
  scale_color_manual(values = c("#984EA3", "#E41A1C","#4DAF4A",  "#FF7F00", "#a38900","#377EB8"))+
  geom_smooth(method = "gam", se = T, fill="lightgrey",color = "grey", alpha=0.3, linetype = "dashed") +  # 全部点的趋势线，黑色  #loess  gam

  xlab('GDP (Log)') +
  ylab('Accuracy') +
  theme_bw() +
  theme(panel.background = element_blank(),
        text = element_text(size=18)) +
  scale_size(range = c(1, 10), breaks=1:10, labels=size_labels) +
  guides(size=guide_legend(title=expression(paste("Human population (" ~ 10^4 ~")"))),
         color=guide_legend(title="Continent"))


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
    legend.position="none"
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
globalSHP2 <- terra::aggregate(globalSHP,'name_r')


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
  hotResult2 <- data.table(countryName=globalSHP2$name_r,hotSum=hotResult$popd2015_30,threshold=i-1)   #threshold=(i-1)*10
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
quan_entr<- 4.107
quan_poul <- quantiles_poul[1,84]
#重分类
pop_new <- ifel(popd2015>quan_pop,1,0)
entr_new <- ifel(Entropy>quan_entr,1,0)
poul_new <- ifel(poul2015>quan_poul,1,0)

hot60Area <- pop_new+entr_new+poul_new
hot60Area2 <- ifel(hot60Area==3,1,0)

plot(hot60Area2)
#writeRaster(hot60Area2,'/root/autodl-tmp/humPoulResult/data/Hot_data/Pop50poulentr.tif')
library(rnaturalearth)
coast <- ne_coastline(scale = "small", returnclass = "sf") %>% vect()
crs <- '+proj=longlat +datum=WGS84'
hot60Area2_df<-as.data.frame(hot60Area2,xy=T) 
names(hot60Area2_df)<-c("x","y","hot")
ggplot(hot60Area2_df) +
  geom_tile(aes(x = x, y = y, fill = factor(hot))) +
  scale_fill_manual(values = c("0" = "grey", "1" = "red")) +
  geom_spatvector(data=coast,fill=NA)+coord_sf(crs = crs,xlim=c(-160,165),ylim=c(-56,90))+
  labs(x = NULL, y = NULL,) +
  #theme_minimal()+
  theme_bw()+
  theme(
    legend.position="none"
  )

# popoulArea <- pop_new+poul_new
# popoulArea2 <- ifel(popoulArea==2,1,0)
# plot(popoulArea2)
#writeRaster(hot60Area2,'/root/autodl-tmp/humPoulResult/data/Hot_data/Pop50poul.tif')

AEpop <- entr_new+pop_new
AEpop2 <- ifel(AEpop==2,1,0)
plot(AEpop2)
AEpop2_df<-as.data.frame(AEpop2,xy=T) 
names(AEpop2_df)<-c("x","y","hot")
ggplot(AEpop2_df) +
  geom_tile(aes(x = x, y = y, fill = factor(hot))) +
  scale_fill_manual(values = c("0" = "grey", "1" = "red")) +
  geom_spatvector(data=coast,fill=NA)+coord_sf(crs = crs,xlim=c(-160,165),ylim=c(-56,90))+
  labs(x = NULL, y = NULL,) +
  #theme_minimal()+
  theme_bw()+
  theme(
    legend.position="none"
  )
#writeRaster(AEpop2,'/root/autodl-tmp/humPoulResult/data/Hot_data/Pop50entr.tif')
AEpoul <- entr_new+poul_new
AEpoul2 <- ifel(AEpoul==2,1,0)
plot(AEpoul2)
AEpoul2_df<-as.data.frame(AEpoul2,xy=T) 
names(AEpoul2_df)<-c("x","y","hot")
ggplot(AEpoul2_df) +
  geom_tile(aes(x = x, y = y, fill = factor(hot))) +
  scale_fill_manual(values = c("0" = "grey", "1" = "red")) +
  geom_spatvector(data=coast,fill=NA)+coord_sf(crs = crs,xlim=c(-160,165),ylim=c(-56,90))+
  labs(x = NULL, y = NULL,) +
  #theme_minimal()+
  theme_bw()+
  theme(
    legend.position="none"
  )
#writeRaster(AEpoul2,'/root/autodl-tmp/humPoulResult/data/Hot_data/Poul50entr.tif')


#求出暴露的人口总数和家禽总数

popd2015<-rast("/root/autodl-tmp/全球人口/GWP_v4/popd2015_30.tif") %>% resample(globalRaster)
Entropy<-rast(paste0(basePath,'AE_data/AE.tif'))%>% mask(globalCountry)
poul2015<-rast("/root/autodl-tmp/zyresult/Poultry_duckchic.tif")%>% resample(globalRaster)

globalSHP <- vect('/root/autodl-tmp/zyresult/Con_popentrpoul_sf_EU.shp')
globalSHP2 <- terra::aggregate(globalSHP,'name_r')
cRaster <- rasterize(globalSHP2,Entropy,field='name_r')

Numpop_hotArea<-hot60Area2*popd2015; plot(Numpop_hotArea)
Numpoul_hotArea<-hot60Area2*poul2015; plot(Numpoul_hotArea)
NumAEpop<-AEpop2*popd2015  ; plot(NumAEpop)
NumAEpoul<-AEpoul2*poul2015  ; plot(NumAEpoul)

All_df <-c(Entropy,Numpop_hotArea,Numpoul_hotArea,NumAEpop,NumAEpoul,cRaster) %>% terra::as.data.frame() %>%na.omit()
head(All_df)
names(All_df)<-c("entr","Numpop_hotArea","Numpoul_hotArea","NumAEpop","NumAEpoul","country")

Num_result <- All_df %>%
  group_by(country) %>%
  summarise(across(everything(), sum, na.rm = TRUE))
Num_result
#fwrite(Num_result,'/root/autodl-tmp/humPoulResult/data/Hot_data/hot_Num_result.csv')


#######(-)重新画一下entr约登指数但60poppoul的热点区域-----------
quan_pop <- quantiles_pop[1,61]
quan_entr<- 4.107
quan_poul <- quantiles_poul[1,61]
#重分类
pop_new <- ifel(popd2015>quan_pop,1,0)
entr_new <- ifel(Entropy>quan_entr,1,0)
poul_new <- ifel(poul2015>quan_poul,1,0)

hot60Area <- pop_new+entr_new+poul_new
hot60Area2 <- ifel(hot60Area==3,1,0)

plot(hot60Area2)
#writeRaster(hot60Area2,'/root/autodl-tmp/humPoulResult/data/Hot_data/Newpoppoul60entr.tif')
library(rnaturalearth)
coast <- ne_coastline(scale = "small", returnclass = "sf") %>% vect()
crs <- '+proj=longlat +datum=WGS84'
hot60Area2_df<-as.data.frame(hot60Area2,xy=T) 
names(hot60Area2_df)<-c("x","y","hot")
ggplot(hot60Area2_df) +
  geom_tile(aes(x = x, y = y, fill = factor(hot))) +
  scale_fill_manual(values = c("0" = "grey", "1" = "red")) +
  geom_spatvector(data=coast,fill=NA)+coord_sf(crs = crs,xlim=c(-160,165),ylim=c(-56,90))+
  labs(x = NULL, y = NULL,) +
  #theme_minimal()+
  theme_bw()+
  theme(
    legend.position="none"
  )

popoulArea <- pop_new+poul_new
popoulArea2 <- ifel(popoulArea==2,1,0)
plot(popoulArea2)

AEpop <- entr_new+pop_new
AEpop2 <- ifel(AEpop==2,1,0)
plot(AEpop2)
writeRaster(hot60Area2,'/root/autodl-tmp/humPoulResult/data/Hot_data/NewAEpop60.tif')
AEpoul <- entr_new+poul_new
AEpoul2 <- ifel(AEpoul==2,1,0)
plot(AEpoul2)
writeRaster(hot60Area2,'/root/autodl-tmp/humPoulResult/data/Hot_data/NewAEpoul60.tif')


#######(-)重新画一下entr约登指数但70poppoul的热点区域-----------
quan_pop <- quantiles_pop[1,71]
quan_entr<- 4.107
quan_poul <- quantiles_poul[1,71]
#重分类
pop_new <- ifel(popd2015>quan_pop,1,0)
entr_new <- ifel(Entropy>quan_entr,1,0)
poul_new <- ifel(poul2015>quan_poul,1,0)

hot60Area <- pop_new+entr_new+poul_new
hot60Area2 <- ifel(hot60Area==3,1,0)

plot(hot60Area2)
#writeRaster(hot60Area2,'/root/autodl-tmp/humPoulResult/data/Hot_data/Newpoppoul70entr.tif')
library(rnaturalearth)
coast <- ne_coastline(scale = "small", returnclass = "sf") %>% vect()
crs <- '+proj=longlat +datum=WGS84'
hot60Area2_df<-as.data.frame(hot60Area2,xy=T) 
names(hot60Area2_df)<-c("x","y","hot")
ggplot(hot60Area2_df) +
  geom_tile(aes(x = x, y = y, fill = factor(hot))) +
  scale_fill_manual(values = c("0" = "grey", "1" = "red")) +
  geom_spatvector(data=coast,fill=NA)+coord_sf(crs = crs,xlim=c(-160,165),ylim=c(-56,90))+
  labs(x = NULL, y = NULL,) +
  #theme_minimal()+
  theme_bw()+
  theme(
    legend.position="none"
  )

popoulArea <- pop_new+poul_new
popoulArea2 <- ifel(popoulArea==2,1,0)
plot(popoulArea2)

AEpop <- entr_new+pop_new
AEpop2 <- ifel(AEpop==2,1,0)
plot(AEpop2)
writeRaster(hot60Area2,'/root/autodl-tmp/humPoulResult/data/Hot_data/NewAEpop70.tif')
AEpoul <- entr_new+poul_new
AEpoul2 <- ifel(AEpoul==2,1,0)
plot(AEpoul2)
writeRaster(hot60Area2,'/root/autodl-tmp/humPoulResult/data/Hot_data/NewAEpoul70.tif')



#######(-)重新画一下entr约登指数但80poppoul的热点区域-----------
quan_pop <- quantiles_pop[1,81]
quan_entr<- 4.107
quan_poul <- quantiles_poul[1,81]
#重分类
pop_new <- ifel(popd2015>quan_pop,1,0)
entr_new <- ifel(Entropy>quan_entr,1,0)
poul_new <- ifel(poul2015>quan_poul,1,0)

hot60Area <- pop_new+entr_new+poul_new
hot60Area2 <- ifel(hot60Area==3,1,0)

plot(hot60Area2)
#writeRaster(hot60Area2,'/root/autodl-tmp/humPoulResult/data/Hot_data/Newpoppoul80entr.tif')
library(rnaturalearth)
coast <- ne_coastline(scale = "small", returnclass = "sf") %>% vect()
crs <- '+proj=longlat +datum=WGS84'
hot60Area2_df<-as.data.frame(hot60Area2,xy=T) 
names(hot60Area2_df)<-c("x","y","hot")
ggplot(hot60Area2_df) +
  geom_tile(aes(x = x, y = y, fill = factor(hot))) +
  scale_fill_manual(values = c("0" = "grey", "1" = "red")) +
  geom_spatvector(data=coast,fill=NA)+coord_sf(crs = crs,xlim=c(-160,165),ylim=c(-56,90))+
  labs(x = NULL, y = NULL,) +
  #theme_minimal()+
  theme_bw()+
  theme(
    legend.position="none"
  )

popoulArea <- pop_new+poul_new
popoulArea2 <- ifel(popoulArea==2,1,0)
plot(popoulArea2)

AEpop <- entr_new+pop_new
AEpop2 <- ifel(AEpop==2,1,0)
plot(AEpop2)
writeRaster(hot60Area2,'/root/autodl-tmp/humPoulResult/data/Hot_data/NewAEpop80.tif')
AEpoul <- entr_new+poul_new
AEpoul2 <- ifel(AEpoul==2,1,0)
plot(AEpoul2)
writeRaster(hot60Area2,'/root/autodl-tmp/humPoulResult/data/Hot_data/NewAEpoul80.tif')



#######(-)重新画一下entr约登指数但90poppoul的热点区域-----------
quan_pop <- quantiles_pop[1,91]
quan_entr<- 4.107
quan_poul <- quantiles_poul[1,91]
#重分类
pop_new <- ifel(popd2015>quan_pop,1,0)
entr_new <- ifel(Entropy>quan_entr,1,0)
poul_new <- ifel(poul2015>quan_poul,1,0)

hot60Area <- pop_new+entr_new+poul_new
hot60Area2 <- ifel(hot60Area==3,1,0)

plot(hot60Area2)
#writeRaster(hot60Area2,'/root/autodl-tmp/humPoulResult/data/Hot_data/Newpoppoul90entr.tif')
library(rnaturalearth)
coast <- ne_coastline(scale = "small", returnclass = "sf") %>% vect()
crs <- '+proj=longlat +datum=WGS84'
hot60Area2_df<-as.data.frame(hot60Area2,xy=T) 
names(hot60Area2_df)<-c("x","y","hot")
ggplot(hot60Area2_df) +
  geom_tile(aes(x = x, y = y, fill = factor(hot))) +
  scale_fill_manual(values = c("0" = "grey", "1" = "red")) +
  geom_spatvector(data=coast,fill=NA)+coord_sf(crs = crs,xlim=c(-160,165),ylim=c(-56,90))+
  labs(x = NULL, y = NULL,) +
  #theme_minimal()+
  theme_bw()+
  theme(
    legend.position="none"
  )

popoulArea <- pop_new+poul_new
popoulArea2 <- ifel(popoulArea==2,1,0)
plot(popoulArea2)

AEpop <- entr_new+pop_new
AEpop2 <- ifel(AEpop==2,1,0)
plot(AEpop2)
writeRaster(hot60Area2,'/root/autodl-tmp/humPoulResult/data/Hot_data/NewAEpop90.tif')
AEpoul <- entr_new+poul_new
AEpoul2 <- ifel(AEpoul==2,1,0)
plot(AEpoul2)
writeRaster(hot60Area2,'/root/autodl-tmp/humPoulResult/data/Hot_data/NewAEpoul90.tif')

#Fig4 Mantel test------------
#c###############
#分功能群计算活动熵
allDf2 <- fread('/root/autodl-tmp/zyresult/allDf786_reclass.csv')
all_rasters <- list.files('/root/autodl-tmp/result/singleModel', pattern = '.tif', full.names = TRUE) 
spName <- basename(all_rasters) %>% str_sub(.,1,-5)
for (fuc in unique(allDf2$`Functional Group`)) {
  fucName <- allDf2[allDf2$`Functional Group`==fuc,]
  paths <- all_rasters[spName%in%fucName$LatName]
  monthNum <- rast(paths) %>% sum(na.rm=T)
  calEntropy <- lapply(paths, function(x){
    r <- rast(x) %>% sum(.,na.rm=T)
    pi <- r/monthNum
    y <- -pi*log(pi)
    names(y) <- str_sub(basename(x),1,-5)
    return(y)
  })
  actEntropy <- rast(calEntropy) %>% sum(na.rm = T)
  writeRaster(actEntropy,paste0('/root/autodl-tmp/humPoulResult/data/AE_data/',fuc,'_','AE.tif'))
}

###合并环境数据和活动熵数据
Global_birds_Entropy<-rast('/root/autodl-tmp/zyresult/Entropy/Global_birds_Entropy.tif')
# pop<-rast("/root/autodl-tmp/zyresult/Mantel_Data/pop_pe-4.tif")
# poul<-rast("/root/autodl-tmp/zyresult/Mantel_Data/poul_pe-4.tif")
pop_log<-rast("/root/autodl-tmp/zyresult/Mantel_Data/popd2015_log.tif")
poul_log<-rast("/root/autodl-tmp/zyresult/Mantel_Data/poul2015_log.tif")
spNum<-rast('/root/autodl-tmp/zyresult/Mantel_Data/spNumTif.tif')
spNum_CV<-rast('/root/autodl-tmp/zyresult/Mantel_Data/sp_cv_than1.tif')
Allmonth<-rast('/root/autodl-tmp/zyresult/Entropy/Global_birds_allMonth.tif')
GPP<-rast('/root/autodl-tmp/zyresult/Mantel_Data/GPP0016_M_mean.tif')
GPP_CV<-rast('/root/autodl-tmp/zyresult/Mantel_Data/GPP0016_M_CV.tif')
below0Days<-rast('/root/autodl-tmp/zyresult/Mantel_Data/below0Days.tif')
bio1temp<-rast('/root/autodl-tmp/zyresult/Mantel_Data/BIO_01.tif');names(bio1temp)
bio12prec<-rast('/root/autodl-tmp/zyresult/Mantel_Data/BIO_12.tif');names(bio12prec)
otherClimateData<-rast('/root/autodl-tmp/otherClimateData/climateBirdData.tif');names(otherClimateData)

allPath <- list.files()
allData <- rast(allPath) %>% resample(globalRaster) %>% crop(globalRaster,mask=T)

countryRaster <- vect("/root/autodl-tmp/zyresult/Con_dissoEU.shp") %>% rasterize(.,allData[[1]],field='name_r')
allData_df <- c(allData,countryRaster) %>% as.data.frame()

plot
mantel <- mantel_test(China_birdEntropy_df_clean,   ## 分类数据 China_birdEntropy_df_clean
                      China_corData_df_clean,    ## 影响因子数据 China_corData_df_clean
                      spec_select = list(
                        seabirds = 1,
                        shorebirds = 2,
                        wadingbirds = 3,
                        waterfowls = 4,
                        unclassified = 5
                      )) %>%
  mutate(rd = cut(r, breaks = c(-Inf, 0 ,0.3, 0.5, 0.7, Inf),
                  labels = c("< 0","< 0.3", "0.3 - 0.5", "0.5 - 0.7", ">= 0.7")),
         pd = cut(p, breaks = c(-Inf, 0.01, 0.05, Inf),
                  labels = c("< 0.01", "0.01 - 0.05", ">= 0.05")))

fwrite(mantel,'/root/autodl-tmp/zyresult/Mantel_Data/China_Mantel.csv')


mantel <- fread('/root/autodl-tmp/zyresult/Mantel_Data/China_Mantel.csv')
# mantel$type <- case_when(mantel$)
#处理数据
head(mantel)
mantel2<-subset(mantel, r>0.4)
mantel2 <- mantel2 %>%
  mutate(rd2 = cut(r, breaks = c(-Inf, 0, 0.4, 0.5, 0.6, 0.7, 0.8, Inf),
                   labels = c("< 0", "< 0.4", "0.4 - 0.5","0.5 - 0.6", "0.6 - 0.7", "0.7 - 0.8", ">= 0.8")))
#改名字
names(China_corData_df_clean)<-
  c('AE','Species Richness', "Cumulative species-months", 'CV of Species Richness',
    'Precipitation', 'Temperature', 'Frost Days',
    'GPP', 'CV of GPP', 'NPP', 
    'NDVI', 'NDWI', 'LAI', 
    'Lontitude', 'Latitude' , 'Light Density', 'Road Density', 
    'Human density(Log)', 'Poultry density(Log)'
  ) 

China_corData_df_clean2<-China_corData_df_clean[,c('AE','Species Richness', "Cumulative species-months", 'CV of Species Richness',
                                                   'Lontitude', 'Latitude' , 'Precipitation', 'Temperature', 'Frost Days',
                                                   'GPP', 'CV of GPP', 'NPP', 
                                                   'NDVI', 'NDWI', 'LAI', 
                                                   'Light Density', 'Road Density', 
                                                   'Human density(Log)', 'Poultry density(Log)')]

mantel2 <- mantel2 %>%
  mutate(env = case_when(
    env == "Waterbird Activity Entropy" ~ "AE",
    env == "Months of Species Richness"~"Cumulative species-months",
    env == "Subzero Temperature Days" ~ "Frost Days",
    env == "Population_Log"~"Human density(Log)",
    env == "Poultry_Log"~"Poultry density(Log)",
    # 添加更多的替换规则
    TRUE ~ env  # 默认情况下保持原值
  ))%>%
  mutate(spec = case_when(
    spec == "seabirds" ~ "Seabirds",
    spec == "shorebirds" ~ "Shorebirds",
    spec == "wadingbirds" ~ "Large wading birds",
    spec == "waterfowls" ~ "Waterfowls",
    spec == "unclassified" ~ "Others",
    # 添加更多的替换规则
    TRUE ~ spec  # 默认情况下保持原值
  ))

qcorrplot(China_corData_df_clean2, type = "lower", diag = FALSE) +  #China_corData_df_clean
  geom_square() +   ## 相关性热图的形状
  ## 
  geom_couple(aes(colour = rd2, linetype = pd), size=1,
              data = mantel2, label.size = 6,
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
  Country = c("China", "India", "European Union", "Nigeria", "United States"),
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
    Attribute == "CV_of_GPP" ~ "CV of GPP",
    TRUE ~ Attribute  # 默认情况下保持原值
  ))


vegetationdata_long2$Attribute<-factor(vegetationdata_long2$Attribute,levels = c("GPP","CV of GPP","NPP","NDVI","NDWI","LAI"))





#####人类活动
humandata <- data.frame(
  Country = c("China", "United States", "India", "European Union", "Nigeria"),
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
  scale_fill_manual(values = c("India"="#377EB8", "China"="#E41A1C", "European Union"="#4DAF4A", "Nigeria"="#984EA3", "United States"="#f7a156")) +  #FF7F00
  #scale_fill_manual(values = c("#377EB8", "#E41A1C", "#4DAF4A", "#984EA3", "#FF7F00")) +  #FF7F00
  #scale_fill_manual(values = c("India"="#4C8BC0", "China"="#C75C64", "European Union"="#32A852", "Nigeria"="#984EA3", "United States"="#F0B57D")) +  #F0B57D
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
  scale_fill_manual(values = c("India"="#377EB8", "China"="#E41A1C", "European Union"="#4DAF4A", "Nigeria"="#984EA3", "United States"="#f7a156")) +  #FF7F00
  #scale_fill_manual(values = c("#377EB8", "#E41A1C", "#4DAF4A", "#984EA3", "#FF7F00")) +  #FF7F00
  #scale_fill_manual(values = c("India"="#4C8BC0", "China"="#C75C64", "European Union"="#32A852", "Nigeria"="#984EA3", "United States"="#F0B57D")) +  #F0B57D
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
globalSHP2 <- terra::aggregate(globalSHP,'name_r')
Con_5 <- subset(globalSHP2, globalSHP2$name_r %in% c("China", "India", "European Union", "Nigeria", "United States"))
cRaster <- rasterize(Con_5,Entropy,field='name_r')
All5_df <- c(Entropy,popd2015,poul2015,cRaster) %>% terra::as.data.frame() %>%na.omit()
names(All5_df)<-c("entr","pop","poul","country")




#####人口线性
# 过滤掉pop值小于0.0001的数据行
All5_df1 <- All5_df[!(All5_df$pop < 0.0001),]

p1<-ggplot(data = All5_df1, aes(x = entr, y = log(pop), color = country)) +
  geom_point(size=1) +  #,shape=21
  scale_color_manual(values = c("India"="#377EB8", "China"="#E41A1C", "European Union"="#4DAF4A", "Nigeria"="#984EA3", "United States"="#f7a156")) +  #FF7F00
  #scale_color_manual(values = c("#377EB8", "#E41A1C", "#4DAF4A", "#984EA3", "#FF7F00")) +
  geom_smooth(method = "lm",linewidth=0.5, se = T, lty=5,aes(color = country)) +
  #geom_smooth(method = "lm", se = T, color = "black", linetype = "dashed") +  # 全部点的趋势线，黑色
  #geom_smooth(method = "lm",linewidth=0.5, se = T, lty=5,color = "black") +  #分界
  # stat_poly_eq(aes(label = after_stat(paste(eq.label, rr.label,p.value.label, sep = "~~~"))),  #,p.value.label  eq.label,
  #              formula = y ~ x,
  #              coef.digits = 3, rr.digits = 2, parse = TRUE, size =3) +
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
EU_count <- nrow(All5_df1[All5_df1$country == 'European Union', ]); EU_count  # 2273
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

model <- lm(log(pop) ~ entr, data = All5_df1[All5_df1$country == 'European Union', ])
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
  scale_color_manual(values = c("India"="#377EB8", "China"="#E41A1C", "European Union"="#4DAF4A", "Nigeria"="#984EA3", "United States"="#f7a156")) +  #FF7F00
  #scale_color_manual(values = c("#377EB8", "#E41A1C", "#4DAF4A", "#984EA3", "#FF7F00")) +
  geom_smooth(method = "lm",linewidth=0.5, se = T, lty=5,aes(color = country)) +
  #geom_smooth(method = "lm", se = T, color = "black", linetype = "dashed") +  # 全部点的趋势线，黑色
  #geom_smooth(method = "lm",linewidth=0.5, se = T, lty=5,color = "black") +  #分界
  # stat_poly_eq(aes(label = after_stat(paste(eq.label, rr.label,p.value.label, sep = "~~~"))),
  #              formula = y ~ x,
  #              coef.digits = 3, rr.digits = 2, parse = TRUE, size = 3) +
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
EU_count <- nrow(All5_df2[All5_df2$country == 'European Union', ]); EU_count  # 2163
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

model <- lm(log(poul) ~ entr, data = All5_df2[All5_df2$country == 'European Union', ])
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





#Fig.5 ############

#a###########


#####饼图，最重要的------------

lineardatalong <-fread('/root/autodl-tmp/zyresult/Entropy/poppulresult/0-Conbine_poppoul_long.csv')
lineardatalong<-lineardatalong %>%
  mutate(Pd = cut(P, breaks = c(-Inf, 0.01, 0.05, Inf),
                  labels = c("< 0.01", "0.01 - 0.05", "> 0.05")))%>%
  mutate(R2d = cut(R2, breaks = c(-Inf, 0.3, 0.4,0.5,0.6, Inf),
                   labels = c("< 0.3", "0.3 - 0.4", "0.4 - 0.5", "0.5 - 0.6", "> 0.6")))
head(lineardatalong)

# 定义一个函数来筛选每个组合的最大R2值对应的行
get_max_r2_rows <- function(dt, type) {
  max_r2 <- dt[Type == type, .(Max_R2 = max(R2)), by = .(Country, Part)]    # 计算每个Country和Part组合的最大R2值
  setkey(dt, Country, Part, R2)                # 将最大R2值合并回原始数据集
  setkey(max_r2, Country, Part, Max_R2)
  max_rows <- max_r2[dt, nomatch=0]      # 筛选出最大R2值对应的行
  #max_rows[, Max_R2 := NULL]    # 从结果中移除Max_R2列
  return(max_rows)
}

# 分别为ALL和CH类型筛选数据
max_r2_rows_ALL <- get_max_r2_rows(lineardatalong, "ALL") 
max_r2_rows_CH <- get_max_r2_rows(lineardatalong, "CH")
lineardata_max <- rbind(max_r2_rows_ALL, max_r2_rows_CH)
print(lineardata_max)

# 使用ifelse语句来创建颜色列
lineardata_max$Color <- ifelse(lineardata_max$Type == "ALL" & lineardata_max$R2d == "< 0.3", "#F5EAD0",
                               ifelse(lineardata_max$Type == "ALL" & lineardata_max$R2d == "0.3 - 0.4", "#E3CB8F",
                                      ifelse(lineardata_max$Type == "ALL" & lineardata_max$R2d == "0.4 - 0.5", "#BF8D44",
                                             ifelse(lineardata_max$Type == "ALL" & lineardata_max$R2d == "0.5 - 0.6", "#8E581C",
                                                    ifelse(lineardata_max$Type == "ALL" & lineardata_max$R2d == "> 0.6", "#59310D",
                                                           ifelse(lineardata_max$Type == "CH" & lineardata_max$R2d == "< 0.3", "#FAE3D7",
                                                                  ifelse(lineardata_max$Type == "CH" & lineardata_max$R2d == "0.3 - 0.4", "#F2B396",
                                                                         ifelse(lineardata_max$Type == "CH" & lineardata_max$R2d == "0.4 - 0.5", "#D9705A",
                                                                                ifelse(lineardata_max$Type == "CH" & lineardata_max$R2d == "0.5 - 0.6", "#B42D34",
                                                                                       ifelse(lineardata_max$Type == "CH" & lineardata_max$R2d == "> 0.6", "#6E0D20",
                                                                                              NA))))))))))
lineardata_max

######——制作图例---------
colors_ALL <- c("< 0.3"="#F5EAD0", "0.3 - 0.4" = "#E3CB8F", "0.4 - 0.5" = "#BF8D44", "0.5 - 0.6" = "#8E581C", "> 0.6" = "#59310D")
#colors_CH <- c("< 0.3"="#E7F2F1", "0.3 - 0.4" = "#AFDCD5", "0.4 - 0.5" = "#6DB2AD", "0.5 - 0.6" = "#347C78", "> 0.6" = "#174A41")
colors_CH <- c("< 0.3"="#FAE3D7", "0.3 - 0.4" = "#F2B396", "0.4 - 0.5" = "#D9705A", "0.5 - 0.6" = "#B42D34", "> 0.6" = "#6E0D20")

data_ALL <- lineardata_max[lineardata_max$Type == "ALL", ]
data_CH <- lineardata_max[lineardata_max$Type == "CH", ]
p1<-ggplot() +
  geom_bar(data=data_ALL, aes(x=Country, y=Max_R2, fill=R2d), stat="identity", position=position_dodge())+
  scale_fill_manual(name="R² (ALL)", values=colors_ALL, breaks=names(colors_ALL), labels=names(colors_ALL))+
  theme(text = element_text(size=20),
        legend.position = "bottom",
        legend.background = element_rect(fill = "transparent")) # 透明的图例背景
#legend.box.background = element_rect(fill = "transparent")) # 透明的图例面板)
p2<-ggplot() +
  geom_bar(data=data_CH, aes(x=Country, y=Max_R2, fill=R2d), stat="identity", position=position_dodge())+
  scale_fill_manual(name="R² (CH)", values=colors_CH, breaks=names(colors_CH), labels=names(colors_CH))+
  theme(text = element_text(size=20),
        legend.position = "bottom",
        legend.background = element_rect(fill = "transparent")) # 透明的图例背景
#legend.box.background = element_rect(fill = "transparent")) # 透明的图例面板)

p1+p2
ggsave("/root/autodl-tmp/zyresult/Entropy/poppulresult/piedata_Lengend.png", p1+p2, bg = "transparent", width = 18, height = 8, units = "in")

######——制作底图---------
#世界地图
world.map <- rnaturalearth::ne_countries(returnclass = "sf") %>% filter(continent != "Antarctica")

Con_popentrpoul_sf <- st_read("/root/autodl-tmp/zyresult/Con_popentrpoul_sf_EU.shp");head(Con_popentrpoul_sf)
Con_ALLina<-subset(Con_popentrpoul_sf,name_r=="China")
Con_India<-subset(Con_popentrpoul_sf,name_r=="India") 
Con_EU<-subset(Con_popentrpoul_sf,name_r=="European Union")
Con_Nigeria<-subset(Con_popentrpoul_sf,name_r=="Nigeria")
Con_US<-subset(Con_popentrpoul_sf,name_r=="United States") 

ggplot() +
  geom_sf(data = world.map, fill = "lightgrey", color = NA) + 
  geom_sf(data = Con_China, fill = "#baa1a1", color = NA) +    #lightcoral
  geom_sf(data = Con_India, fill = "#baa1a1", color = NA) + 
  geom_sf(data = Con_EU, fill = "#baa1a1", color = NA) + 
  geom_sf(data = Con_Nigeria, fill = "#baa1a1", color = NA) + 
  geom_sf(data = Con_US, fill = "#baa1a1", color = NA) + 
  coord_sf(crs = "+proj=longlat +datum=WGS84", xlim = c(-160, 165), ylim = c(-56, 90)) +
  theme_bw() + # 使用简洁主题
  theme(legend.position = "none")



#####【China饼图】------
piedata_China <- lineardata_max[lineardata_max$Country == "China", ]
piedata_China$Value<-c(1,1,1,1)

# 顺序
desired_order <- c("pop_ALL", "pop_CH", "poul_CH", "poul_ALL")
piedata_China$Order <- factor(paste(piedata_China$Part, piedata_China$Type, sep = "_"), levels = desired_order)
piedata_China <- piedata_China[order(piedata_China$Order), ]
piedata_China$Color <- factor(piedata_China$Color, levels = unique(piedata_China$Color))
#c("#8E581C","#D9705A","#B42D34","#59310D")

p <- ggplot(piedata_China, aes(x = "", y = Value, fill = Color)) +
  geom_bar(stat = "identity", width = 1) +
  coord_polar("y", start = 0) +
  scale_fill_manual(values = levels(piedata_China$Color)) +
  theme_void() +
  theme(
    panel.background = element_rect(fill = "transparent"), # 透明的面板背景
    plot.background = element_rect(fill = "transparent", color = NA), # 透明的绘图背景
    legend.position = "none"
  )
p
# 保存为 PNG 图像，背景透明
ggsave("/root/autodl-tmp/zyresult/Entropy/poppulresult/piedata_China.png", p, bg = "transparent", width = 10, height = 8, units = "in")


#####【India 饼图】------
piedata_India <- lineardata_max[lineardata_max$Country == "India", ]
piedata_India$Value<-c(1,1,1,1)

# 顺序
desired_order <- c("pop_ALL", "pop_CH", "poul_CH", "poul_ALL")
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


#####【EU 饼图】------
piedata_EU <- lineardata_max[lineardata_max$Country == "EU", ]
piedata_EU$Value<-c(1,1,1,1)
color_EU<-c("pop_ALL"="#BF8D44","pop_CH"="#FAE3D7","poul_CH"="#FAE3D7","poul_ALL"="#BF8D44")
# 顺序
desired_order <- c("pop_ALL", "pop_CH", "poul_CH", "poul_ALL")
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

#####【Nigeria 饼图】------
piedata_Nigeria <- lineardata_max[lineardata_max$Country == "Nigeria", ]
piedata_Nigeria$Value<-c(1,1,1,1)

color_Nigeria<-c("pop_ALL"="#F5EAD0","pop_CH"="#FAE3D7","poul_CH"="#FAE3D7","poul_ALL"="#F5EAD0")

# 顺序
desired_order <- c("pop_ALL", "pop_CH", "poul_CH", "poul_ALL")
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


#####【US 饼图】------
piedata_US <- lineardata_max[lineardata_max$Country == "US", ]
piedata_US$Value<-c(1,1,1,1)

color_US<-c("pop_ALL"="#E3CB8F","pop_CH"="#F2B396","poul_CH"="#FAE3D7","poul_ALL"="#E3CB8F")

# 顺序
desired_order <- c("pop_ALL", "pop_CH", "poul_CH", "poul_ALL")
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










#b############
#####弦图，ALL和CH----------

lineardatalong <-fread('/root/autodl-tmp/zyresult/Entropy/poppulresult/0-Conbine_poppoul_long.csv')
lineardatalong<-lineardatalong %>%
  mutate(Pd = cut(P, breaks = c(-Inf, 0.01, 0.05, Inf),
                  labels = c("< 0.01", "0.01 - 0.05", "> 0.05")))%>%
  mutate(R2d = cut(R2, breaks = c(-Inf, 0.3, 0.4,0.5,0.6, Inf),
                   labels = c("< 0.3", "0.3 - 0.4", "0.4 - 0.5", "0.5 - 0.6", "> 0.6")))%>%
  subset(., Type!="SH")
head(lineardatalong)

#逐个画：
lineardata_ALLpop<-subset(lineardatalong, Type=="ALL" & Part == "pop")
lineardata_ALLpoul<-subset(lineardatalong, Type=="ALL" & Part == "poul")
lineardata_CHpop<-subset(lineardatalong, Type=="CH" & Part == "pop")
lineardata_CHpoul<-subset(lineardatalong, Type=="CH" & Part == "poul")

library(circlize)

######__ ALLpop创建矩阵格式--------
# mat_ALLpop <- dcast(lineardata_ALLpop, Variable ~ Country, value.var = "R2")
# mat_ALLpop2 <- as.matrix(mat_ALLpop[, -1])
# rownames(mat_ALLpop2) <- mat_ALLpop$Variable
# mat_ALLpop2 <- mat_ALLpop2[c("Seabirds", "Shorebirds", "Waterfowls", "Wading bird", "unclassified"),]
# rownames(mat_ALLpop2)<-c("Seabirds", "Shorebirds", "Waterfowls", "Large Wading birds","unclassified")
# colnames(mat_ALLpop2)<-c("China", "India", "EU", "Nigeria","US")
# chordDiagram(mat_ALLpop2)

#国家在下面
mat_ALLpop <- dcast(lineardata_ALLpop, Country ~ Variable, value.var = "R2")
mat_ALLpop2 <- as.matrix(mat_ALLpop[, -1])
rownames(mat_ALLpop2) <- mat_ALLpop$Country
#mat_ALLpop2 <- mat_ALLpop2[c("Seabirds", "Shorebirds", "Waterfowls", "Wading bird", "unclassified"),]
colnames(mat_ALLpop2)<-c("Seabirds", "Shorebirds", "Waterfowls", "Large Wading birds","Others")
rownames(mat_ALLpop2)<-c("China", "India", "EU", "Nigeria","US")
chordDiagram(mat_ALLpop2)

#设定颜色
grid.col = NULL
grid.col[c("Seabirds", "Shorebirds", "Waterfowls", "Large Wading birds","Others")] = c("#E3CB8F", "#E3CB8F","#E3CB8F",  "#E3CB8F", "#E3CB8F") #c("#E41A1C","#4DAF4A","#377EB8","#984EA3","#f7a156") 
grid.col[rownames(mat_ALLpop2)] = c("#E41A1C","#377EB8","#4DAF4A","#984EA3","#f7a156")  #c("#FF6699","#33CCFF","#99FFC9","#E600E6","#FFD700")
#scale_fill_manual(values = c("#FF7F00", "#984EA3", "#4DAF4A", "#E41A1C", "#377EB8")) +
#此处为分段设置，因此为行或列数目-1
# circos.par(gap.degree = c(rep(2, nrow(mat_ALLpop2)-1), 5, rep(2, ncol(mat_ALLpop2)-1), 5), 
#            start.degree = 0)  #从什么角度开始绘图，本设定为从0度开始，则上方为分组，下方为组分
# 设置全局字体大小
# 设置全局字体大小
par(cex = 2.4) # 这里的1.5是字体大小的倍数
chordDiagram(mat_ALLpop2,
             #annotationTrack = "grid", # 不创建表示标签的轨迹
             big.gap = 10,
             diffHeight = 0.06, # 外圈与连线间隔高度
             grid.col = grid.col, # 线条颜色
             link.lwd = 0.02, # 线条宽度
             transparency = 0.5) # 连接颜色透明度

circos.clear()



######__ ALLpoul创建矩阵格式--------
#国家在下面
mat_ALLpoul <- dcast(lineardata_ALLpoul, Country ~ Variable, value.var = "R2")
mat_ALLpoul2 <- as.matrix(mat_ALLpoul[, -1])
rownames(mat_ALLpoul2) <- mat_ALLpoul$Country
#mat_ALLpoul2 <- mat_ALLpoul2[c("Seabirds", "Shorebirds", "Waterfowls", "Wading bird", "unclassified"),]
colnames(mat_ALLpoul2)<-c("Seabirds", "Shorebirds", "Waterfowls", "Large Wading birds","Others")
rownames(mat_ALLpoul2)<-c("China", "India", "EU", "Nigeria","US")
chordDiagram(mat_ALLpoul2)

#设定颜色
grid.col = NULL
grid.col[c("Seabirds", "Shorebirds", "Waterfowls", "Large Wading birds","Others")] = c("#E3CB8F", "#E3CB8F","#E3CB8F",  "#E3CB8F", "#E3CB8F") #c("#E41A1C","#4DAF4A","#377EB8","#984EA3","#f7a156") 
grid.col[rownames(mat_ALLpoul2)] = c("#E41A1C","#377EB8","#4DAF4A","#984EA3","#f7a156")  #c("#FF6699","#33CCFF","#99FFC9","#E600E6","#FFD700")
#scale_fill_manual(values = c("#FF7F00", "#984EA3", "#4DAF4A", "#E41A1C", "#377EB8")) +
#此处为分段设置，因此为行或列数目-1
# circos.par(gap.degree = c(rep(2, nrow(mat_ALLpoul2)-1), 5, rep(2, ncol(mat_ALLpoul2)-1), 5), 
#            start.degree = 0)  #从什么角度开始绘图，本设定为从0度开始，则上方为分组，下方为组分
# 设置全局字体大小
# 设置全局字体大小
par(cex = 2.3) # 这里的1.5是字体大小的倍数
chordDiagram(mat_ALLpoul2,
             #annotationTrack = "grid", # 不创建表示标签的轨迹
             big.gap = 10,
             diffHeight = 0.06, # 外圈与连线间隔高度
             grid.col = grid.col, # 线条颜色
             link.lwd = 0.02, # 线条宽度
             transparency = 0.5) # 连接颜色透明度

circos.clear()



######__ CHpop创建矩阵格式--------
#国家在下面
mat_CHpop <- dcast(lineardata_CHpop, Country ~ Variable, value.var = "R2")
mat_CHpop2 <- as.matrix(mat_CHpop[, -1])
rownames(mat_CHpop2) <- mat_CHpop$Country
#mat_CHpop2 <- mat_CHpop2[c("Seabirds", "Shorebirds", "Waterfowls", "Wading bird", "unclassified"),]
colnames(mat_CHpop2)<-c("Seabirds", "Shorebirds", "Waterfowls", "Large Wading birds","Others")
rownames(mat_CHpop2)<-c("China", "India", "EU", "Nigeria","US")
chordDiagram(mat_CHpop2)

#设定颜色
grid.col = NULL
grid.col[c("Seabirds", "Shorebirds", "Waterfowls", "Large Wading birds","Others")] = c("#D9705A", "#D9705A","#D9705A",  "#D9705A", "#D9705A") #c("#E41A1C","#4DAF4A","#377EB8","#984EA3","#f7a156") 
grid.col[rownames(mat_CHpop2)] = c("#E41A1C","#377EB8","#4DAF4A","#984EA3","#f7a156")  #c("#FF6699","#33CCFF","#99FFC9","#E600E6","#FFD700")
#scale_fill_manual(values = c("#FF7F00", "#984EA3", "#4DAF4A", "#E41A1C", "#377EB8")) +
#此处为分段设置，因此为行或列数目-1
# circos.par(gap.degree = c(rep(2, nrow(mat_CHpop2)-1), 5, rep(2, ncol(mat_CHpop2)-1), 5), 
#            start.degree = 0)  #从什么角度开始绘图，本设定为从0度开始，则上方为分组，下方为组分
# 设置全局字体大小
# 设置全局字体大小
par(cex = 2.3) # 这里的1.5是字体大小的倍数
chordDiagram(mat_CHpop2,
             #annotationTrack = "grid", # 不创建表示标签的轨迹
             big.gap = 10,
             diffHeight = 0.06, # 外圈与连线间隔高度
             grid.col = grid.col, # 线条颜色
             link.lwd = 0.02, # 线条宽度
             transparency = 0.5) # 连接颜色透明度

circos.clear()


######__ CHpoul创建矩阵格式--------
#国家在下面
mat_CHpoul <- dcast(lineardata_CHpoul, Country ~ Variable, value.var = "R2")
mat_CHpoul2 <- as.matrix(mat_CHpoul[, -1])
rownames(mat_CHpoul2) <- mat_CHpoul$Country
#mat_CHpoul2 <- mat_CHpoul2[c("Seabirds", "Shorebirds", "Waterfowls", "Wading bird", "unclassified"),]
colnames(mat_CHpoul2)<-c("Seabirds", "Shorebirds", "Waterfowls", "Large Wading birds","Others")
rownames(mat_CHpoul2)<-c("China", "India", "EU", "Nigeria","US")
chordDiagram(mat_CHpoul2)

#设定颜色
grid.col = NULL
grid.col[c("Seabirds", "Shorebirds", "Waterfowls", "Large Wading birds","Others")] = c("#D9705A", "#D9705A","#D9705A",  "#D9705A", "#D9705A") #c("#E41A1C","#4DAF4A","#377EB8","#984EA3","#f7a156") 
grid.col[rownames(mat_CHpoul2)] = c("#E41A1C","#377EB8","#4DAF4A","#984EA3","#f7a156")  #c("#FF6699","#33CCFF","#99FFC9","#E600E6","#FFD700")
#scale_fill_manual(values = c("#FF7F00", "#984EA3", "#4DAF4A", "#E41A1C", "#377EB8")) +
#此处为分段设置，因此为行或列数目-1
# circos.par(gap.degree = c(rep(2, nrow(mat_CHpoul2)-1), 5, rep(2, ncol(mat_CHpoul2)-1), 5), 
#            start.degree = 0)  #从什么角度开始绘图，本设定为从0度开始，则上方为分组，下方为组分
# 设置全局字体大小
# 设置全局字体大小
par(cex = 2.2) # 这里的1.5是字体大小的倍数
chordDiagram(mat_CHpoul2,
             #annotationTrack = "grid", # 不创建表示标签的轨迹
             big.gap = 10,
             diffHeight = 0.06, # 外圈与连线间隔高度
             grid.col = grid.col, # 线条颜色
             link.lwd = 0.02, # 线条宽度
             transparency = 0.5) # 连接颜色透明度

circos.clear()



