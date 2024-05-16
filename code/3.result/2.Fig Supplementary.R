
basePath<-"/root/autodl-tmp/humPoulResult/data/"

world.map <- rnaturalearth::ne_countries(returnclass = "sf") |>filter(continent != "Antarctica")
globalCountry <- vect(world.map) 
globalRaster <- rast(vals=1:259200,nrows=360, ncols=720,xmin=-180, xmax=180,ymin=-90, ymax=90,crs=crs)
coast <- ne_coastline(scale = "small", returnclass = "sf")
crs <- '+proj=longlat +datum=WGS84'
allDf <- fread(paste0(basePath,'allDf786_reclass.csv'))
speciesPixelNumPath <- list.files('/root/autodl-tmp/humPoulResult/data/single_model',pattern = '.tif',full.names = T)
spName <- basename(speciesPixelNumPath) %>% str_sub(.,1,-5)                      
speciesPixelNumPath2 <- speciesPixelNumPath[spName%in%allDf$LatName]



#Fig S2 CV of SR&GPP#########
sp_cv<-rast(paste0(basePath,'AE_data/sp_cv.tif'))%>% resample(globalRaster) %>% crop(globalRaster,mask=T)
gpp_cv<-rast(paste0(basePath,'AE_data/GPP0016_M_CV.tif'))  %>% resample(globalRaster) %>% crop(.,globalRaster,mask=T)
crs(gpp_cv) <- crs

Vborder <- vect('/root/autodl-tmp/Wallace_zoogeographic/newValisBorder.shp') 
vRaster <- rasterize(Vborder,sp_cv,field='name')
valisCV <- c(sp_cv,gpp_cv,vRaster) %>% terra::as.data.frame(xy=T) %>%na.omit()
names(valisCV)<-c("lon","lat","sp_cv","gpp_cv","name")


#归一化
valisCV$sp <- scale(valisCV$sp_cv, center = min(valisCV$sp_cv), scale = max(valisCV$sp_cv) - min(valisCV$sp_cv))
valisCV$gpp <- scale(valisCV$gpp_cv, center = min(valisCV$gpp_cv), scale = max(valisCV$gpp_cv) - min(valisCV$gpp_cv))
valisCV <- valisCV[valisCV$name != "Antarctica",]   #去掉不需要的南极界
summary(valisCV)


library(ggplot2)
library(ggpmisc)
library(ggpointdensity)


###分南北颜色分各界形状mian图-------
library(ggExtra)
#z2$Hemisphere <- ifelse(z2$Vborder %in% c("Neotropical", "Ethiopian", "Australian"), "Southern", "Northern")
valisCV$Hemisphere <- ifelse(valisCV$lat<0, "Southern", "Northern")
head(valisCV)
# 然后，使用这个新变量来设置散点的颜色
p<-ggplot(data = valisCV, aes(x = gpp, y = sp, color = Hemisphere)) +
  geom_smooth(data = subset(valisCV, Hemisphere == "Northern"), aes(group = 1), method = "lm", color = "#E08784", linetype = "dashed") + # 北半球线性拟合
  geom_point(aes(shape = name))+
  geom_smooth(data = subset(valisCV, Hemisphere == "Southern"), aes(group = 1), method = "lm", color = "#33969E", linetype = "dashed") + # 南半球线性拟合
  #geom_pointdensity(adjust = 2, size = 0.5) + 
  scale_color_manual(values = c("Southern" = "#00BFC4", "Northern" = "#E08784")) + # 南半球为蓝色，北半球为红色
  scale_shape_manual(values = c(16, 1, 0, 17, 15, 4)) +
  scale_size_area(max_size = 10) +
  #geom_line(data = new_data, aes(x = gpp, y = sp), color = "grey",lty=2) +
  stat_poly_eq(aes(label = paste(..eq.label.., ..rr.label.., sep = "~~")), formula = y ~ x, parse = TRUE) + # 添加公式和R^2
  labs(x = "CV of GPP", y = "CV of Species Richness") +
  theme(
    # axis.title.x = element_text(size=20),
    # axis.title.y = element_text(size=20),
    # axis.text = element_text(size=18),
    text=element_text(size=14),
    panel.background = element_rect(fill = "white"),
    panel.border = element_rect(color = "black", linewidth = 1, fill = NA),
    panel.grid = element_blank(),
    legend.position = "bottom",
    legend.title = element_blank(),
    plot.title = element_text(size = 14, hjust = 0.5),
    axis.title = element_text(size = 14)
  )
p
#边缘直方图
library(ggExtra)
ggMarginal(p, type = "histogram", groupColour = TRUE, groupFill = TRUE)





#####大界分面展示supply图-------
ggplot(data = valisCV, aes(x = gpp, y = sp, color=Hemisphere)) +
  facet_wrap(~name,ncol=2,scales ="free_y")+ #依据Attribute进行分面绘制
  geom_point(size=0.3)+
  geom_smooth(data = subset(valisCV, Hemisphere == "Northern"), aes(group = 1), method = "lm", color = "#E08784", linetype = "dashed",linewidth=0.5) + # 北半球线性拟合
  geom_smooth(data = subset(valisCV, Hemisphere == "Southern"), aes(group = 1), method = "lm", color = "#33969E", linetype = "dashed",linewidth=0.5) + # 南半球线性拟合
  #geom_pointdensity(adjust = 2, size = 0.5) + 
  coord_cartesian(ylim = c(0, 1))+ #限定X轴的显示范围
  scale_color_manual(values = c("Southern" = "#00BFC4", "Northern" = "#E08784")) + # 南半球为蓝色，北半球为红色
  scale_size_area(max_size = 10) +
  #geom_line(data = new_data, aes(x = gpp, y = sp), color = "grey",lty=2) +
  stat_poly_eq(aes(label = paste(..eq.label.., ..rr.label.., sep = "~~")), formula = y ~ x, parse = TRUE,size=2) + # 添加公式和R^2
  # stat_poly_eq(aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")), 
  #              label.x.npc = "right", label.y.npc = 60,
  #              coef.digits = 3, rr.digits = 2, parse = TRUE, size = 3, color = "black")+ # 绘制线性回归方程和决定系数
  labs(x = "CV of GPP", y = "CV of Species Richness") +
  theme(
    panel.background = element_rect(fill = "white"),
    panel.border = element_rect(color = "black", linewidth = 1, fill = NA),
    panel.grid = element_blank(),
    legend.position = "none",
    legend.title = element_blank(),
    plot.title = element_text(size = 14, hjust = 0.5),
    axis.title = element_text(size = 12)
  )


