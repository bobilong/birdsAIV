


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




#Fig S1 流程图#########
allDf <- fread(paste0(basePath,'Supplementary Tables.csv'))
head(allDf)
list_of_tables <- lapply(allDf, table)
list_of_tables$Order
list_of_tables$`Functional Group`

functional_order_table <- table(allDf$`Functional Group`, allDf$Order)
functional_order_table




#Fig S2 鸿雁的逐月物种分布#########

ggplot() +
  geom_spatraster(data = waterPer,show.legend=NA) +
  facet_wrap(.~lyr,nrow = 4,ncol = 3)+
  labs(fill = "Water proportion")+
  #scale_fill_gradientn(colours = paletteer_c("grDevices::Zissou 1", 30) ,na.value='white')+
  scale_fill_gradient(low = "lightblue", high = "darkblue", na.value = "white")+
  theme(
    panel.background = element_rect(fill = "white"),#背景设置
    strip.background.x = element_rect(fill = "white"), 
    # legend.position = c(1.2,0.2),#设置图例与主图距离
  )


#allDf <- fread(paste0(basePath,'allDf786_reclass.csv'))
globalCountry <- vect(world.map) 

speciesPixelNumPath <- list.files('/root/autodl-tmp/result/singleModel_water_xgb',pattern = '.tif',full.names = T)
spName <- basename(speciesPixelNumPath) %>% str_sub(.,1,-5)
AnserPath <- speciesPixelNumPath[spName=='Anser cygnoides']

AnserNum <- rast(AnserPath) %>% mask(.,globalCountry)
names(AnserNum)
length(names(AnserNum))
plot(AnserNum)

#用原来的分布
AnserMaxent<-rast('/root/autodl-tmp/result/model0.5/Anser cygnoides/maxent_1/probability_1.0.tif')

plotData <- AnserMaxent[[c('X1','X2','X3','X4','X5','X6','X7','X8','X9','X10','X11')]]
names(plotData) <- c('Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov')
ggplot() +
  geom_spatraster(data = plotData,show.legend=NA) +
  facet_wrap(.~lyr,nrow = 4,ncol = 3)+
  labs(fill = "Species numbers")+
  scale_fill_gradientn(colours = paletteer_c("grDevices::Zissou 1", 30) ,na.value='white')+
  theme(
    panel.background = element_rect(fill = "white"),#背景设置
    strip.background.x = element_rect(fill = "white"), 
    # legend.position = c(1.2,0.2),#设置图例与主图距离
  )





#Fig S2 Seasonal analysis#########


library(broom)
library(dplyr)
datacityError<-read.csv('/root/autodl-tmp/zyresult/dataurbanSp4_cityName.csv')
mydata<-as.data.frame(datacityError)
head(mydata)
datawin <- dplyr::select(mydata, cityName , Jan:Dec, lon, lat, CV)  #lon, lat

#####城市City----------
#!!思路：写一个循环，第一层是逐个读取每个城市的数据，第二层是对某个城市的12个月数据进行如下处理：取出1-6月的数据作为a1，接着2-7月数据作为a2，以此类推，最后是12，1-5月数据。对得到的这组a1,,,,,an数据进行方差分析，得出结果。
cities <- unique(mydata$cityName)
length(cities)

anova_winresults <- data.frame() # 创建一个空的数据框，用于存储方差分析的结果

# 第一层循环，逐个读取每个城市的数据
for (i in 1:318) {
  # 提取第i个城市的数据，只保留1-12月的丰富度数据
  city_data <- datawin[i, 2:13]
  city_data
  # 创建一个空的数据框，用来存储每个城市的a1, a2, ..., an数据
  city_a <- data.frame()
  city_a <- rbind(city_a, rep(0, 6))
  names(city_a)<-c("x1","x2","x3","x4","x5","x6")
  # 第二层循环，对某个城市的12个月数据进行如下处理
  for (j in 1:12) {
    # 取出第j到第j+5个月的数据，作为aj
    aj_index <- sapply(j:(j+5), function(x) x %% 12)   #超过12的模运算
    aj_index[aj_index == 0] <- 12
    aj <- city_data[, aj_index]
    names(aj)<-names(city_a)
    # 将aj添加到city_a数据框中
    city_a <- rbind(city_a, aj)
  }
  
  city_a <- city_a[-1, ]  #去掉一开始设置的第一行0000
  city_a_vec <- as.vector(t(city_a))
  city_a_group <- rep(1:12, each = 6)
  model_i  <- aov(city_a_vec ~ city_a_group)
  
  anova_i <- tidy(model_i, conf.int = TRUE)
  anova_i$cityName <- datawin$cityName[i]
  anova_winresults <- bind_rows(anova_winresults, anova_i)
}

# 查看results列表的结构
str(anova_winresults)


# 选出显著的
names(anova_winresults)<-c("term","df","sumsq","meansq","statistic","p_value","cityName")
anova_winresults$p_value
sig_wcities <- anova_winresults %>%
  dplyr::filter(`p_value` < 0.01) %>%
  distinct(cityName) %>%
  pull(cityName)

#sig_cities
length(sig_wcities)  #p0.05~270  ;  p0.01~231

#不显著的：
nosig_wcities <- anova_winresults %>% 
  dplyr::filter(p_value > 0.01) %>% 
  distinct(cityName) %>%  
  pull(cityName) 
#sig_cities
length(nosig_wcities)  #p0.05~270  ;  p0.01~231



#####赋予maxMonth-------
monthVariation <- fread(paste0(basePath,'AE_data/monthVariation.csv')); names(monthVariation)<-c("turn","direction","maxmonth","maxvalue","lon","lat")
# monthVariation2<-select(monthVariation,c('lon','lat','maxmonth','maxvalue'))
# monthVariation2 <- na.omit(monthVariation2);head(monthVariation2)
# monthVariationR<-rast(monthVariation2)
monthVariationR<-rast('/root/autodl-tmp/zyresult/monthVariationR.tif')
datacityError<-read.csv('/root/autodl-tmp/zyresult/dataurbanSp4_cityName.csv')
mydata<-as.data.frame(datacityError)
head(mydata$cityName)

mydata_vect <- vect(mydata, geom=c("lon", "lat"), crs=crs(monthVariationR));plot(mydata_vect)
extracted_values <- terra::extract(monthVariationR, mydata_vect) %>%as.data.frame()   #提取值
mydata$maxmonth <- extracted_values$maxmonth
mydata$maxvalue <- extracted_values$maxvalue
names(mydata)



#b上---------------------
#cities <- c("Shanghai", "Moscow", "New York", "Sao Paulo", "Abuja", "Sydney")
#mydata2 <- mydata[mydata$cityName %in% cities, ]

cities_geom <- c(29, 155, 26, 266, 171, 313, 114, 306)  #代表型城市
mydata2 <- mydata[mydata$geom %in% cities_geom, ]
#mydata2$maxmonth[mydata2$cityName == "Nairobi"] <- 1

#mydata2转换为长数据
mydata2_long <- pivot_longer(mydata2, cols = Jan:Dec, names_to = "month", values_to = "richness")
str(mydata2_long)
mydata2_long$month_num <- as.numeric(match(mydata2_long$month, month.abb)) ;mydata2_long$month_num# 创建月份序号变量


names(mydata2_long)

ggplot(mydata2_long # 筛选出有显著差异的城市的数据
       , aes(x = month_num, y = richness, group = lat, color = lat)) +   #Entropy_sum
  geom_line(size = 0.5, linetype = "solid",na.rm = FALSE) +
  scale_x_continuous(limits = c(1, 12), breaks = c(1:12)) +
  scale_y_continuous(limits = c(0, 140), breaks = seq(0,140,by=10)) +
  scale_color_gradient(low = "green", high = "darkblue") +  #连续变量设置颜色
  #scale_color_brewer(palette = "Set1")  #字符类型设置颜色
  labs(x = "Month", y = "Richness",title = "Significent")


ggplot(mydata2_long, aes(x = month_num, y = richness, group = cityName, color = cityName)) +
  geom_line(size = 0.5, linetype = "solid", na.rm = FALSE) +
  facet_wrap(~cityName, scales = 'free_y') + # 使用facet_wrap来创建分面
  scale_x_continuous(limits = c(1, 12), breaks = c(1:12)) +
  scale_y_continuous() +
  scale_color_manual(values = rainbow(n = length(unique(mydata2_long$cityName)))) + # 为每个城市分配不同的颜色
  labs(x = "Month", y = "Richness", title = "Richness by City and Month") +
  theme_minimal() +
  theme(legend.position = "none") # 隐藏图例

# 使用facet_grid创建图
library(RColorBrewer)

#自定义颜色
colors <- c("1" = "#F8766D", "2" = "#DE8C00", "3" = "#B79F00", "4" = "#7CAE00", "5" = "#00BA38",
            "6" = "#00C08B", "7" = "#00BFC4", "8" = "#00B4F0", "9" = "#619CFF", "10" = "#C77CFF",
            "11" = "#F564E3", "12" = "#FF64B0")   #"NA"='#7F7F7F'   "10" = "#C77CFF",

mydata2_long$cityName <- factor(mydata2_long$cityName, levels = c("Calgary", "Moscow", "Houston", "Guangzhou", "Kampala", "Adelaide", "Paris", "Makassar"))

# 可以在ggplot中使用这个颜色向量
ggplot(mydata2_long, aes(x = month_num, y = richness, color = as.factor(maxmonth), fill= as.factor(maxmonth))) +
  geom_line(size = 0.5, linetype = "solid", na.rm = FALSE) +
  #facet_grid(cityName ~ .) +   # 设置每个分面的y轴范围为自由 是  #, scales = 'free_y'
  facet_wrap(~cityName, ncol = 2) +   # 使用facet_wrap并设置每行两个图表
  geom_area(position = 'stack') +
  scale_color_manual(values = colors) +
  scale_fill_manual(values = colors) +   #alpha(colors, 0.5)
  #scale_y_continuous(limits = c(0, 150)) + # 设置y轴起始点为50
  theme_minimal() +
  scale_x_continuous(limits = c(1, 12),breaks = seq(1,12,2))+
  labs(x = "Month", y = "Species Richness")+
  theme(panel.grid = element_blank(), 
        panel.background = element_rect(fill = 'transparent', color = "black"), 
        legend.position = 'none',
        axis.text = element_text(size = 10),
        text = element_text(size = 10),
        #strip.text.y = element_text(size = 12, angle = 0) # 将分面标签文本横向展示
        strip.text = element_text(size = 12, angle = 0) # 将分面标签文本横向展示
  ) 



#a----------
coast <- ne_coastline(scale = "small", returnclass = "sf") %>% vect()
crs <- '+proj=longlat +datum=WGS84'

ggplot() +
  geom_raster(data = monthVariation, aes(x = lon, y = lat, fill = factor(maxmonth))) +
  geom_point(data = mydata2, aes(x = lon, y = lat), shape = 1, size = 4, color = "black", stroke = 1) +
  geom_spatvector(data=coast,fill=NA)+coord_sf(crs = crs,xlim=c(-160,165),ylim=c(-56,90))+
  geom_text(data = mydata2, aes(x = lon, y = lat, label = cityName), nudge_x = 0, nudge_y = -5, check_overlap = TRUE, size = 5) +
  scale_fill_discrete(labels = c("Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec")) +
  labs(x = NULL, y = NULL,) +
  theme_bw()+
  theme(
    #axis.text  = 'none',
    plot.title = element_text(hjust=0.5),
    # axis.line = element_line(color = 'black'),
    # panel.background = element_rect(fill = "white"),#背景设置
    legend.title = element_blank(),
    legend.position = "none"
  )

df<-as.data.frame(monthVariation)
#环形饼图
# install.packages("dplyr")
library(dplyr)
# 计算各月份的百分比
df <- df %>%
  group_by(maxmonth) %>%
  summarise(count = n()) %>%
  mutate(percentage = count / sum(count) * 100)

# install.packages("ggplot2")
library(ggplot2)
# 绘制环形饼图
ggplot(df, aes(x = "", y = percentage, fill = factor(maxmonth))) +
  geom_bar(stat = "identity", width = 1) +
  coord_polar(theta = "y") +
  geom_text(aes(label = paste0(round(percentage, 1), "%")), position = position_stack(vjust = 0.5)) +
  labs(x = NULL, y = NULL, fill = "Maxmonth") +
  scale_fill_brewer(palette = "Set3") +
  theme_void()


ggplot(df, aes(x = "", y = percentage, fill = factor(maxmonth))) +
  geom_bar(stat = "identity", width = 1) +
  coord_polar(theta = "y") +
  geom_text(aes(label = paste0(round(percentage, 1), "%")), position = position_stack(vjust = 0.5)) +
  labs(x = NULL, y = NULL, fill = "Maxmonth") +
  scale_fill_discrete() +
  theme_void()


#c--------
head(monthVariation)
monthVariation$season <- ifelse(monthVariation$lat > 0 & monthVariation$maxmonth %in% c(5, 6, 7, 8), "Breeding Season Risk (NH)",
                                ifelse(monthVariation$lat > 0 & monthVariation$maxmonth %in% c(11, 12, 1, 2), "Overwinter Season Risk (NH)",
                                       ifelse(monthVariation$lat < 0 & monthVariation$maxmonth %in% c(11, 12, 1, 2), "Overwinter Season Risk (SH)",
                                              "Other")))
monthVariation$season <- factor(monthVariation$season, levels = c("Breeding Season Risk (NH)", "Overwinter Season Risk (NH)", "Overwinter Season Risk (SH)", "Other"))

colors_season <- c("Breeding Season Risk (NH)" = "#fd8a8a",
                   "Overwinter Season Risk (NH)" = "#71959a", 
                   "Overwinter Season Risk (SH)" = "#73a8aa",
                   "Other" = "lightgrey")   #"NA"='#7F7F7F'   "10" = "#f49cba",


coast <- ne_coastline(scale = "small", returnclass = "sf") %>% vect()
crs <- '+proj=longlat +datum=WGS84'

ggplot() +
  geom_raster(data = monthVariation, aes(x = lon, y = lat, fill = season)) +
  #geom_point(data = mydata2, aes(x = lon, y = lat), shape = 1, size = 4, color = "black", stroke = 1) +
  geom_spatvector(data=coast,fill=NA)+coord_sf(crs = crs,xlim=c(-160,165),ylim=c(-56,90))+
  #geom_text(data = mydata2, aes(x = lon, y = lat, label = cityName), nudge_x = 0, nudge_y = -5, check_overlap = TRUE, size = 5) +
  scale_fill_manual(values=colors_season) +
  labs(x = NULL, y = NULL,) +
  theme_bw()+
  theme(
    #axis.text  = 'none',
    plot.title = element_text(hjust=0.5),
    # axis.line = element_line(color = 'black'),
    # panel.background = element_rect(fill = "white"),#背景设置
    legend.title = element_blank(),
    legend.position = c(0.13, 0.15),#设置图例与主图距离
    legend.key.width = unit(1,'cm'), #图例宽度
    legend.key.height = unit(1,'cm'),
    legend.text = element_text(size = 16),
  )





#b下-----------
city_northsummer<-c(29, 33, 40, 45, 56, 60, 72, 81, 144, 149, 151, 155, 158, 164, 172, 179, 189, 192, 241, 242, 248, 249)
city_northwinter<-c(4, 5, 7, 22, 24, 26, 27, 32, 37, 43, 46, 48, 49, 50, 110, 113, 116, 134, 138, 163, 167, 171, 175, 185, 187, 191, 198, 199, 202, 206, 207, 208, 209, 210, 213, 214, 215, 217, 221, 224, 226, 228, 237, 253, 260, 266, 271, 273, 277, 281, 283, 286, 287, 288, 291, 294, 295, 298, 300, 301, 302, 303, 304, 305, 307, 311, 312)
city_southsummer<-c(87, 90, 154, 159, 165, 166, 173, 176, 290, 297, 313, 171)
city_stable<-c(98, 100, 103, 104, 105, 108, 109, 114, 115, 117, 119, 120, 122, 123, 129, 133, 136, 139, 140, 153, 181, 183, 190, 233,38, 95, 264, 267, 270, 279, 306)

mydata_long <- pivot_longer(mydata, cols = Jan:Dec, names_to = "month", values_to = "richness")
str(mydata_long)
mydata_long$month_num <- as.numeric(match(mydata_long$month, month.abb)) ;mydata_long$month_num# 创建月份序号变量

mydata5_long <- mydata_long[mydata_long$geom %in% city_northsummer, ]
mydata6_long <- mydata_long[mydata_long$geom %in% city_northwinter, ]
mydata7_long <- mydata_long[mydata_long$geom %in% city_southsummer, ]
mydata8_long <- mydata_long[mydata_long$geom %in% city_stable, ]


library(ggplot2)
p5<-ggplot(mydata5_long # 筛选出北半球夏高城市的数据
           , aes(x = month_num, y = richness, group = lat)) +   #Entropy_sum   
  geom_line(size = 0.5, linetype = "solid",na.rm = FALSE, color="#fd8a8a") +
  scale_x_continuous(limits = c(1, 12), breaks = c(1:12)) +
  #scale_y_continuous(limits = c(0, 120), breaks = seq(0,120,by=10)) +
  #scale_color_gradient(low = "#feffe8", high = "#fd8a8a") +  #连续变量设置颜色
  #scale_color_brewer(palette = "Set1")  #字符类型设置颜色
  labs(x = NULL, y = "Species Richness",title = "Breeding Season Risk (NH)")+
  theme_bw()+
  theme(axis.text = element_text(size = 10),
        text = element_text(size = 10),
        plot.title = element_text(hjust=0.5))

p6<-ggplot(mydata6_long # 筛选出北半球冬高城市的数据
           , aes(x = month_num, y = richness, group = lat)) +   #Entropy_sum
  geom_line(size = 0.5, linetype = "solid",na.rm = FALSE, color="#71959a") +
  scale_x_continuous(limits = c(1, 12), breaks = c(1:12)) +
  #scale_y_continuous(limits = c(0, 120), breaks = seq(0,120,by=10)) +
  #scale_color_gradient(low = "#feffe8", high = "#71959a") +  #连续变量设置颜色
  #scale_color_brewer(palette = "Set1")  #字符类型设置颜色
  labs(x = NULL, y = "Species Richness",title = "Overwinter Season Risk (NH)")+
  theme_bw()+
  theme(        axis.text = element_text(size = 10),
                text = element_text(size = 10),
                plot.title = element_text(hjust=0.5))

p7<-ggplot(mydata7_long # 筛选出南半球高城市的数据
           , aes(x = month_num, y = richness, group = lat)) +   #Entropy_sum
  geom_line(size = 0.5, linetype = "solid",na.rm = FALSE, color="#73a8aa") +
  scale_x_continuous(limits = c(1, 12), breaks = c(1:12)) +
  #scale_y_continuous(limits = c(0, 120), breaks = seq(0,120,by=10)) +
  #scale_color_gradient(low = "#feffe8", high = "#f7a156") +  #连续变量设置颜色
  #scale_color_brewer(palette = "Set1")  #字符类型设置颜色
  labs(x = NULL, y = "Species Richness",title = "Overwinter Season Risk (SH)")+
  theme_bw()+
  theme(        axis.text = element_text(size = 10),
                text = element_text(size = 10),
                plot.title = element_text(hjust=0.5))

p8<-ggplot(mydata8_long # 筛选出南半球高城市的数据
           , aes(x = month_num, y = richness, group = lat)) +   #Entropy_sum
  geom_line(size = 0.5, linetype = "solid",na.rm = FALSE, color="grey") +
  scale_x_continuous(limits = c(1, 12), breaks = c(1:12)) +
  #scale_y_continuous(limits = c(0, 120), breaks = seq(0,120,by=10)) +
  #scale_color_gradient(low = "#feffe8", high = "lightgrey") +  #连续变量设置颜色
  #scale_color_brewer(palette = "Set1")  #字符类型设置颜色
  labs(x = 'Month', y = "Species Richness",title = "Other")+
  theme_bw()+
  theme(        axis.text = element_text(size = 10),
                text = element_text(size = 10),
                plot.title = element_text(hjust=0.5))

library(patchwork)
p5 + p6 + p7 + p8 + plot_layout(ncol = 1)






#——————————————————————————————————————————————————————————-----------------

#Fig S3 CV of SR&GPP#########
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


###a 分南北颜色分各界形状mian图-------
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
  stat_poly_eq(aes(label = paste(..eq.label.., ..rr.label..,..p.value.label.., sep = "~~")), formula = y ~ x, parse = TRUE) + # 添加公式和R^2
  labs(x = "GPP CV", y = "Species Richness CV") +
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





###b 大界分面展示supply图-------
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
  stat_poly_eq(aes(label = paste(..eq.label.., ..rr.label.., ..p.value.label.., sep = "~~")), formula = y ~ x, parse = TRUE,size=2) + # 添加公式和R^2
  # stat_poly_eq(aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")), 
  #              label.x.npc = "right", label.y.npc = 60,
  #              coef.digits = 3, rr.digits = 2, parse = TRUE, size = 3, color = "black")+ # 绘制线性回归方程和决定系数
  labs(x = "GPP CV", y = "Species Richness CV") +
  theme(
    panel.background = element_rect(fill = "white"),
    panel.border = element_rect(color = "black", linewidth = 1, fill = NA),
    panel.grid = element_blank(),
    legend.position = "none",
    legend.title = element_blank(),
    plot.title = element_text(size = 14, hjust = 0.5),
    axis.title = element_text(size = 12)
  )

#——————————————————————————————————————————————————————————-----------------
#Fig S4 分类accuracy#########

#1.分毒株
result_st2 <- fread('/root/autodl-tmp/humPoulResult/data/result_different_Serotype.csv')
ggplot(data = result_st2, aes(x = label, y = recall, fill = label)) +
  geom_col(position = 'dodge2') +
  theme_bw() +  # 使用黑白主题作为基础
  labs( y = "Accuracy") +  
  scale_fill_brewer(palette = "Pastel1") +
  ylim(c(0, 0.9))+
  theme(
    axis.title.x = element_blank(),
    text = element_text(size=18),
    plot.title = element_text(hjust = 0.5),  #
    legend.position = "none"  
  )
result_st3 <- mutate(result_st2, df_label = cut(cases, breaks = c(1500, 2000, 10000, 40000, 50000),
                                                     labels = c("1500 - 2000", "5000 - 6000", "35000 - 40000", "> 50000")))
size_values <- c("1500 - 2000" = 4, "5000 - 6000" = 8, "35000 - 40000" = 14, "> 50000" =20)
result_st3$Group<-c("Pathogenicity","Pathogenicity","Subtype","Subtype")
#result_st3 <- result_st3 %>%  mutate(recall = percent(recall))


# 假设result_st2是您的数据框，其中包含病例数量的列名为cases
result_st3$Group<-factor(result_st3$Group,levels = c("Subtype","Pathogenicity"))
result_st3$label<-factor(result_st3$label,levels = c("HPAI","LPAI","H5N1","H9N2"))

p1 <- ggplot(data = result_st3, aes(x = label, y = recall)) +
  geom_point(aes(size = df_label, color = Group), alpha = 0.5) +  # 使用df_label列来决定气泡大小
  geom_hline(aes(yintercept = 0.732), color = "grey", linetype = "dashed", alpha = 1) +
  scale_size_manual(values = size_values) +  # 手动设置气泡大小的值
  scale_color_manual(values = c("orange", "red")) +  # 手动设置颜色的值
  ylim(c(0.6, 0.85)) +
  labs(size = "Number of Cases", y = "Recall") +  # 添加气泡大小的图例标题
  theme_bw() +
  theme(
    axis.text = element_text(size = 12),
    axis.title.y = element_text(size = 16),
    axis.title.x = element_blank()
  )

p1

#2.分感染对象
`%notin%` <- Negate(`%in%`)
result_type2 <- fread('/root/autodl-tmp/humPoulResult/data/result_different_type.csv')
#result_type2$df[7]<-sum(result_type2$df)

#result_type2[type == "Captive", type := "Wild"]

# 计算合并后的新值
# result_type2 <- result_type2[, .(df = sum(df),
#                                  TP = sum(TP),
#                                  FN = sum(FN),
#                                  FP = sum(FP),
#                                  TN = sum(TN),
#                                  accuracy = mean(accuracy)), by = .(type)]

#result_type2 <- result_type2 %>%  mutate(recall = percent(recall))
result_type2 <- mutate(result_type2, df_label = cut(cases, breaks = c(0, 100, 500, 3000, 15000, 40000),
                                                labels = c("< 100", "100 - 500","2000 - 3000","10000 - 15000", "> 35000")))
size_values <- c("< 100" = 1, "100 - 500" = 3, "2000 - 3000" = 6, "10000 - 15000"=10, "> 35000" =14)
# result_type2 <- result_type2 %>% 
#   rename(Total = overall)

result_type2$type<-factor(result_type2$type,levels = c("Wild bird","Poultry","Wild mammal","Domestic mammal","Human affected","Environmental sample"))
p2<-ggplot(data=result_type2)+   #subset(result_type2,result_type2$type%notin%c('Captive_mammal',''))
  geom_point(aes(x = type, y = recall,size = df_label) , color="darkblue",alpha=0.5) + 
  #geom_col(aes(type,accuracy,fill=type),position = 'dodge2')+
  geom_hline(aes(yintercept =0.732),color='blue',linetype = "dashed",alpha=0.5)+
  scale_size_manual(values = size_values) +
  ylim(c(0.65,0.9))+
  labs(size = "Number of Cases",y="Recall",x=NULL) + 
  theme_bw()+
  #geom_text(aes(type,accuracy,label=df),size=4,vjust=-0.5)+
  theme(
    axis.text.y = element_text(size=12),
    axis.text.x = element_text(angle = 45,hjust=1,size=12),   #angle = 45,,hjust=1
    axis.title.y = element_text(size=16)
  )
p2

p1+p2

#Fig S4 H5N1#########
#a##########
library(terra)
library(data.table)
library(dplyr)

overEntropy <- rast(paste0(basePath,'AE_data/AE.tif'))%>% mask(globalCountry)
global(overEntropy,quantile,probs=seq(0, 1, 0.05),na.rm=T)


`%notin%` <- Negate(`%in%`)
library(lubridate)
outBreak1 <- fread('/root/autodl-tmp/YANZHENG/point/allData.csv') %>% subset(Diagnosis.status=='Confirmed'&Animal.type%in%c('Domestic','Wild'))

outBreak1$label <- str_extract(outBreak1$Serotype,'HPAI|LPAI')
outBreak1$h_label <- str_extract(outBreak1$Serotype,'H[0-9]N[0-9]|H[0-9]')
outBreak1 <- subset(outBreak1,outBreak1$h_label%notin%c('H9N2','H5N6'))   #这一句是用来运行全部的

outBreak1 <- subset(outBreak1,outBreak1$h_label=='H5N1'&outBreak1$label=='HPAI'|outBreak1$label=='LPAI')  #这一句是用来仅运行H5N1-HPAI的

outBreak1$Longitude <- as.numeric(outBreak1$Longitude)

addGeom <- cellFromXY(globalRaster,outBreak1[,c('Longitude','Latitude')]) %>% 
  cbind(id=.,outBreak1)
thinData <- unique(data.table(addGeom),by='id') %>% dplyr::select(.,-id)
outBreak2 <- vect(thinData,geom=c('Longitude','Latitude'),crs=crs)


plotEntropy <- ifel(overEntropy>=4.19,1,0)
#加载海岸线数据
coast <- rnaturalearth::ne_coastline(scale = "small", returnclass = "sf")
crs <- '+proj=longlat +datum=WGS84'

ggplot() +
  geom_spatraster(data = plotEntropy) +
  geom_spatvector(data=coast,fill=NA)+
  #coord_sf(crs = crs,xlim=c(-160,180),ylim=c(-56,90))+
  #geom_spatvector(data=outBreak2,size=2, shape = 4,color='#ff4f4c',fill='#f26c6a', alpha = 0.8,stroke = 0.3)+   #shape = 21是圆圈，shape = 4是×
  geom_spatvector(data = outBreak2, aes(color = label), size = 3, shape = 4, alpha = 0.8, stroke = 0.3) +
  scale_color_manual(values = c("HPAI" = "#370d51", "LPAI" = "#e040fb")) +
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

#b##########
# `%notin%` <- Negate(`%in%`)
# library(lubridate)
# #outBreak <- fread('/root/autodl-tmp/YANZHENG/point/allData.csv')
# outBreak1 <- fread('/root/autodl-tmp/YANZHENG/point/allData.csv') %>% subset(Diagnosis.status=='Confirmed'&Animal.type%in%c('Domestic','Wild'))
# # outBreak <- subset(outBreak,outBreak$Diagnosis.status=='Confirmed')
# # r <- subset(outBreak,outBreak$Species=='Unspecified Bird'&outBreak$Animal.type=='	
# # Domestic')
# outBreak1$label <- str_extract(outBreak1$Serotype,'HPAI|LPAI')
# outBreak1$h_label <- str_extract(outBreak1$Serotype,'H[0-9]N[0-9]|H[0-9]')
# # r <- outBreak[outBreak$Longitude>60&outBreak$Longitude<120&outBreak$Latitude>35&outBreak$Latitude<60,]
# #outBreak1 <- subset(outBreak1,outBreak1$h_label%in%c('H9N2','H5N6'))
# outBreak1 <- subset(outBreak1,outBreak1$h_label%notin%c('H9N2','H5N6'))
# outBreak1 <- subset(outBreak1,outBreak1$h_label=='H5N1'&outBreak1$label=='HPAI'|outBreak1$label=='LPAI')  #这一句是用来仅运行H5N1-HPAI的
# 
# 
# outBreak1$Longitude <- as.numeric(outBreak1$Longitude)
# 
# # outBreak <- subset(outBreak,outBreak$year%in%2021:2023)
# addGeom <- cellFromXY(globalRaster,outBreak1[,c('Longitude','Latitude')]) %>% 
#   cbind(id=.,outBreak1)
# thinData <- unique(data.table(addGeom),by='id') %>% dplyr::select(.,-id)
# outBreak2 <- vect(thinData,geom=c('Longitude','Latitude'),crs=crs)
# # r <- terra::extract(v,outBreak2)
# 
# overEntropy <- rast(paste0(basePath,'AE_data/AE.tif'))%>% crop(globalCountry,mask=T)
# # global(overEntropy,quantile,probs=seq(0, 1, 0.05),na.rm=T)
# 
# # overEntropy2 <- ifel(overEntropy>4.152,overEntropy,NA)
# 
# 
# df <- terra::extract(overEntropy,outBreak2,mean,na.rm=T,cells=T)%>% as.data.frame()
# allDf <- terra::as.data.frame(overEntropy,na.rm=T,cells=T)
# 
# library(pROC)
# 
# df$label <- 'good'
# rocDf <- left_join(allDf,df[,c('cell','label')])
# rocDf$label <- ifelse(is.na(rocDf$label),'poor',rocDf$label)
# dfroc1<- roc(rocDf$label, rocDf$sum)
# plot(dfroc1,col="red",#颜色
#      legacy.axes=T,#y轴格式更改
#      print.auc=TRUE,#显示AUC面积
#      expand=c(0,0),
#      print.thres=TRUE,#添加截点和95%CI
#      grid=c(0.2,0.2),grid.col=c("lightgrey","lightgrey"),
#      #cex.text=1.5,
#      cex.main=1,  # 主标题字体放大1.5倍
#      cex.sub=1,   # 副标题字体放大1.5倍
#      cex.axis=1,  # 坐标轴刻度字体放大1.5倍
#      cex.lab=1)#网格线设置
# 
# 
# 
# threshold <- 4.107
# 
# TP <- sum(df$sum>=threshold,na.rm = T)  ;TP
# FN <- sum(df$sum<threshold,na.rm = T)   ;FN
# FP <- sum(allDf$sum>=threshold,na.rm = T)-TP  ;FP
# TN <- nrow(allDf)-TP-FN-FP  ;TN
# 
# accuracy <- (TP+TN)/(TP+TN+FN+FP)  ;accuracy
# precise <- (TP)/(TP+FP)  ;precise
# recall <- (TP)/(TP+FN)  ;recall
# FI <- 2*((recall*precise)/(recall+precise))  ;FI

#H5N1的最新的ROC##############
`%notin%` <- Negate(`%in%`)
crs <- '+proj=longlat +datum=WGS84'
globalRaster <- rast(vals=1:259200,nrows=360, ncols=720,xmin=-180, xmax=180,ymin=-90, ymax=90,crs=crs)
globalSHP <- vect('/root/autodl-tmp/worldBorder/ne_10m_admin_0_countries.shp')
# calculate  threshold
outBreakData <- fread('/root/autodl-tmp/YANZHENG/point/allData.csv')%>% subset(Diagnosis.status=='Confirmed')

outBreakData$label <- str_extract(outBreakData$Serotype,'HPAI|LPAI')
outBreakData$h_label <- str_extract(outBreakData$Serotype,'H[0-9]N[0-9]|H[0-9]')
outBreakData <- subset(outBreakData,outBreakData$h_label =='H5N1')
#outBreakData <- subset(outBreakData,outBreakData$h_label %in% c('H5N1','H5N2','H5N8','H5','H5N3','H5N6','H5N5','H5N9','H5N4'))


outBreakData2 <- na.omit(outBreakData[,c('Country','Longitude','Latitude')])
addGeom <- cellFromXY(globalRaster,outBreakData2[,c('Longitude','Latitude')]) %>% 
  cbind(id=.,outBreakData2)
thinData <- unique(data.table(addGeom),by='id') %>% dplyr::select(.,-id)
breakVect <- vect(thinData,geom=c('Longitude','Latitude'),crs=crs)
breakRast <- rasterize(breakVect,globalRaster)

TP_points <- vect(outBreakData2,geom=c('Longitude','Latitude'),crs=crs)
AE <- rast('/root/autodl-tmp/humPoulResult/data/AE_data/AE.tif')   #
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


threshold <- dfroc1$auc

WAE<- rast(paste0(basePath,'AE_data/AE_Waterfowls.tif'))%>% mask(globalCountry)
plotWAE <- ifel(WAE>=4.19,1,0)   #2.559

coast <- rnaturalearth::ne_coastline(scale = "small", returnclass = "sf")
crs <- '+proj=longlat +datum=WGS84'
ggplot() +
  geom_spatraster(data = plotWAE) +
  geom_spatvector(data=coast,fill=NA)+
  #coord_sf(crs = crs,xlim=c(-160,180),ylim=c(-56,90))+
  geom_spatvector(data=breakVect,size=2, shape = 4,color='#ff4f4c',fill='#f26c6a', alpha = 0.8,stroke = 0.3) +   #shape = 21是圆圈，shape = 4是×
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




#试试只用水禽WAE预测##############
`%notin%` <- Negate(`%in%`)
crs <- '+proj=longlat +datum=WGS84'
globalRaster <- rast(vals=1:259200,nrows=360, ncols=720,xmin=-180, xmax=180,ymin=-90, ymax=90,crs=crs)
globalSHP <- vect('/root/autodl-tmp/worldBorder/ne_10m_admin_0_countries.shp')
# calculate  threshold
outBreakData <- fread('/root/autodl-tmp/YANZHENG/point/allData.csv')%>% subset(Diagnosis.status=='Confirmed')

outBreakData$label <- str_extract(outBreakData$Serotype,'HPAI|LPAI')
outBreakData$h_label <- str_extract(outBreakData$Serotype,'H[0-9]N[0-9]|H[0-9]')
outBreakData <- subset(outBreakData,outBreakData$h_label =='H5N1')
#outBreakData <- subset(outBreakData,outBreakData$h_label %in% c('H5N1','H5N2','H5N8','H5','H5N3','H5N6','H5N5','H5N9','H5N4'))


outBreakData2 <- na.omit(outBreakData[,c('Country','Longitude','Latitude')])
addGeom <- cellFromXY(globalRaster,outBreakData2[,c('Longitude','Latitude')]) %>% 
  cbind(id=.,outBreakData2)
thinData <- unique(data.table(addGeom),by='id') %>% dplyr::select(.,-id)
breakVect <- vect(thinData,geom=c('Longitude','Latitude'),crs=crs)
breakRast <- rasterize(breakVect,globalRaster)

TP_points <- vect(outBreakData2,geom=c('Longitude','Latitude'),crs=crs)
AE <- rast('/root/autodl-tmp/humPoulResult/data/AE_data/AE_Waterfowls.tif')   #
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

WAE<- rast(paste0(basePath,'AE_data/AE.tif'))%>% mask(globalCountry)
plotWAE <- ifel(WAE>=2.559,1,0)   #2.559

coast <- rnaturalearth::ne_coastline(scale = "small", returnclass = "sf")
crs <- '+proj=longlat +datum=WGS84'
ggplot() +
  geom_spatraster(data = plotWAE) +
  geom_spatvector(data=coast,fill=NA)+
  #coord_sf(crs = crs,xlim=c(-160,180),ylim=c(-56,90))+
  geom_spatvector(data=breakVect,size=2, shape = 4,color='#ff4f4c',fill='#f26c6a', alpha = 0.8,stroke = 0.3) +   #shape = 21是圆圈，shape = 4是×
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




#——————————————————————————————————————————————————————————-----------------
#Fig S5 WAE-Human & WAE-Poultry#########


#a : 人口和活动熵栅格二元--------
popd2015<-rast("/root/autodl-tmp/全球人口/GWP_v4/popd2015_30.tif")
Entropy<-rast(paste0(basePath,'AE_data/AE.tif'))%>% mask(globalCountry)
world.map <- ne_countries(returnclass = "sf") |>filter(continent != "Antarctica")
globalCountry <- vect(world.map) 

#预处理
popd2015_m <- trim(mask(popd2015, globalCountry))
Entropy_m <- trim(mask(Entropy, globalCountry))  ; plot(Entropy_m)   ##writeRaster(Entropy_m,"/root/autodl-tmp/zyresult/EntropySum_Contry.tif")

popd2015_r <- resample(popd2015_m, Entropy_m, method="bilinear")
Entropy_m <- trim(mask(Entropy_m, popd2015_r))
Entropy_e <- extend(Entropy_m, popd2015_r)

popd2015_df <- as.data.frame(popd2015_r, xy=TRUE)
Entropy_df <- as.data.frame(Entropy_e, xy=TRUE)


#画一下人口图
library(paletteer)
ggplot() +
  geom_spatraster(data = popd2015_r) +  #sp_cv
  #  geom_spatvector(data=coast,fill=NA)+coord_sf(crs = crs,xlim=c(-160,165),ylim=c(-56,90))+
  theme_bw()+
  #scale_fill_gradient(low = "grey", high = "red",na.value='white',values = c(0, 0.2, 1))
  scale_fill_gradientn(colours = paletteer_c("grDevices::Zissou 1", 30) ,na.value='NA',values = c(0, 0.01, 1))+
  labs(title =  "Population")+
  theme(
    plot.title = element_text(hjust=0.5),
    legend.title = element_text("Population"),
    legend.title.align = -10,
    legend.position = c(0.15, 0.05),#设置图例与主图距离
    legend.direction='horizontal',#图例水平放置vertical，垂直horizontal，
    legend.key.width = unit(1.2,'cm'), #图例宽度
    legend.key.height = unit(0.3,'cm')
  )

# 定义归一化函数
normalize <- function(x) {
  return((x - min(x, na.rm = TRUE)) / (max(x, na.rm = TRUE) - min(x, na.rm = TRUE)))
}


#####归一化处理
head(popd2015_df)
popd2015_df$popd2015_30_normalized <- normalize(popd2015_df$popd2015_30)

head(Entropy_df)
Entropy_df$sum_normalized <- normalize(Entropy_df$sum)

popentr_df <- left_join(popd2015_df, Entropy_df, by = c("x", "y"))
popentr_df <- na.omit(popentr_df)
head(popentr_df)
names(popentr_df)<-c("x","y","pop","pop_n","entr","entr_n")
max(popentr_df$pop_n)
max(popentr_df$entr_n)

#####(1)设置分位数和10*10颜色变量
# create 10 buckets for gini
quantiles_pop <- popentr_df %>%
  pull(pop_n) %>%
  quantile(probs = seq(0, 1, length.out = 11))
quantiles_pop

quantiles_entr <- quantile(popentr_df$entr_n, probs = seq(0, 1, length.out = 11), na.rm = TRUE)
quantiles_entr 


bivariate_color_scale <- tibble(
  "1 - 10" ="#00FFFF","2 - 10" = "#11E6FD", "3 - 10" ="#23CDFB","4 - 10" = "#35B4FA","5 - 10" ="#479BF8","6 - 10"= "#5883F6","7 - 10" ="#6A6AF5","8 - 10" ="#7C51F3","9 - 10" ="#8E38F1","10 - 10" ="#A020F0",
  "1 - 9" = "#1BFDFB", "2 - 9" = "#2AE4F7", "3 - 9" = "#3ACCF4", "4 - 9" = "#4AB4F0","5 - 9" = "#5A9CED","6 - 9" = "#6A83E9","7 - 9" = "#7A6BE6","8 - 9" = "#8A53E2","9 - 9" = "#9A3BDF","10 - 9" = "#AA23DC",
  "1 - 8" = "#36FCF7", "2 - 8" = "#44E4F1", "3 - 8" = "#52CCEC", "4 - 8" = "#60B5E7","5 - 8" = "#6E9DE2","6 - 8" = "#7C85DC","7 - 8" = "#8A6ED7","8 - 8" = "#9856D2","9 - 8" = "#A63ECD","10 - 8" = "#B527C8",
  "1 - 7" = "#51FBF3", "2 - 7" = "#5DE3EC", "3 - 7" = "#69CCE5", "4 - 7" = "#75B5DE","5 - 7" = "#819ED7","6 - 7" = "#8E86D0","7 - 7" = "#9A6FC9","8 - 7" = "#A658C2","9 - 7" = "#B241BB","10 - 7" = "#BF2AB5",
  "1 - 6" = "#6CFAEF", "2 - 6" = "#76E3E6", "3 - 6" = "#80CCDD", "4 - 6" = "#8BB6D5","5 - 6" = "#959FCC","6 - 6" = "#A088C3","7 - 6" = "#AA72BB","8 - 6" = "#B55BB2","9 - 6" = "#BF44A9","10 - 6" = "#CA2EA1",
  "1 - 5" = "#88F9EB", "2 - 5" = "#90E2E0", "3 - 5" = "#98CCD6", "4 - 5" = "#A1B6CB","5 - 5" = "#A9A0C1","6 - 5" = "#B289B7","7 - 5" = "#BA73AD","8 - 5" = "#C35DA2","9 - 5" = "#CB4798","10 - 5" = "#D4318E",
  "1 - 4" = "#A3F8E7", "2 - 4" = "#A9E2DA", "3 - 4" = "#B0CCCE", "4 - 4" = "#B7B7C2","5 - 4" = "#BDA1B6","6 - 4" = "#C48BAA","7 - 4" = "#CB769E","8 - 4" = "#D16092","9 - 4" = "#D84A86","10 - 4" = "#DF357A",
  "1 - 3" = "#BEF7E3", "2 - 3" = "#C2E1D5", "3 - 3" = "#C7CCC7", "4 - 3" = "#CCB7B9","5 - 3" = "#D1A2AB","6 - 3" = "#D58C9E","7 - 3" = "#DA7790","8 - 3" = "#DF6282","9 - 3" = "#E44D74","10 - 3" = "#E93867",
  "1 - 2" = "#D9F6DF", "2 - 2" = "#D9F6DF", "3 - 2" = "#D9F6DF", "4 - 2" = "#DBE1CF","5 - 2" = "#DFCCBF","6 - 2" = "#DFCCBF","7 - 2" = "#E2B8B0","8 - 2" = "#E5A3A0","9 - 2" = "#E88E91","10 - 2" = "#EE6572",
  "1 - 1" = "#F5F5DC", "2 - 1" = "#F5F5DC", "3 - 1" = "#F5F5DC", "4 - 1" = "#F6E0CA","5 - 1" = "#F7CCB9","6 - 1" = "#F7CCB9","7 - 1" = "#F7CCB9","8 - 1" = "#F8B8A8","9 - 1" = "#F9A496","10 - 1" = "#FA9085",
  #  "1 - 1" = "#F5F5DC", "2 - 1" = "#F6E0CA", "3 - 1" = "#F6E0CA", "4 - 1" = "#F7CCB9","5 - 1" = "#F7CCB9","6 - 1" = "#F8B8A8","7 - 1" = "#F9A496","8 - 1" = "#FA9085","9 - 1" = "#FB7C74","10 - 1" = "#FC6862",
  #"1 - 1" = "#F5F5DC", "2 - 1" = "#F6E0CA", "3 - 1" = "#F7CCB9", "4 - 1" = "#F8B8A8","5 - 1" = "#F9A496","6 - 1" = "#FA9085","7 - 1" = "#FB7C74","8 - 1" = "#FC6862","9 - 1" = "#FD5451","10 - 1" = "#FF4040",
  #"1 - 2" = "#D9F6DF", "2 - 2" = "#DBE1CF", "3 - 2" = "#DFCCBF", "4 - 2" = "#E2B8B0","5 - 2" = "#E5A3A0","6 - 2" = "#E88E91","7 - 2" = "#EB7981","8 - 2" = "#EE6572","9 - 2" = "#F15062","10 - 2" = "#F43C53",
  
) %>%
  gather("group", "fill")

#####(2)分配到不同的分位数组，并且给每个组分配一个颜色
# 使用cut函数将pop和entr的值分配到相应的分位数桶中
# 需要保证分位数唯一

popentr_df1 <- popentr_df %>%
  mutate(
    pop_quantiles = cut(pop_n, breaks = quantiles_pop, include.lowest = TRUE, labels = FALSE),
    entr_quantiles = cut(entr_n, breaks = quantiles_entr, include.lowest = TRUE, labels = FALSE)
  ) %>%
  # 将分位数转换为数字，并创建一个新的组合列
  mutate(
    group = paste(as.numeric(pop_quantiles), "-", as.numeric(entr_quantiles))
    #group = paste(as.numeric(entr_quantiles), "-", as.numeric(pop_quantiles))
  ) %>%
  # 假设 已经有了一个颜色表bivariate_color_scale，其中包含group和对应的颜色值
  left_join(bivariate_color_scale, by = "group") # 将颜色值与主数据框合并

#####(3)绘制地图
library(rnaturalearth)
coast <- ne_coastline(scale = "small", returnclass = "sf") %>% vect()
crs <- '+proj=longlat +datum=WGS84'
map<-ggplot(popentr_df1) +
  geom_tile(aes(x = x, y = y,fill = fill)) +
  scale_fill_identity() +
  geom_spatvector(data=coast,fill=NA)+coord_sf(crs = crs,xlim=c(-160,165),ylim=c(-56,90))+
  #geom_sf(data = canton_geo, fill = NA, color = "white", size = 0.5) +    # 为行政区划绘制边界
  labs(x = NULL, y = NULL,) +
  #theme_minimal()+
  theme_bw()+
  theme(
    legend.position="none"
  )
print(map)

#绘制图例
# 确保bivariate_color_scale是 的颜色表，并且包含group和对应的颜色值
# 分离group列为pop和entr两个独立的列，并转换为整数类型
bivariate_color_scale1 <- bivariate_color_scale %>%
  separate(group, into = c("pop", "entr"), sep = " - ", convert = TRUE) %>%
  mutate(
    pop = as.integer(pop),
    entr = as.integer(entr)
  )


# 创建图例
legend <- ggplot(data = bivariate_color_scale1) +
  geom_tile(mapping = aes(x = pop, y = entr, fill = fill)) +
  scale_fill_identity() +
  labs(x = "Higher population density ⟶️",   #population density
       y = "Higher entropy ⟶️") +   #entropy
  theme_void() +  # 使用theme_void()来移除多余的元素
  theme(
    axis.title.x = element_text(size = 16, angle = 0),  # 旋转X轴标题
    axis.title.y = element_text(size = 16, angle = 90),   # 旋转Y轴标题
    plot.margin = margin(t = 10, r = 10, b = 30, l = 10)  # 调整图表边距
  ) +
  # Set discrete values for y axis
  coord_fixed()  # 保持坐标比例固定
print(legend)


#b :  家禽和活动熵栅格二元-------------
# duck2015_r<-raster(duck2015_r)
# chic2015_r<-raster(chic2015_r)
# # 定义一个自定义函数，该函数将两个数值相加，并在其中一个为NA时保留另一个值
# sum_na <- function(x, y) {
#   return(ifelse(is.na(x), y, ifelse(is.na(y), x, x + y)))
# }
# # 使用overlay函数应用自定义函数
# poul_raster <- overlay(duck2015_r, chic2015_r, fun = sum_na)
# # 查看结果
# print(poul_raster)
# #poul_raster <- calc(duck2015_r + chic2015_r, fun = sum)
# writeRaster(poul_raster,"/root/autodl-tmp/zyresult/Poultry_duckchic.tif",overwrite=T)

poul2015<-rast("/root/autodl-tmp/zyresult/Poultry_duckchic.tif")
Entropy<-rast(paste0(basePath,'AE_data/AE.tif'))%>% mask(globalCountry)
hist(Entropy)
world.map <- ne_countries(returnclass = "sf") |>filter(continent != "Antarctica")
globalCountry <- vect(world.map) 

#预处理
poul2015_m <- trim(mask(poul2015, globalCountry))
Entropy_m <- trim(mask(Entropy, globalCountry))  ;plot(Entropy_m)

poul2015_r <- resample(poul2015_m, Entropy_m, method="bilinear");plot(poul2015_r)
Entropy_m <- trim(mask(Entropy_m, poul2015_r))
Entropy_e <- extend(Entropy_m, poul2015_r);plot(Entropy_e)

poul2015_df <- as.data.frame(poul2015_r, xy=TRUE)
Entropy_df <- as.data.frame(Entropy_e, xy=TRUE)

#对数据进行归一化处理，统一度量方便展示：
# 定义归一化函数
normalize <- function(x) {
  return((x - min(x, na.rm = TRUE)) / (max(x, na.rm = TRUE) - min(x, na.rm = TRUE)))
}
#####归一化处理
head(poul2015_df)   #查看数据列的名字为popd2015_30
names(poul2015_df)<-c("x","y","poul")
poul2015_df$poul2015_30_normalized <- normalize(poul2015_df$poul)
head(Entropy_df)   #查看数据列的名字为sum
Entropy_df$sum_normalized <- normalize(Entropy_df$sum)

poulentr_df <- left_join(poul2015_df, Entropy_df, by = c("x", "y"))
poulentr_df <- na.omit(poulentr_df)
head(poulentr_df)
names(poulentr_df)<-c("x","y","poul","poul_n","entr","entr_n")


# poulentr_df$color_n <- with(poulentr_df, pouln + entrn)
# 
# ggplot(poulentr_df, aes(x = x, y = y, fill = color_n)) +
#   geom_tile() +
#   scale_fill_gradientn(colours  = c("blue", "green", "yellow", "red")) +
#   theme_minimal() +
#   labs(fill = "pouln + entrn") +
#   coord_fixed()

#####(1)设置10个分位数
quantiles_poul <- poulentr_df %>%
  pull(poul_n) %>%
  quantile(probs = seq(0, 1, length.out = 11))
quantiles_poul


quantiles_entr <- poulentr_df %>%
  pull(entr_n) %>%
  quantile(probs = seq(0, 1, length.out = 11))

quantiles_entr 




#设置10*10颜色变量
bivariate_color_scale <- tibble(
  "1 - 10" ="#00FFFF","2 - 10" = "#11E6FD", "3 - 10" ="#23CDFB","4 - 10" = "#35B4FA","5 - 10" ="#479BF8","6 - 10"= "#5883F6","7 - 10" ="#6A6AF5","8 - 10" ="#7C51F3","9 - 10" ="#8E38F1","10 - 10" ="#A020F0",
  "1 - 9" = "#1BFDFB", "2 - 9" = "#2AE4F7", "3 - 9" = "#3ACCF4", "4 - 9" = "#4AB4F0","5 - 9" = "#5A9CED","6 - 9" = "#6A83E9","7 - 9" = "#7A6BE6","8 - 9" = "#8A53E2","9 - 9" = "#9A3BDF","10 - 9" = "#AA23DC",
  "1 - 8" = "#36FCF7", "2 - 8" = "#44E4F1", "3 - 8" = "#52CCEC", "4 - 8" = "#60B5E7","5 - 8" = "#6E9DE2","6 - 8" = "#7C85DC","7 - 8" = "#8A6ED7","8 - 8" = "#9856D2","9 - 8" = "#A63ECD","10 - 8" = "#B527C8",
  "1 - 7" = "#51FBF3", "2 - 7" = "#5DE3EC", "3 - 7" = "#69CCE5", "4 - 7" = "#75B5DE","5 - 7" = "#819ED7","6 - 7" = "#8E86D0","7 - 7" = "#9A6FC9","8 - 7" = "#A658C2","9 - 7" = "#B241BB","10 - 7" = "#BF2AB5",
  "1 - 6" = "#6CFAEF", "2 - 6" = "#76E3E6", "3 - 6" = "#80CCDD", "4 - 6" = "#8BB6D5","5 - 6" = "#959FCC","6 - 6" = "#A088C3","7 - 6" = "#AA72BB","8 - 6" = "#B55BB2","9 - 6" = "#BF44A9","10 - 6" = "#CA2EA1",
  "1 - 5" = "#88F9EB", "2 - 5" = "#90E2E0", "3 - 5" = "#98CCD6", "4 - 5" = "#A1B6CB","5 - 5" = "#A9A0C1","6 - 5" = "#B289B7","7 - 5" = "#BA73AD","8 - 5" = "#C35DA2","9 - 5" = "#CB4798","10 - 5" = "#D4318E",
  "1 - 4" = "#A3F8E7", "2 - 4" = "#A9E2DA", "3 - 4" = "#B0CCCE", "4 - 4" = "#B7B7C2","5 - 4" = "#BDA1B6","6 - 4" = "#C48BAA","7 - 4" = "#CB769E","8 - 4" = "#D16092","9 - 4" = "#D84A86","10 - 4" = "#DF357A",
  "1 - 3" = "#BEF7E3", "2 - 3" = "#C2E1D5", "3 - 3" = "#C7CCC7", "4 - 3" = "#CCB7B9","5 - 3" = "#D1A2AB","6 - 3" = "#D58C9E","7 - 3" = "#DA7790","8 - 3" = "#DF6282","9 - 3" = "#E44D74","10 - 3" = "#E93867",
  "1 - 2" = "#D9F6DF", "2 - 2" = "#DBE1CF", "3 - 2" = "#DFCCBF", "4 - 2" = "#E2B8B0","5 - 2" = "#E5A3A0","6 - 2" = "#E88E91","7 - 2" = "#EB7981","8 - 2" = "#EE6572","9 - 2" = "#F15062","10 - 2" = "#F43C53",
  "1 - 1" = "#F5F5DC", "2 - 1" = "#F6E0CA", "3 - 1" = "#F7CCB9", "4 - 1" = "#F8B8A8","5 - 1" = "#F9A496","6 - 1" = "#FA9085","7 - 1" = "#FB7C74","8 - 1" = "#FC6862","9 - 1" = "#FD5451","10 - 1" = "#FF4040",
  
) %>%
  gather("group", "fill")

#####(2)分配
poulentr_df1 <- poulentr_df %>%
  mutate(
    poul_quantiles = cut(poul_n, breaks = quantiles_poul, include.lowest = TRUE, labels = FALSE),
    entr_quantiles = cut(entr_n, breaks = quantiles_entr, include.lowest = TRUE, labels = FALSE)
  ) %>%
  # 将分位数转换为数字，并创建一个新的组合列
  mutate(
    group = paste(as.numeric(poul_quantiles), "-", as.numeric(entr_quantiles))
  ) %>%
  # 假设 已经有了一个颜色表bivariate_color_scale，其中包含group和对应的颜色值
  left_join(bivariate_color_scale, by = "group") # 将颜色值与主数据框合并



#####(3)绘制地图
library(rnaturalearth)
coast <- ne_coastline(scale = "small", returnclass = "sf") %>% vect()
crs <- '+proj=longlat +datum=WGS84'

map<-ggplot(poulentr_df1) +
  geom_tile(aes(x = x, y = y,fill = fill)) +
  scale_fill_identity() +
  geom_spatvector(data=coast,fill=NA)+coord_sf(crs = crs,xlim=c(-160,165),ylim=c(-56,90))+
  #geom_sf(data = canton_geo, fill = NA, color = "white", size = 0.5) +    # 为行政区划绘制边界
  labs(x = NULL, y = NULL,) +
  #theme_minimal()+
  theme_bw()+
  theme(
    legend.position="none"
  )
print(map)


# bivariate_color_scale是颜色表，并且包含group和对应的颜色值
# 分离group列为poul和entr两个独立的列，并转换为整数类型
bivariate_color_scale1 <- bivariate_color_scale %>%
  separate(group, into = c("poul", "entr"), sep = " - ", convert = TRUE) %>%
  mutate(
    poul = as.integer(poul),
    entr = as.integer(entr)
  )

#绘制图例
legend <- ggplot(data = bivariate_color_scale1) +
  geom_tile(mapping = aes(x = poul, y = entr, fill = fill)) +
  scale_fill_identity() +
  labs(x = "Higher poultry(duck&chicken) density ⟶️",
       y = "Higher entropy ⟶️") +
  theme_void() +  # 使用theme_void()来移除多余的元素
  theme(
    axis.title.x = element_text(size = 16, angle = 0),  # 旋转X轴标题
    axis.title.y = element_text(size = 16, angle = 90),   # 旋转Y轴标题
    plot.margin = margin(t = 10, r = 10, b = 30, l = 10)  # 调整图表边距
  ) +
  coord_fixed()  # 保持坐标比例固定
print(legend)




#new: 牲畜 牛cattle-------------
globalRaster <- rast(vals=1:259200,nrows=360, ncols=720,xmin=-180, xmax=180,ymin=-90, ymax=90,crs=crs)

cattle2015<-rast("/root/autodl-tmp/全球家禽/Cattle/5_Ct_2015_Da.tif")
Entropy<-rast("/root/autodl-tmp/humPoulResult/data/AE_data/AE.tif")

world.map <- ne_countries(returnclass = "sf") |>filter(continent != "Antarctica")
globalCountry <- vect(world.map) 

#预处理
cattle2015_m <- trim(mask(cattle2015, globalCountry))
Entropy_m <- trim(mask(Entropy, globalCountry))  ;plot(Entropy_m)

cattle2015_r <- resample(cattle2015_m, Entropy_m, method="bilinear");plot(cattle2015_r)
Entropy_m <- trim(mask(Entropy_m, cattle2015_r))
Entropy_e <- extend(Entropy_m, cattle2015_r);plot(Entropy_e)

cattle2015_df <- as.data.frame(cattle2015_r, xy=TRUE)
Entropy_df <- as.data.frame(Entropy_e, xy=TRUE)

#对数据进行归一化处理，统一度量方便展示：
# 定义归一化函数
normalize <- function(x) {
  return((x - min(x, na.rm = TRUE)) / (max(x, na.rm = TRUE) - min(x, na.rm = TRUE)))
}
#####归一化处理
head(cattle2015_df)   #查看数据列的名字为popd2015_30
names(cattle2015_df)<-c("x","y","cattle")
cattle2015_df$cattle2015_30_normalized <- normalize(cattle2015_df$cattle)
head(Entropy_df)   #查看数据列的名字为sum
Entropy_df$sum_normalized <- normalize(Entropy_df$sum)

cattleentr_df <- left_join(cattle2015_df, Entropy_df, by = c("x", "y"))
cattleentr_df <- na.omit(cattleentr_df)
head(cattleentr_df)
names(cattleentr_df)<-c("x","y","cattle","cattle_n","entr","entr_n")


# cattleentr_df$color_n <- with(cattleentr_df, cattlen + entrn)
# 
# ggplot(cattleentr_df, aes(x = x, y = y, fill = color_n)) +
#   geom_tile() +
#   scale_fill_gradientn(colours  = c("blue", "green", "yellow", "red")) +
#   theme_minimal() +
#   labs(fill = "cattlen + entrn") +
#   coord_fixed()

#####(1)设置10个分位数
quantiles_cattle <- cattleentr_df %>%
  pull(cattle_n) %>%
  quantile(probs = seq(0, 1, length.out = 11))
quantiles_cattle


quantiles_entr <- cattleentr_df %>%
  pull(entr_n) %>%
  quantile(probs = seq(0, 1, length.out = 11))

quantiles_entr 




#设置10*10颜色变量
bivariate_color_scale <- tibble(
  "1 - 10" ="#00FFFF","2 - 10" = "#11E6FD", "3 - 10" ="#23CDFB","4 - 10" = "#35B4FA","5 - 10" ="#479BF8","6 - 10"= "#5883F6","7 - 10" ="#6A6AF5","8 - 10" ="#7C51F3","9 - 10" ="#8E38F1","10 - 10" ="#A020F0",
  "1 - 9" = "#1BFDFB", "2 - 9" = "#2AE4F7", "3 - 9" = "#3ACCF4", "4 - 9" = "#4AB4F0","5 - 9" = "#5A9CED","6 - 9" = "#6A83E9","7 - 9" = "#7A6BE6","8 - 9" = "#8A53E2","9 - 9" = "#9A3BDF","10 - 9" = "#AA23DC",
  "1 - 8" = "#36FCF7", "2 - 8" = "#44E4F1", "3 - 8" = "#52CCEC", "4 - 8" = "#60B5E7","5 - 8" = "#6E9DE2","6 - 8" = "#7C85DC","7 - 8" = "#8A6ED7","8 - 8" = "#9856D2","9 - 8" = "#A63ECD","10 - 8" = "#B527C8",
  "1 - 7" = "#51FBF3", "2 - 7" = "#5DE3EC", "3 - 7" = "#69CCE5", "4 - 7" = "#75B5DE","5 - 7" = "#819ED7","6 - 7" = "#8E86D0","7 - 7" = "#9A6FC9","8 - 7" = "#A658C2","9 - 7" = "#B241BB","10 - 7" = "#BF2AB5",
  "1 - 6" = "#6CFAEF", "2 - 6" = "#76E3E6", "3 - 6" = "#80CCDD", "4 - 6" = "#8BB6D5","5 - 6" = "#959FCC","6 - 6" = "#A088C3","7 - 6" = "#AA72BB","8 - 6" = "#B55BB2","9 - 6" = "#BF44A9","10 - 6" = "#CA2EA1",
  "1 - 5" = "#88F9EB", "2 - 5" = "#90E2E0", "3 - 5" = "#98CCD6", "4 - 5" = "#A1B6CB","5 - 5" = "#A9A0C1","6 - 5" = "#B289B7","7 - 5" = "#BA73AD","8 - 5" = "#C35DA2","9 - 5" = "#CB4798","10 - 5" = "#D4318E",
  "1 - 4" = "#A3F8E7", "2 - 4" = "#A9E2DA", "3 - 4" = "#B0CCCE", "4 - 4" = "#B7B7C2","5 - 4" = "#BDA1B6","6 - 4" = "#C48BAA","7 - 4" = "#CB769E","8 - 4" = "#D16092","9 - 4" = "#D84A86","10 - 4" = "#DF357A",
  "1 - 3" = "#BEF7E3", "2 - 3" = "#C2E1D5", "3 - 3" = "#C7CCC7", "4 - 3" = "#CCB7B9","5 - 3" = "#D1A2AB","6 - 3" = "#D58C9E","7 - 3" = "#DA7790","8 - 3" = "#DF6282","9 - 3" = "#E44D74","10 - 3" = "#E93867",
  "1 - 2" = "#D9F6DF", "2 - 2" = "#DBE1CF", "3 - 2" = "#DFCCBF", "4 - 2" = "#E2B8B0","5 - 2" = "#E5A3A0","6 - 2" = "#E88E91","7 - 2" = "#EB7981","8 - 2" = "#EE6572","9 - 2" = "#F15062","10 - 2" = "#F43C53",
  "1 - 1" = "#F5F5DC", "2 - 1" = "#F6E0CA", "3 - 1" = "#F7CCB9", "4 - 1" = "#F8B8A8","5 - 1" = "#F9A496","6 - 1" = "#FA9085","7 - 1" = "#FB7C74","8 - 1" = "#FC6862","9 - 1" = "#FD5451","10 - 1" = "#FF4040",
  
) %>%
  gather("group", "fill")

#####(2)分配
cattleentr_df1 <- cattleentr_df %>%
  mutate(
    cattle_quantiles = cut(cattle_n, breaks = quantiles_cattle, include.lowest = TRUE, labels = FALSE),
    entr_quantiles = cut(entr_n, breaks = quantiles_entr, include.lowest = TRUE, labels = FALSE)
  ) %>%
  # 将分位数转换为数字，并创建一个新的组合列
  mutate(
    group = paste(as.numeric(cattle_quantiles), "-", as.numeric(entr_quantiles))
  ) %>%
  # 假设您已经有了一个颜色表bivariate_color_scale，其中包含group和对应的颜色值
  left_join(bivariate_color_scale, by = "group") # 将颜色值与主数据框合并



#####(3)绘制地图
library(rnaturalearth)
coast <- ne_coastline(scale = "small", returnclass = "sf") %>% vect()
crs <- '+proj=longlat +datum=WGS84'

map<-ggplot(cattleentr_df1) +
  geom_tile(aes(x = x, y = y,fill = fill)) +
  scale_fill_identity() +
  geom_spatvector(data=coast,fill=NA)+coord_sf(crs = crs,xlim=c(-160,165),ylim=c(-56,90))+
  #geom_sf(data = canton_geo, fill = NA, color = "white", size = 0.5) +    # 为行政区划绘制边界
  labs(x = NULL, y = NULL,) +
  #theme_minimal()+
  theme_bw()+
  theme(
    legend.position="none"
  )
print(map)


# bivariate_color_scale是颜色表，并且包含group和对应的颜色值
# 分离group列为cattle和entr两个独立的列，并转换为整数类型
bivariate_color_scale1 <- bivariate_color_scale %>%
  separate(group, into = c("cattle", "entr"), sep = " - ", convert = TRUE) %>%
  mutate(
    cattle = as.integer(cattle),
    entr = as.integer(entr)
  )

#绘制图例
legend <- ggplot(data = bivariate_color_scale1) +
  geom_tile(mapping = aes(x = cattle, y = entr, fill = fill)) +
  scale_fill_identity() +
  labs(x = "Higher cattle density ⟶️",
       y = "Higher WAE ⟶️") +
  theme_void() +  # 使用theme_void()来移除多余的元素
  theme(
    axis.title.x = element_text(size = 16, angle = 0),  # 旋转X轴标题
    axis.title.y = element_text(size = 16, angle = 90),   # 旋转Y轴标题
    plot.margin = margin(t = 10, r = 10, b = 30, l = 10)  # 调整图表边距
  ) +
  coord_fixed()  # 保持坐标比例固定
print(legend)

#hotsopt-----
##二元的：活动熵和牛
globalRaster <- rast(vals=1:259200,nrows=360, ncols=720,xmin=-180, xmax=180,ymin=-90, ymax=90,crs=crs)

popd2015<-rast("/root/autodl-tmp/全球人口/GWP_v4/popd2015_30.tif") %>% resample(globalRaster)
cattle2015<-rast("/root/autodl-tmp/全球家禽/Cattle/5_Ct_2015_Da.tif") %>% resample(globalRaster)
Entropy<-rast("/root/autodl-tmp/humPoulResult/data/AE_data/AE.tif")

#####国家矢量
globalSHP <- vect('/root/autodl-tmp/zyresult/Con_popentrpoul_sf_EU.shp')
globalSHP2 <- terra::aggregate(globalSHP,'name_ec')


###取阈值
quantiles_pop <- global(popd2015,quantile,probs=seq(0, 1, 0.01),na.rm=T)
quantiles_entr<- global(Entropy,quantile,probs=seq(0, 1, 0.01),na.rm=T)
quantiles_cattle<- global(cattle2015,quantile,probs=seq(0, 1, 0.01),na.rm=T)

quan_pop <- 50
quan_entr<- 4.19
quan_cattle <- quantiles_cattle[1,84]
#重分类
pop_new <- ifel(popd2015 > quan_pop,1,0)
entr_new <- ifel(Entropy > quan_entr,10,0)
cattle_new <- ifel(cattle2015 > quan_cattle,2,0)

hotentrpopcattle <- pop_new+entr_new+cattle_new
hotAND <- ifel(hotentrpopcattle==13,1,0)  ; plot(hotAND)
hotOR_pop <- ifel(hotentrpopcattle==11|hotentrpopcattle==13,1,0) ; plot(hotOR_pop)
hotOR_cattle <- ifel(hotentrpopcattle==12|hotentrpopcattle==13,1,0)  ; plot(hotOR_cattle)

#二元作图
coast <- ne_coastline(scale = "small", returnclass = "sf") %>% vect()
crs <- '+proj=longlat +datum=WGS84'

plot(hotOR_cattle)
hotOR_cattle_df<-as.data.frame(hotOR_cattle,xy=T) 
names(hotOR_cattle_df)<-c("x","y","hot")
ggplot(hotOR_cattle_df) +
  geom_tile(aes(x = x, y = y, fill = factor(hot))) +
  scale_fill_manual(values = c("0" = "grey", "1" = "#885c15")) +
  geom_spatvector(data=coast,fill=NA)+coord_sf(crs = crs,xlim=c(-160,165),ylim=c(-56,90))+
  labs(x = NULL, y = NULL,) +
  #theme_minimal()+
  theme_bw()+
  theme(
    legend.position="none"
  )








#c、d-----------
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


quan_pop <- 50
quan_entr<- 4.107
quan_poul <- quantiles_poul[1,84]
#重分类
pop_new <- ifel(popd2015>quan_pop,1,0)
entr_new <- ifel(Entropy>quan_entr,10,0)
poul_new <- ifel(poul2015>quan_poul,2,0)

hotentrpoppoul <- pop_new+entr_new+poul_new
hotAND <- ifel(hotentrpoppoul==13,1,0)  ; plot(hotAND)
hotOR_pop <- ifel(hotentrpoppoul==11|hotentrpoppoul==13,1,0) ; plot(hotOR_pop)
hotOR_poul <- ifel(hotentrpoppoul==12|hotentrpoppoul==13,1,0)  ; plot(hotOR_poul)



#作图
plot(hotOR_pop)
hotOR_pop_df<-as.data.frame(hotOR_pop,xy=T) 
names(hotOR_pop_df)<-c("x","y","hot")
ggplot(hotOR_pop_df) +
  geom_tile(aes(x = x, y = y, fill = factor(hot))) +
  scale_fill_manual(values = c("0" = "grey", "1" = "#f58a2c")) +
  geom_spatvector(data=coast,fill=NA)+coord_sf(crs = crs,xlim=c(-160,165),ylim=c(-56,90))+
  labs(x = NULL, y = NULL,) +
  #theme_minimal()+
  theme_bw()+
  theme(
    legend.position="none"
  )
#writeRaster(AEpop2,'/root/autodl-tmp/humPoulResult/data/Hot_data/Pop50entr.tif')

plot(hotOR_poul)
hotOR_poul_df<-as.data.frame(hotOR_poul,xy=T) 
names(hotOR_poul_df)<-c("x","y","hot")
ggplot(hotOR_poul_df) +
  geom_tile(aes(x = x, y = y, fill = factor(hot))) +
  scale_fill_manual(values = c("0" = "grey", "1" = "yellow")) +
  geom_spatvector(data=coast,fill=NA)+coord_sf(crs = crs,xlim=c(-160,165),ylim=c(-56,90))+
  labs(x = NULL, y = NULL,) +
  #theme_minimal()+
  theme_bw()+
  theme(
    legend.position="none"
  )
#writeRaster(AEpoul2,'/root/autodl-tmp/humPoulResult/data/Hot_data/Poul50entr.tif')



#e :  Continent (population)-------------
popd2015<-rast("/root/autodl-tmp/全球人口/GWP_v4/popd2015_30.tif") %>% resample(globalRaster)
Entropy<-rast(paste0(basePath,'AE_data/AE.tif'))%>% mask(globalCountry)
poul2015<-rast("/root/autodl-tmp/zyresult/Poultry_duckchic.tif")%>% resample(globalRaster)

Continent_vect<- vect ("/root/autodl-tmp/worldBorder/continentNew.shp") |> dplyr::filter(CONTINENT != "Antarctica")
ContinentRaster <- rasterize(Continent_vect,Entropy,field='CONTINENT')

Numpop_hotOR_pop<-hotOR_pop*popd2015; plot(Numpop_hotOR_pop)
Numpoul_hotOR_poul<-hotOR_poul*poul2015; plot(Numpoul_hotOR_poul)

AllhotOR_df <-c(Entropy,Numpop_hotOR_pop,Numpoul_hotOR_poul,ContinentRaster) %>% terra::as.data.frame() %>%na.omit()
head(AllhotOR_df)
names(AllhotOR_df)<-c("entr","Numpop_hotOR_pop","Numpoul_hotOR_poul","continent")

NumhotOR_result <- AllhotOR_df %>%
  group_by(continent) %>%
  summarise(across(everything(), ~ sum(.x, na.rm = TRUE)))
NumhotOR_result
#fwrite(NumhotOR_result,'/root/autodl-tmp/humPoulResult/data/Hot_data/hotOR_Num_result.csv')


#g : Country (poulation)---------------

popd2015<-rast("/root/autodl-tmp/全球人口/GWP_v4/popd2015_30.tif") %>% resample(globalRaster)
Entropy<-rast(paste0(basePath,'AE_data/AE.tif'))%>% mask(globalCountry)
poul2015<-rast("/root/autodl-tmp/zyresult/Poultry_duckchic.tif")%>% resample(globalRaster)

globalSHP <- vect('/root/autodl-tmp/zyresult/Con_popentrpoul_sf_EU.shp')
globalSHP2 <- terra::aggregate(globalSHP,'name_ec')
cRaster <- rasterize(globalSHP2,Entropy,field='name_ec')

Numpop_hotArea<-hotAND*popd2015; plot(Numpop_hotArea)
Numpoul_hotArea<-hotAND*poul2015; plot(Numpoul_hotArea)
NumAEpop<-hotOR_pop*popd2015  ; plot(NumAEpop)
NumAEpoul<-hotOR_poul*poul2015  ; plot(NumAEpoul)

All_df <-c(Entropy,Numpop_hotArea,Numpoul_hotArea,NumAEpop,NumAEpoul,cRaster) %>% terra::as.data.frame() %>%na.omit()
head(All_df)
names(All_df)<-c("entr","Numpop_hotArea","Numpoul_hotArea","NumAEpop","NumAEpoul","country")

Num_result <- All_df %>%
  group_by(country) %>%
  summarise(across(everything(), ~ sum(.x, na.rm = TRUE)))
Num_result
#fwrite(Num_result,'/root/autodl-tmp/humPoulResult/data/Hot_data/hot_Num_result.csv')



#——————————————————————————————————————————————————————————-----------------
#Fig S6 Mantel-test#########





#——————————————————————————————————————————————————————————-----------------

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

#Fig S7 China_Famliy#########

#计算中国的分科计算活动熵--------------
globalSHP <- vect('/root/autodl-tmp/zyresult/Con_popentrpoul_sf_EU.shp')
ChinaSHP <- subset(globalSHP,globalSHP$name_ec=="China")

allDf2 <- fread('/root/autodl-tmp/zyresult/allDf786_reclass.csv')
fam5China<-c("Anatidae","Scolopacidae","Laridae","Ardeidae","Charadriidae")
spName <- basename(speciesPixelNumPath2) %>% str_sub(.,1,-5)
for (fam in fam5China) {          #unique(allDf2$`Family`)
  famName <- allDf2[allDf2$`Family`==fam,]
  paths <- speciesPixelNumPath2[spName%in%famName$LatName]
  monthNum <- rast(paths) %>% sum(na.rm=T)
  calEntropy <- lapply(paths, function(x){
    r <- rast(x) %>% sum(.,na.rm=T)
    pi <- r/monthNum
    y <- -pi*log(pi)
    names(y) <- str_sub(basename(x),1,-5)
    return(y)
  })
  actEntropy <- rast(calEntropy) %>% sum(na.rm = T)
  China_actEntropy<-crop(actEntropy, ChinaSHP) %>% mask(., ChinaSHP)
  #writeRaster(China_actEntropy,paste0('/root/autodl-tmp/humPoulResult/data/AE_data/ChinaAE_',fam,'.tif'))
}

#记录中国的每个科的种数：
pathLengths <- data.table(family = character(), pathLength = numeric(), stringsAsFactors = FALSE)
for (fam in fam5China) {
  famName <- allDf2[allDf2$`Family` == fam,]
  paths <- speciesPixelNumPath2[spName %in% famName$LatName]
  
  # 初始化变量来存储中国范围内的路径长度
  pathLengthInChina <- 0
  for (path in paths) {
    rastData <- rast(path) %>% sum(na.rm=T)
    rastDataChina <- crop(rastData, ChinaSHP) %>% mask(., ChinaSHP)
    if (!all(is.na(values(rastDataChina)))) {   # 如果数据中有非NA值
      if (max(values(rastDataChina), na.rm = TRUE) > 0) {   # 检查是否有非零值（即存在鸟类的数据）
        pathLengthInChina <- pathLengthInChina + 1
      }
    }
  }
  # 将结果添加到数据表
  pathLengths <- rbind(pathLengths, data.table(family = fam, pathLength = pathLengthInChina))
}
print(pathLengths)

#计算中国的分科CH计算活动熵---------------
globalSHP <- vect('/root/autodl-tmp/zyresult/Con_popentrpoul_sf_EU.shp')
ChinaSHP <- subset(globalSHP,globalSHP$name_ec=="China")

allDf2 <- fread('/root/autodl-tmp/zyresult/allDf786_reclass.csv')
allDf2CH <- allDf2[allDf2$`Host`== "Confirmed Host",]
fam5China<-c("Anatidae","Scolopacidae","Laridae","Ardeidae","Charadriidae")
spName <- basename(speciesPixelNumPath2) %>% str_sub(.,1,-5)

for (fam in fam5China) {          #unique(allDf2CH$`Family`)
  famName <- allDf2CH[allDf2CH$`Family`==fam,]
  paths <- speciesPixelNumPath2[spName%in%famName$LatName]
  monthNum <- rast(paths) %>% sum(na.rm=T)
  calEntropy <- lapply(paths, function(x){
    r <- rast(x) %>% sum(.,na.rm=T)
    pi <- r/monthNum
    y <- -pi*log(pi)
    names(y) <- str_sub(basename(x),1,-5)
    return(y)
  })
  actEntropy <- rast(calEntropy) %>% sum(na.rm = T)
  China_actEntropy<-crop(actEntropy, ChinaSHP) %>% mask(., ChinaSHP)
  #writeRaster(China_actEntropy,paste0('/root/autodl-tmp/humPoulResult/data/AE_data/ChinaAECH_',fam,'.tif'))
}

#记录中国的每个科CH的种数：
pathLengths <- data.table(family = character(), pathLength = numeric(), stringsAsFactors = FALSE)
for (fam in fam5China) {
  famName <- allDf2CH[allDf2CH$`Family` == fam,]
  paths <- speciesPixelNumPath2[spName %in% famName$LatName]
  
  # 初始化变量来存储中国范围内的路径长度
  pathLengthInChina <- 0
  for (path in paths) {
     rastData <- rast(path) %>% sum(na.rm=T)
     rastDataChina <- crop(rastData, ChinaSHP) %>% mask(., ChinaSHP)
     if (!all(is.na(values(rastDataChina)))) {   # 如果数据中有非NA值
       if (max(values(rastDataChina), na.rm = TRUE) > 0) {   # 检查是否有非零值（即存在鸟类的数据）
         pathLengthInChina <- pathLengthInChina + 1
       }
     }
  }
  # 将结果添加到数据表
  pathLengths <- rbind(pathLengths, data.table(family = fam, pathLength = pathLengthInChina))
}
print(pathLengths)


#（1）总 不同科 的水鸟熵与中国人口家禽做回归-----------
globalSHP <- vect('/root/autodl-tmp/zyresult/Con_popentrpoul_sf_EU.shp')
ChinaSHP <- subset(globalSHP,globalSHP$name_ec=="China")

AEPath <- list.files('/root/autodl-tmp/humPoulResult/data/AE_data',pattern = 'ChinaAE_',full.names = T)
ChinaAE <- rast(AEPath)
China_Pop<-rast("/root/autodl-tmp/全球人口/GWP_v4/popd2015_30.tif") %>% mask(., ChinaSHP) %>% resample(.,ChinaAE, method="bilinear")
China_Poul<-rast("/root/autodl-tmp/zyresult/Poultry_duckchic.tif") %>% mask(., ChinaSHP) %>% resample(.,ChinaAE, method="bilinear")

ChinaAEpoppul_df <- c(ChinaAE,China_Pop,China_Poul) %>% as.data.frame(xy=T)
names(ChinaAEpoppul_df)<-c('x','y','Anatidae','Ardeidae','Charadriidae','Laridae','Scolopacidae','pop','poul')  #,'Other'  ,'Charadriidae'
head(ChinaAEpoppul_df)

#出图
ChinaAEpoppul_df_long <- reshape2::melt(ChinaAEpoppul_df, id.vars = c("x", "y", "pop", "poul"), 
                                     measure.vars = c('Anatidae','Ardeidae','Charadriidae','Laridae','Scolopacidae'))  #,'Other'  #,'Charadriidae'
ChinaAEpoppul_df_long2<-subset(ChinaAEpoppul_df_long, value>0 & pop>0.0001 & poul>0.0001)


#####【水鸟】提琴图---------

x2<-factor(ChinaAEpoppul_df_long2$variable,levels = c('Charadriidae','Ardeidae','Laridae','Scolopacidae','Anatidae')) #自定义展示顺序

ggplot(ChinaAEpoppul_df_long2, aes(x = x2, y = value,fill=x2)) + #,
  scale_fill_manual(values = c("#FF7F00", "#984EA3", "#4DAF4A", "#E41A1C", "#377EB8")) +
  #scale_fill_manual(values = c("#377EB8", "#E41A1C", "#4DAF4A", "#984EA3", "#FF7F00")) +
  geom_hline(yintercept = mean(ChinaAEpoppul_df_long2$value), size = 1, color = "grey", linetype = "dashed")+
  #scale_fill_manual(values = c("#E0312D","#F6BA26","#8204FB","#3ca1f9","#45F3DD","#49FB18"))+
  #geom_jitter(width = 0.1) + #
  #geom_jitter(alpha=0.5,position=position_jitterdodge(jitter.width = 0.5, jitter.height = 0, dodge.width = 0.8))+
  geom_violin(alpha=0.5,position = position_dodge(width = 0.01), scale = 'width') +
  geom_boxplot(alpha=0.8,width=0.45,position=position_dodge(width=0.1),size=0.75,outlier.colour = NA)+
  #scale_x_continuous(position = "top")+
  scale_y_continuous(position = "right")+
  labs(y = "", x = " ",title = "Entropy")+
  theme_bw()+
  theme(panel.grid = element_blank(), 
        panel.background = element_rect(fill = 'transparent', color = "grey"), 
        legend.position = 'none',
        axis.text = element_text(size = 26),
        text = element_text(size = 26),
  ) +
  coord_flip()#转换横纵方向

##### 【人口】散点图并添加线性回归拟合线-------
p1<-ggplot(data = ChinaAEpoppul_df_long2, aes(x = value, y = log(pop), color = variable)) +
  geom_point(size=1.5) +  #,shape=21
  scale_color_manual(values = c("#377EB8", "#E41A1C", "#4DAF4A", "#984EA3", "#FF7F00")) +#,"#FA9085"
  geom_smooth(method = "lm",linewidth=0.5, se = T, lty=5,aes(color = variable)) +  #,aes(color = variable)   , color="black"
  #geom_smooth(method = "lm", se = T, color = "black", linetype = "dashed") +  # 全部点的趋势线，黑色
  #geom_smooth(method = "lm",linewidth=0.5, se = T, lty=5,color = "black") +  #分界
  # stat_poly_eq(aes(label = after_stat(paste(eq.label, rr.label, sep = "~~~"))),
  #              formula = y ~ x,  label.x = Inf, # 右边界
  #              #label.y=-Inf,
  #              coef.digits = 3, rr.digits = 2, parse = TRUE, size = 3) +
  #facet_wrap(~ variable) +  #如果要分界
  labs(x = "WAE", y = "Human density (Log)",title = "All Waterbirds") +
  #ylim(-7,10)+
  theme_bw() +
  theme(legend.position = "none",
        legend.title = element_blank(),
        plot.title = element_text(size = 24, hjust = 0.5),
        #axis.title.x = element_blank(),
        #axis.title.y = element_blank(),
        text = element_text(size = 18))
p1


#计算参数
head(ChinaAEpoppul_df_long2)
# 初始化一个空的数据框来存储结果
popresults_df <- data.frame(  Variable = character(),  b = numeric(),  SE = numeric(), b_SE = character(),  z = numeric(),  P = numeric(),  R2 = numeric(),  stringsAsFactors = FALSE)
variables <- c('Charadriidae','Ardeidae','Laridae','Scolopacidae','Anatidae')

# 循环每个变量，进行回归分析，并将结果添加到数据框
for (var in variables) {
  model <- lm(log(pop) ~ value, data = ChinaAEpoppul_df_long2[ChinaAEpoppul_df_long2$variable == var, ])
  model_summary <- summary(model)
  
  # 提取所需的参数
  b <- model_summary$coefficients[2, 1]
  SE <- model_summary$coefficients[2, 2]
  b_SE <- sprintf("%.2f ± %.2f", b, SE)
  z <- b / SE
  P <- model_summary$coefficients[2, 4]
  R2 <- model_summary$r.squared
  
  # 将结果添加到数据框
  popresults_df <- rbind(popresults_df, data.frame(Variable = var, b_pop = b, SE_pop = SE, b_SE_pop = b_SE, z_pop = z, P_pop = P, R2_pop = R2))
}
head(popresults_df)

##### 【家禽】散点图并添加线性回归拟合线---------
p2<-ggplot(data = ChinaAEpoppul_df_long2, aes(x = value, y = log(poul), color = variable)) +
  geom_point(size=1.5) +  #,shape=21
  scale_color_manual(values = c("#377EB8", "#E41A1C", "#4DAF4A", "#984EA3", "#FF7F00")) +  #,"#FA9085"
  geom_smooth(method = "lm",linewidth=0.5, se = T, lty=5,aes(color = variable)) +  #,aes(color = variable)   , color="black"
  #geom_smooth(method = "lm", se = T, color = "black", linetype = "dashed") +  # 全部点的趋势线，黑色
  #geom_smooth(method = "lm",linewidth=0.5, se = T, lty=5,color = "black") +  #分界
  # stat_poly_eq(aes(label = after_stat(paste(eq.label, rr.label, sep = "~~~"))),
  #              formula = y ~ x,  label.x = Inf, # 右边界
  #              #label.y=-Inf,
  #              coef.digits = 3, rr.digits = 2, parse = TRUE, size = 3) +
  #facet_wrap(~ variable) +  #如果要分界
  labs(x = "WAE", y = "Poultry density(Log)",title = "All Waterbirds") +
  #ylim(-7,10)+
  theme_bw() +
  theme(legend.position = "none",
        legend.title = element_blank(),
        plot.title = element_text(size = 24, hjust = 0.5),
        #axis.title.x = element_blank(),
        #axis.title.y = element_blank(),
        text = element_text(size = 18))
p2

#计算参数
head(ChinaAEpoppul_df_long2)
# 初始化一个空的数据框来存储结果
poulresults_df <- data.frame(  Variable = character(),  b = numeric(),  SE = numeric(), b_SE = character(),  z = numeric(),  P = numeric(),  R2 = numeric(),  stringsAsFactors = FALSE)
variables <- c('Charadriidae','Ardeidae','Laridae','Scolopacidae','Anatidae')

# 循环每个变量，进行回归分析，并将结果添加到数据框
for (var in variables) {
  model <- lm(log(poul) ~ value, data = ChinaAEpoppul_df_long2[ChinaAEpoppul_df_long2$variable == var, ])
  model_summary <- summary(model)
  # 提取所需的参数
  b <- model_summary$coefficients[2, 1]
  SE <- model_summary$coefficients[2, 2]
  b_SE <- sprintf("%.2f ± %.2f", b, SE)
  z <- b / SE
  P <- model_summary$coefficients[2, 4]
  R2 <- model_summary$r.squared
  poulresults_df <- rbind(poulresults_df, data.frame(Variable = var, b_poul = b, SE_poul = SE, b_SE_poul = b_SE, z_poul = z, P_poul = P, R2_poul = R2))
}
head(poulresults_df)

#合并人口和家禽的线性拟合结果
ChinaFam_poppoul_r<-left_join(popresults_df, poulresults_df, by="Variable")
#fwrite(ChinaFam_poppoul_r,"/root/autodl-tmp/humPoulResult/data/AE_data/ChinaFam_poppoul_r.csv")



#（2）宿主CH 不同科 的水鸟熵与中国人口家禽做回归-----------
envariables <- ls()
envariables<- envariables[!(envariables %in% c("p1", "p2"))]
rm(list = envariables)

globalSHP <- vect('/root/autodl-tmp/zyresult/Con_popentrpoul_sf_EU.shp')
ChinaSHP <- subset(globalSHP,globalSHP$name_ec=="China")

AECHPath <- list.files('/root/autodl-tmp/humPoulResult/data/AE_data',pattern = 'ChinaAECH_',full.names = T)
ChinaAE <- rast(AECHPath)
China_Pop<-rast("/root/autodl-tmp/全球人口/GWP_v4/popd2015_30.tif") %>% mask(., ChinaSHP) %>% resample(.,ChinaAE, method="bilinear")
China_Poul<-rast("/root/autodl-tmp/zyresult/Poultry_duckchic.tif") %>% mask(., ChinaSHP) %>% resample(.,ChinaAE, method="bilinear")

ChinaAEpoppul_df <- c(ChinaAE,China_Pop,China_Poul) %>% as.data.frame(xy=T)
names(ChinaAEpoppul_df)<-c('x','y','Anatidae','Ardeidae','Charadriidae','Laridae','Scolopacidae','pop','poul')  #,'Other'  ,'Charadriidae'
head(ChinaAEpoppul_df)

#出图
ChinaAEpoppul_df_long <- reshape2::melt(ChinaAEpoppul_df, id.vars = c("x", "y", "pop", "poul"), 
                                        measure.vars = c('Anatidae','Ardeidae','Charadriidae','Laridae','Scolopacidae'))  #,'Other'  #,'Charadriidae'
ChinaAEpoppul_df_long2<-subset(ChinaAEpoppul_df_long, value>0 & pop>0.0001 & poul>0.0001)


#####【水鸟】提琴图---------
# 
# x2<-factor(ChinaAEpoppul_df_long2$variable,levels = c('Charadriidae','Ardeidae','Laridae','Scolopacidae','Anatidae')) #自定义展示顺序
# 
# ggplot(ChinaAEpoppul_df_long2, aes(x = x2, y = value,fill=x2)) + #,
#   scale_fill_manual(values = c("#FF7F00", "#984EA3", "#4DAF4A", "#E41A1C", "#377EB8")) +
#   #scale_fill_manual(values = c("#377EB8", "#E41A1C", "#4DAF4A", "#984EA3", "#FF7F00")) +
#   geom_hline(yintercept = mean(ChinaAEpoppul_df_long2$value), size = 1, color = "grey", linetype = "dashed")+
#   #scale_fill_manual(values = c("#E0312D","#F6BA26","#8204FB","#3ca1f9","#45F3DD","#49FB18"))+
#   #geom_jitter(width = 0.1) + #
#   #geom_jitter(alpha=0.5,position=position_jitterdodge(jitter.width = 0.5, jitter.height = 0, dodge.width = 0.8))+
#   geom_violin(alpha=0.5,position = position_dodge(width = 0.01), scale = 'width') +
#   geom_boxplot(alpha=0.8,width=0.45,position=position_dodge(width=0.1),size=0.75,outlier.colour = NA)+
#   #scale_x_continuous(position = "top")+
#   scale_y_continuous(position = "right")+
#   labs(y = "", x = " ",title = "Entropy")+
#   theme_bw()+
#   theme(panel.grid = element_blank(), 
#         panel.background = element_rect(fill = 'transparent', color = "grey"), 
#         legend.position = 'none',
#         axis.text = element_text(size = 26),
#         text = element_text(size = 26),
#   ) +
#   coord_flip()#转换横纵方向

##### 【人口】散点图并添加线性回归拟合线-------
p3<-ggplot(data = ChinaAEpoppul_df_long2, aes(x = value, y = log(pop), color = variable)) +
  geom_point(size= 1.5) +  #,shape=21
  scale_color_manual(values = c("#377EB8", "#E41A1C", "#4DAF4A", "#984EA3", "#FF7F00")) +#,"#FA9085"
  geom_smooth(method = "lm",linewidth=0.5, se = T, lty=5,aes(color = variable)) +  #,aes(color = variable)   , color="black"
  #geom_smooth(method = "lm", se = T, color = "black", linetype = "dashed") +  # 全部点的趋势线，黑色
  #geom_smooth(method = "lm",linewidth=0.5, se = T, lty=5,color = "black") +  #分界
  # stat_poly_eq(aes(label = after_stat(paste(eq.label, rr.label, sep = "~~~"))),
  #              formula = y ~ x,  label.x = Inf, # 右边界
  #              #label.y=-Inf,
  #              coef.digits = 3, rr.digits = 2, parse = TRUE, size = 3) +
  #facet_wrap(~ variable) +  #如果要分界
  labs(x = "WAE", y = "Human density (Log)",title = "Confirmed Host") +
  #ylim(-7,10)+
  theme_bw() +
  theme(legend.position = "none",
        legend.title = element_blank(),
        plot.title = element_text(size = 24, hjust = 0.5),
        #axis.title.x = element_blank(),
        #axis.title.y = element_blank(),
        text = element_text(size = 18))
p3


#计算参数
head(ChinaAEpoppul_df_long2)
# 初始化一个空的数据框来存储结果
popresults_df <- data.frame(  Variable = character(),  b = numeric(),  SE = numeric(), b_SE = character(),  z = numeric(),  P = numeric(),  R2 = numeric(),  stringsAsFactors = FALSE)
variables <- c('Charadriidae','Ardeidae','Laridae','Scolopacidae','Anatidae')

# 检查每个类别的数据点数量发现Charadriidae为0
table(ChinaAEpoppul_df_long2$variable)
variables <- c('Ardeidae','Laridae','Scolopacidae','Anatidae')

# 循环每个变量，进行回归分析，并将结果添加到数据框
for (var in variables) {
  model <- lm(log(pop) ~ value, data = ChinaAEpoppul_df_long2[ChinaAEpoppul_df_long2$variable == var, ])
  model_summary <- summary(model)
  
  # 提取所需的参数
  b <- model_summary$coefficients[2, 1]
  SE <- model_summary$coefficients[2, 2]
  b_SE <- sprintf("%.2f ± %.2f", b, SE)
  z <- b / SE
  P <- model_summary$coefficients[2, 4]
  R2 <- model_summary$r.squared
  
  # 将结果添加到数据框
  popresults_df <- rbind(popresults_df, data.frame(Variable = var, b_pop = b, SE_pop = SE, b_SE_pop = b_SE, z_pop = z, P_pop = P, R2_pop = R2))
}
head(popresults_df)

##### 【家禽】散点图并添加线性回归拟合线---------
p4<-ggplot(data = ChinaAEpoppul_df_long2, aes(x = value, y = log(poul), color = variable)) +
  geom_point(size=1.5) +  #,shape=21
  scale_color_manual(values = c("#377EB8", "#E41A1C", "#4DAF4A", "#984EA3", "#FF7F00")) +  #,"#FA9085"
  geom_smooth(method = "lm",linewidth=0.5, se = T, lty=5,aes(color = variable)) +  #,aes(color = variable)   , color="black"
  #geom_smooth(method = "lm", se = T, color = "black", linetype = "dashed") +  # 全部点的趋势线，黑色
  #geom_smooth(method = "lm",linewidth=0.5, se = T, lty=5,color = "black") +  #分界
  # stat_poly_eq(aes(label = after_stat(paste(eq.label, rr.label, sep = "~~~"))),
  #              formula = y ~ x,  label.x = Inf, # 右边界
  #              #label.y=-Inf,
  #              coef.digits = 3, rr.digits = 2, parse = TRUE, size = 3) +
  #facet_wrap(~ variable) +  #如果要分界
  labs(x = "WAE", y = "Poultry density(Log)",title = "Confirmed Host") +
  #ylim(-7,10)+
  theme_bw() +
  theme(legend.position = "none",
        legend.title = element_blank(),
        plot.title = element_text(size = 24, hjust = 0.5),
        #axis.title.x = element_blank(),
        #axis.title.y = element_blank(),
        
        text = element_text(size = 18))
p4

#计算参数
head(ChinaAEpoppul_df_long2)
# 初始化一个空的数据框来存储结果
poulresults_df <- data.frame(  Variable = character(),  b = numeric(),  SE = numeric(), b_SE = character(),  z = numeric(),  P = numeric(),  R2 = numeric(),  stringsAsFactors = FALSE)
variables <- c('Charadriidae','Ardeidae','Laridae','Scolopacidae','Anatidae')
# 检查每个类别的数据点数量发现Charadriidae为0
table(ChinaAEpoppul_df_long2$variable)
variables <- c('Ardeidae','Laridae','Scolopacidae','Anatidae')

# 循环每个变量，进行回归分析，并将结果添加到数据框
for (var in variables) {
  model <- lm(log(poul) ~ value, data = ChinaAEpoppul_df_long2[ChinaAEpoppul_df_long2$variable == var, ])
  model_summary <- summary(model)
  # 提取所需的参数
  b <- model_summary$coefficients[2, 1]
  SE <- model_summary$coefficients[2, 2]
  b_SE <- sprintf("%.2f ± %.2f", b, SE)
  z <- b / SE
  P <- model_summary$coefficients[2, 4]
  R2 <- model_summary$r.squared
  poulresults_df <- rbind(poulresults_df, data.frame(Variable = var, b_poul = b, SE_poul = SE, b_SE_poul = b_SE, z_poul = z, P_poul = P, R2_poul = R2))
}
head(poulresults_df)


#合并人口和家禽的线性拟合结果
ChinaFamCH_poppoul_r<-left_join(popresults_df, poulresults_df, by="Variable")
#fwrite(ChinaFamCH_poppoul_r,"/root/autodl-tmp/humPoulResult/data/AE_data/ChinaFamCH_poppoul_r.csv")

#全部图：
library(patchwork)
p3+p4+p1+p2 + plot_layout(ncol = 4)

