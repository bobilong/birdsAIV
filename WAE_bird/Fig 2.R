
#——————————————————————————————————————————————————————————-----------------

#2.Fig2. AE validation#########################
basePath<-"/root/autodl-tmp/root/autodl-tmp/WAEdata_new_y/result/"
crs <- '+proj=longlat +datum=WGS84'

world.map <- rnaturalearth::ne_countries(returnclass = "sf") |> dplyr::filter(continent != "Antarctica")
globalCountry <- vect(world.map) 
globalRaster <- rast(vals=1:259200,nrows=360, ncols=720,xmin=-180, xmax=180,ymin=-90, ymax=90,crs=crs)
coast <- rnaturalearth::ne_coastline(scale = "small", returnclass = "sf")

###a##########
overEntropy <- rast(paste0(basePath,'WAE.tif'))%>% mask(globalCountry)
#overEntropy<-AE
global(overEntropy,quantile,probs=seq(0, 1, 0.05),na.rm=T)


`%notin%` <- Negate(`%in%`)
library(lubridate)
#outBreak1 <- fread('/root/autodl-tmp/YANZHENG/point/allData.csv') %>% subset(Diagnosis.status=='Confirmed'&Animal.type%in%c('Domestic','Wild'))
outBreak1 <- fread('/root/autodl-tmp/root/autodl-tmp/WAEdata_new_y/result/otherdata/outBreakData.csv') %>% subset(Diagnosis.status=='Confirmed'&Animal.type%in%c('Domestic','Wild'))
outBreak1$label <- str_extract(outBreak1$Serotype,'HPAI|LPAI')
outBreak1$h_label <- str_extract(outBreak1$Serotype,'H[0-9]N[0-9]|H[0-9]')
outBreak1 <- subset(outBreak1,outBreak1$h_label%notin%c('H9N2','H5N6'))   #这一句是用来运行全部的
#outBreak1 <- subset(outBreak1,outBreak1$h_label=='H5N1'&outBreak1$label=='HPAI')  #这一句是用来仅运行H5N1-HPAI的

outBreak1$Longitude <- as.numeric(outBreak1$Longitude)

addGeom <- cellFromXY(globalRaster,outBreak1[,c('Longitude','Latitude')]) %>% 
  cbind(id=.,outBreak1)
thinData <- unique(data.table(addGeom),by='id') %>% dplyr::select(.,-id)
outBreak2 <- vect(thinData,geom=c('Longitude','Latitude'),crs=crs)

thinData_H5N1 <- thinData[thinData$h_label == 'H5N1', ]
outBreak2_H5N1 <- vect(thinData_H5N1,geom=c('Longitude','Latitude'),crs=crs)

type_counts <- outBreak2 %>% as.data.frame() %>% 
  group_by(h_label) %>% 
  summarize(count = n())
type_counts

plotEntropy <-  ifel(overEntropy>=4.11,1,0)


coast <- rnaturalearth::ne_coastline(scale = "small", returnclass = "sf")
crs <- '+proj=longlat +datum=WGS84'
ggplot() +
  geom_spatraster(data = plotEntropy) +
  scale_fill_gradient(low = "lightgrey",high = "yellow" ,space = "Lab",n.break=2,
                      labels=c('Low risk','High risk'),na.value='white')+
  geom_spatvector(data=coast,fill=NA)+
  #coord_sf(crs = crs,xlim=c(-160,180),ylim=c(-56,90))+
  geom_spatvector(data=outBreak2,size=2, shape = 4,color='#ff4f4c',fill='#f26c6a', alpha = 0.8,stroke = 0.3)+   #shape = 21是圆圈，shape = 4是×
  geom_spatvector(data = outBreak2_H5N1, size = 2, shape = 1, color = '#370d51', fill = '#e040fb', alpha = 0.8, stroke = 0.3) +
  #geom_spatvector(data = outBreak2_H5N1, size = 2, shape = 4, color = '#370d51', fill = '#e040fb', alpha = 0.8, stroke = 0.3) + # shape = 4 表示全部点 
  theme_bw()+
  #scale_color_manual(values='#fc4e4e', labels='Outbreak of navian influenza') +
  guides(
    fill=guide_legend()  )+
  labs(title = NULL)+
  theme(
    text = element_text(size = 15),
    #axis.text = element_text(size = 15),
    plot.title = element_text(hjust=0.5),
    # axis.line = element_line(color = 'black'),
    # panel.background = element_rect(fill = "white"),#
    legend.title = element_blank(),
    legend.text = element_text(size = 15),
    legend.position = "none",
    #legend.position = c(0.15, 0.15),#
    legend.direction='horizontal',#
    legend.key.width = unit(0.4,'cm'), #
    legend.key.height = unit(3,'cm')
  )



popd2015<-rast("/root/autodl-tmp/root/autodl-tmp/WAEdata_new_y/result/otherdata/pop.tif") %>% resample(globalRaster)
poul2015<-rast("/root/autodl-tmp/root/autodl-tmp/WAEdata_new_y/result/otherdata/poultry.tif")%>% resample(globalRaster)
cattle <- rast("/root/autodl-tmp/root/autodl-tmp/WAEdata_new_y/result/otherdata/cattle.tif")%>% resample(globalRaster)
plot(plotEntropy);global(plotEntropy, sum, na.rm=T)

plotEntropy_df<-c(plotEntropy,popd2015,cattle,poul2015) %>% terra::as.data.frame() %>%na.omit()
head(plotEntropy_df)
names(plotEntropy_df)<-c("plotEntropy","Numpop","Numcattle","Numpoul")

plotEntropy_result <- plotEntropy_df %>%
  group_by(plotEntropy) %>%
  summarise(across(everything(), ~ sum(.x, na.rm = TRUE)))
plotEntropy_result
#fwrite(plotEntropy_result,'/root/autodl-tmp/root/autodl-tmp/WAEdata_new_y/result/hotOR_Num_result.csv')

# 
ones <- sum(values(plotEntropy) == 1, na.rm = TRUE); ones
# 
zeros <- sum(values(plotEntropy) == 0, na.rm = TRUE); zeros


###a little##########
###ROC##############
library(pROC)
`%notin%` <- Negate(`%in%`)
crs <- '+proj=longlat +datum=WGS84'
globalRaster <- rast(vals=1:259200,nrows=360, ncols=720,xmin=-180, xmax=180,ymin=-90, ymax=90,crs=crs)
globalSHP <- vect('/root/autodl-tmp/allData/worldBorder/continentNew.shp')
# calculate  threshold
outBreakData <- fread('/root/autodl-tmp/root/autodl-tmp/WAEdata_new_y/result/otherdata/outBreakData.csv')


outBreakData$label <- str_extract(outBreakData$Serotype,'HPAI|LPAI')
outBreakData$h_label <- str_extract(outBreakData$Serotype,'H[0-9]N[0-9]|H[0-9]')

outBreakData2 <- na.omit(outBreakData[,c('Country','Longitude','Latitude')])
addGeom <- cellFromXY(globalRaster,outBreakData2[,c('Longitude','Latitude')]) %>% 
  cbind(id=.,outBreakData2)
thinData <- unique(data.table(addGeom),by='id') %>% dplyr::select(.,-id) #
breakVect <- vect(thinData,geom=c('Longitude','Latitude'),crs=crs)
breakRast <- rasterize(breakVect,globalRaster)

TP_points <- vect(outBreakData2,geom=c('Longitude','Latitude'),crs=crs)
AE <- rast('/root/autodl-tmp/root/autodl-tmp/WAEdata_new_y/result/WAE.tif')   #


#writeRaster(AE, '/root/autodl-tmp/allData/result/re_result/AE.tif',overwrite=T)
sample_TP <- terra::extract(AE,TP_points,min,bind=T) %>% as.data.frame() %>% na.omit()
#sample_TP <- subset(sample_TP,sample_TP$sum>3.6)
sample_TP$label <- 'good'


trueNgRast <- mask(AE,breakRast,inverse=T)

#set.seed(1234)
set.seed(1)
sample_TN<- terra::spatSample(trueNgRast,dim(sample_TP)[1],method='random',exhaustive=T,na.rm=T)
sample_TN$label <- 'poor'
sampleData <- rbind(sample_TN,sample_TP[,c('sum','label')])
dfroc1<- roc(sampleData$label, sampleData$sum)
dfroc1$auc

#litter
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


threshold <- dfroc1$auc
basePath<-'/root/autodl-tmp/root/autodl-tmp/WAEdata_new_y/result/'
WAE<- rast(paste0(basePath,'WAE.tif'))%>% mask(globalCountry)
plotWAE <- ifel(WAE>=4.11,1,0)   
plot(plotWAE)

ggplot() +
  geom_spatraster(data = plotWAE) +
  scale_fill_gradient(low = "lightgrey",high = "yellow" ,space = "Lab",n.break=2,
                      labels=c('Low risk','High risk'),na.value='white')+
  geom_spatvector(data=coast,fill=NA)




###differnt type of AIV#####
crs <- '+proj=longlat +datum=WGS84'
world.map <- rnaturalearth::ne_countries(returnclass = "sf") |> dplyr::filter(continent != "Antarctica")
globalCountry <- vect(world.map) 
globalRaster <- rast(vals=1:259200,nrows=360, ncols=720,xmin=-180, xmax=180,ymin=-90, ymax=90,crs=crs)
coast <- rnaturalearth::ne_coastline(scale = "small", returnclass = "sf")
`%notin%` <- Negate(`%in%`)

AE <- rast('/root/autodl-tmp/root/autodl-tmp/WAEdata_new_y/result/WAE.tif')
outBreakData <- fread('/root/autodl-tmp/root/autodl-tmp/WAEdata_new_y/result/otherdata/outBreakData.csv') #%>% subset(Diagnosis.status=='Confirmed'&Animal.type%in%c('Domestic','Wild'))

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

H5N1_df <- data.table(
  label='H5N1',
  cases=nrow(h5n1_df),
  recall=nrow(subset(h5n1_value,h5n1_value$sum>4.00))/nrow(h5n1_value)
)

Hx_df <- data.table(
  label='all AIV',
  cases=nrow(hx_df),
  recall=nrow(subset(hx_value,hx_value$sum>4.11))/nrow(hx_value)
)

result_st <- rbind(H5N1_df,Hx_df,fill=T)
result_st


h_df <- subset(outBreakData,outBreakData$label_hl%in%c('HPAI'))
l_df <- subset(outBreakData,outBreakData$label_hl%in%c('LPAI'))

h_vect <- vect(h_df,geom=c('Longitude','Latitude'),crs=crs)
l_vect <- vect(l_df,geom=c('Longitude','Latitude'),crs=crs)

h_value <- terra::extract(AE,h_vect,mean,na.rm=T) %>% 
  as.data.frame() %>% na.omit()
l_value <- terra::extract(AE,l_vect,mean,na.rm=T) %>% 
  as.data.frame() %>% na.omit()

HPAI_df <- data.table(
  label='HPAI',
  cases=nrow(h_df),
  recall=nrow(subset(h_value,h_value$sum>4.00))/nrow(h_value)
)

LPAI_df <- data.table(
  label='LPAI',
  cases=nrow(l_df),
  recall=nrow(subset(l_value,l_value$sum>4.15))/nrow(l_value)
)

result_st4 <- rbind(result_st,HPAI_df,fill=T) %>% rbind(.,LPAI_df,fill=T)
result_st4

#fwrite(result_st4, '/root/autodl-tmp/root/autodl-tmp/WAEdata_new_y/result/Serotype_edit.csv')





###d: serotype----------
result_st2 <- fread('/root/autodl-tmp/root/autodl-tmp/WAEdata_new_y/result/Serotype_edit.csv')
result_st2


ggplot(data = result_st2, aes(x = label, y = recall, fill = label)) +
  geom_col(position = 'dodge2') +
  theme_bw() +  # 
  labs( y = "Sensitivity") +  
  scale_fill_brewer(palette = "Pastel1") +
  ylim(c(0, 0.9))+
  theme(
    axis.title.x = element_blank(),
    text = element_text(size=18),
    plot.title = element_text(hjust = 0.5),  #
    legend.position = "none"  
  )
result_st3 <- mutate(result_st2, df_label = cut(cases, breaks = c(1500, 2000, 10000, 40000, 60000),
                                                labels = c("1500 - 2000", "5000 - 6000", "35000 - 40000", "> 50000")))
size_values <- c("1500 - 2000" = 4, "5000 - 6000" = 8, "35000 - 40000" = 14, "> 50000" =20)
result_st3$Group<-c("Subtype","Subtype","Pathogenicity","Pathogenicity")
#result_st3 <- result_st3 %>%  mutate(recall = percent(recall))


# 
result_st3$Group<-factor(result_st3$Group,levels = c("Subtype","Pathogenicity"))
result_st3$label<-factor(result_st3$label,levels = c("H5N1","all AIV","HPAI","LPAI"))

ggplot(data = result_st3, aes(x = label, y = recall)) +
  geom_point(aes(size = df_label, color = Group), alpha = 0.5) +  #
  geom_hline(aes(yintercept = 0.78), color = "grey", linetype = "dashed", alpha = 1) +
  scale_size_manual(values = size_values) +  # 
  scale_color_manual(values = c("orange", "red")) +  # 
  ylim(c(0.73, 0.82)) +
  labs(size = "Number of Cases", y = "Sensitivity") +  # 
  theme_bw() +
  theme(
    axis.text = element_text(size = 12),
    axis.title.y = element_text(size = 16),
    axis.title.x = element_blank()
  )






#b、c###########
#accuracy
result_country <- fread('/root/autodl-tmp/root/autodl-tmp/WAEdata_new_y/result/result_country.csv')
globalSHP <- vect('/root/autodl-tmp/root/autodl-tmp/WAEdata_new_y/result/otherdata/ne_10m_admin_0_countries.shp')
countryData <- as.data.frame(globalSHP) %>% .[,c('NAME_LONG','POP_EST','POP_RANK','ECONOMY','INCOME_GRP','CONTINENT','GDP_MD')]
result_df4 <- left_join(result_country,countryData,by=c('country'='NAME_LONG'))
result_df4$per <- (result_df4$df/result_df4$alldf)


#GDP
# POP_EST
pop_quantiles <- quantile(result_df4$POP_EST, probs=seq(0, 1, by=0.1), na.rm = TRUE)

# POP_EST to 1-10
result_df4$POP_EST_size <- findInterval(result_df4$POP_EST, vec = pop_quantiles)
result_df4$POP_EST <- as.numeric(result_df4$POP_EST)
result_df4 <- mutate(result_df4, POP_EST_label = cut(POP_EST, breaks = c(0, 10000000, 20000000, 50000000, 100000000, 500000000, 1000000000, 1397800000),
                                                     labels = c("<1000", "1000 - 2000", "2000-5000", "5000 - 10000", "10000 - 50000", "50000 - 100000", ">100000")))

# 
result_df4 <- result_df4[order(result_df4$CONTINENT, -result_df4$POP_EST), ]

#
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
# # 
# linear_model <- lm(accuracy ~ log(POP_EST), data = result_df4)
# summary(linear_model)
# 
#size_labels <- c("<1000", "1000 - 2000", "2000-5000", "5000 - 10000", "10000 - 50000", "50000 - 100000", ">100000")
size_values <- c("<1000" = 3, "1000 - 2000" = 4, "2000-5000" = 6, "5000 - 10000" = 8, "10000 - 50000" = 10, "50000 - 100000" = 12, ">100000" = 18)

my_formula <- y ~ x

#GDP
p1<-ggplot(result_df4, aes(log(GDP_MD), accuracy, color=CONTINENT)) +
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

p1


#per!!
result_country <- fread('/root/autodl-tmp/humPoulResult/data/result_country.csv')
#healthData <- fread('/root/autodl-tmp/humPoulResult/data/health expenditure.csv',header = T,drop=c(5:44,66:69))
globalSHP <- vect('/root/autodl-tmp/worldBorder/ne_10m_admin_0_countries.shp')
countryData <- as.data.frame(globalSHP) %>% .[,c('NAME_LONG','POP_EST','POP_RANK','ECONOMY','INCOME_GRP','CONTINENT','GDP_MD')]
result_df3 <- left_join(result_country,countryData,by=c('country'='NAME_LONG'))
#result_df4 <- left_join(result_df3,healthData,by=c('country'='Country Name')) %>% na.omit()
result_df4$per <- (result_df4$df/result_df4$alldf)


# 
pop_quantiles <- quantile(result_df4$POP_EST, probs=seq(0, 1, by=0.1), na.rm = TRUE)
# 
result_df4$POP_EST_size <- findInterval(result_df4$POP_EST, vec = pop_quantiles)
result_df4$POP_EST <- as.numeric(result_df4$POP_EST)
result_df4 <- mutate(result_df4, POP_EST_label = cut(POP_EST, breaks = c(0, 10000000, 20000000, 50000000, 100000000, 500000000, 1000000000, 1397800000),
                                                     labels = c("<1000", "1000 - 2000", "2000-5000", "5000 - 10000", "10000 - 50000", "50000 - 100000", ">100000")))

# 
result_df4 <- result_df4[order(result_df4$CONTINENT, -result_df4$POP_EST), ]

# 
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
p2<-ggplot(result_df5, aes(per, accuracy, color=CONTINENT)) +
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


