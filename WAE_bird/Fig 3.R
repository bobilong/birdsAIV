
#######Fig3. new cattle#########################
#a-------------
crs <- '+proj=longlat +datum=WGS84'

globalRaster <- rast(vals=1:259200,nrows=360, ncols=720,xmin=-180, xmax=180,ymin=-90, ymax=90,crs=crs)

popd2015<-rast("/root/autodl-tmp/root/autodl-tmp/WAEdata_new_y/result/otherdata/pop.tif") %>% resample(globalRaster)
Entropy<-rast("/root/autodl-tmp/root/autodl-tmp/WAEdata_new_y/result/WAE.tif")
poul2015<-rast("/root/autodl-tmp/root/autodl-tmp/WAEdata_new_y/result/otherdata/poultry.tif")%>% resample(globalRaster)
cattle2015<-rast("/root/autodl-tmp/root/autodl-tmp/WAEdata_new_y/result/otherdata/cattle.tif") %>% resample(globalRaster)


quantiles_pop <- global(popd2015,quantile,probs=seq(0, 1, 0.01),na.rm=T)
quantiles_entr<- global(Entropy,quantile,probs=seq(0, 1, 0.01),na.rm=T)
quantiles_poul<- global(poul2015,quantile,probs=seq(0, 1, 0.01),na.rm=T)
quantiles_cattle<- global(cattle2015,quantile,probs=seq(0, 1, 0.01),na.rm=T)

quan_pop <- 60  #15%
quan_entr<- 4.11  #
quan_poul <- 10000  #
quan_cattle <- 700 #

pop_new <- ifel(popd2015>quan_pop,1,0)
entr_new <- ifel(Entropy>quan_entr,10,0)
poul_new <- ifel(poul2015>quan_poul,2,0)
cattle_new <- ifel(cattle2015>quan_cattle,4,0)

hotentrpoppoulcat <- pop_new+entr_new+poul_new+cattle_new   ;plot(hotentrpoppoulcat)
names(hotentrpoppoulcat) <-'type'
hotAND <- ifel(hotentrpoppoulcat==17,1,0)  ; plot(hotAND); global(hotAND,sum,na.rm=T)
#writeRaster(hotentrpoppoulcat, "/root/autodl-tmp/root/autodl-tmp/WAEdata_new_y/result/Hot_data/hotentrpoppoulcat.tif")

hot_11 <- ifel(hotentrpoppoulcat==11,1,0)  ; plot(hot_11); global(hot_11,sum,na.rm=T)
hot_12 <- ifel(hotentrpoppoulcat==12,1,0)  ; plot(hot_12); global(hot_12,sum,na.rm=T)
hot_13 <- ifel(hotentrpoppoulcat==13,1,0)  ; plot(hot_13); global(hot_13,sum,na.rm=T)
hot_14 <- ifel(hotentrpoppoulcat==14,1,0)  ; plot(hot_14); global(hot_14,sum,na.rm=T)
hot_15 <- ifel(hotentrpoppoulcat==15,1,0)  ; plot(hot_15); global(hot_15,sum,na.rm=T)
hot_16 <- ifel(hotentrpoppoulcat==16,1,0)  ; plot(hot_16); global(hot_16,sum,na.rm=T)
hot_17 <- ifel(hotentrpoppoulcat==17,1,0)  ; plot(hot_17); global(hot_17,sum,na.rm=T)


countrys <- vect ('/root/autodl-tmp/root/autodl-tmp/WAEdata_new_y/result/otherdata/Con_EU_dis.shp')

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
# 
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
    legend.position = "none",  # 
    legend.text = element_text(size= 12),
    axis.text = element_text(size = 12),
    panel.grid.major = element_blank()
  )


#b.-------------
plot(hotentrpoppoulcat)

hot4_df<-c(hotentrpoppoulcat, popd2015 ,cattle2015, poul2015) %>% terra::as.data.frame(.,xy=T) %>%na.omit()
head(hot4_df)
names(hot4_df)<-c("x","y","type","pop","cattle","poul")

# 
type_counts <- hot4_df %>% 
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




custom_colors <- c(#"0" = "white", "1" = "white","2" = "white","3" = "white","4" = "white",
  #"5" = "white","6" = "white","7" = "white","10" = "white",
  "0" = "lightgrey",
  "11" = "#FE97A4", "12" = "#71A3F1", "14" = "#EBBF00",  "13" = "#EA68A2", 
  "15" = "#EB7100", "16" = "#9BB31D", "17" = "#E53341")
labels <- c("0" = " ",
            "12" = "Poultry to WAE", 
            "11" = "Human to WAE", 
            "14" = "Cattle to WAE", 
            "16" = "Cattle-Poultry to WAE", 
            "13" = "Human-Poultry to WAE", 
            "15" = "Human-Cattle to WAE", 
            "17" = "Human-Poultry-Cattle to WAE")
# 
ggplot(Num_result, aes(x = factor(1), y = pop, fill = factor(type))) +
  geom_bar(width = 1, stat = "identity") +
  coord_polar(theta = "y") +
  scale_fill_manual(values = custom_colors, labels = labels) +
  theme_void() +
  labs(fill = "Type")+
  theme(
    legend.position = "none")  # 

ggplot(Num_result, aes(x = factor(1), y = poul, fill = factor(type))) +
  geom_bar(width = 1, stat = "identity") +
  coord_polar(theta = "y") +
  scale_fill_manual(values = custom_colors, labels = labels) +
  theme_void() +
  labs(fill = "Type") +
  ggtitle("Poultry Distribution")

ggplot(Num_result, aes(x = factor(1), y = cattle, fill = factor(type))) +
  geom_bar(width = 1, stat = "identity") +
  coord_polar(theta = "y") +
  scale_fill_manual(values = custom_colors, labels = labels) +
  theme_void() +
  labs(fill = "Type") +
  ggtitle("Cattle Distribution")










