basePath<-"/root/autodl-tmp/root/autodl-tmp/WAEdata_new_y/result/"
crs <- '+proj=longlat +datum=WGS84'
world.map <- rnaturalearth::ne_countries(returnclass = "sf") |> dplyr::filter(continent != "Antarctica")
globalCountry <- vect(world.map) 
globalRaster <- rast(vals=1:259200,nrows=360, ncols=720,xmin=-180, xmax=180,ymin=-90, ymax=90,crs=crs)
coast <- ne_coastline(scale = "small", returnclass = "sf")

# speciesPixelNumPath <- list.files('/root/autodl-tmp/humPoulResult/data/single_model',pattern = '.tif',full.names = T)
# spName <- basename(speciesPixelNumPath) %>% str_sub(.,1,-5)                      
# speciesPixelNumPath2 <- speciesPixelNumPath[spName%in%allDf$LatName]

allDf <- fread('/root/autodl-tmp/root/autodl-tmp/WAEdata_new_y/result/otherdata/allDf779_reclass.csv')
speciesPixelNumPath <- list.files('/root/autodl-tmp/root/autodl-tmp/WAEdata_new_y/singleRast',pattern = '.tif',full.names = T)

spName <- basename(speciesPixelNumPath) %>% str_sub(.,1,-5)                      
speciesPixelNumPath5 <- speciesPixelNumPath[spName%in%allDf$LatName]


#Fig.4 Calculate AE for functional groups###############
#CH
allDf2 <- fread('/root/autodl-tmp/root/autodl-tmp/WAEdata_new_y/result/otherdata/allDf779_reclass.csv')
allDf2CH <- allDf2[allDf2$`Host`== "Confirmed Host",]
spName <- basename(speciesPixelNumPath5) %>% str_sub(.,1,-5)

for (fuc in unique(allDf2CH$`Functional Group`)) {
  fucName <- allDf2CH[allDf2CH$`Functional Group`==fuc,]
  paths <- speciesPixelNumPath5[spName%in%fucName$LatName]
  monthNum <- rast(paths) %>% sum(na.rm=T)
  calEntropy <- lapply(paths, function(x){
    r <- rast(x) %>% sum(.,na.rm=T)
    pi <- r/monthNum
    y <- -pi*log(pi)
    names(y) <- str_sub(basename(x),1,-5)
    return(y)
  })
  actEntropy <- rast(calEntropy) %>% sum(na.rm = T)
  writeRaster(actEntropy,paste0('/root/autodl-tmp/root/autodl-tmp/WAEdata_new_y/result/WAE_data/WAE_',fuc,'_CH.tif'))
}


#All
allDf2 <- fread('/root/autodl-tmp/root/autodl-tmp/WAEdata_new_y/result/otherdata/allDf779_reclass.csv')
#allDf2CH <- allDf2[allDf2$`Host`== "Confirmed Host",]
spName <- basename(speciesPixelNumPath5) %>% str_sub(.,1,-5)

for (fuc in unique(allDf2$`Functional Group`)) {
  fucName <- allDf2[allDf2$`Functional Group`==fuc,]
  paths <- speciesPixelNumPath5[spName%in%fucName$LatName]
  monthNum <- rast(paths) %>% sum(na.rm=T)
  calEntropy <- lapply(paths, function(x){
    r <- rast(x) %>% sum(.,na.rm=T)
    pi <- r/monthNum
    y <- -pi*log(pi)
    names(y) <- str_sub(basename(x),1,-5)
    return(y)
  })
  actEntropy <- rast(calEntropy) %>% sum(na.rm = T)
  writeRaster(actEntropy,paste0('/root/autodl-tmp/root/autodl-tmp/WAEdata_new_y/result/WAE_data/WAE_',fuc,'_ALL.tif'))
}


#data prepare-----
#ALL
listFunc <- list.files('/root/autodl-tmp/root/autodl-tmp/WAEdata_new_y/result/WAE_data/', pattern = 'ALL', full.names = T)
AE_Func<-rast(listFunc)%>% resample(globalRaster)%>% mask(globalCountry) %>% crop(ext(globalRaster))

pop<-rast("/root/autodl-tmp/root/autodl-tmp/WAEdata_new_y/result/otherdata/pop.tif") %>% resample(globalRaster)%>% mask(globalCountry) %>% crop(ext(globalRaster))
poul<-rast("/root/autodl-tmp/root/autodl-tmp/WAEdata_new_y/result/otherdata/poultry.tif")%>% resample(globalRaster)%>% mask(globalCountry) %>% crop(ext(globalRaster))
cattle<-rast("/root/autodl-tmp/root/autodl-tmp/WAEdata_new_y/result/otherdata/cattle.tif")%>% resample(globalRaster)%>% mask(globalCountry) %>% crop(ext(globalRaster))

pop_log<- rast("/root/autodl-tmp/root/autodl-tmp/WAEdata_new_y/result/otherdata/pop_log.tif")%>% resample(globalRaster)%>% mask(globalCountry) %>% crop(ext(globalRaster))
poul_log<- rast("/root/autodl-tmp/root/autodl-tmp/WAEdata_new_y/result/otherdata/poul_log.tif")%>% resample(globalRaster)%>% mask(globalCountry) %>% crop(ext(globalRaster))
cattle_log<- rast("/root/autodl-tmp/root/autodl-tmp/WAEdata_new_y/result/otherdata/cattle_log.tif")%>% resample(globalRaster)%>% mask(globalCountry) %>% crop(ext(globalRaster))
values(pop_log)[is.infinite(values(pop_log))] <- NA
values(poul_log)[is.infinite(values(poul_log))] <- NA
values(cattle_log)[is.infinite(values(cattle_log))] <- NA
# globalSHP <- vect('/root/autodl-tmp/zyresult/Con_popentrpoul_sf_EU.shp')
# globalSHP2 <- terra::aggregate(globalSHP,'name_ec')
countrys <- vect ('/root/autodl-tmp/root/autodl-tmp/WAEdata_new_y/result/otherdata/Con_EU_dis.shp')
Con_5Raster <- subset(countrys, countrys$EU753 %in% c("United States", "EU", "China", "India"))%>% rasterize(.,globalRaster,field='EU753')


Global_birdEntropy_df <- c(AE_Func,pop_log,poul_log,cattle_log,pop,poul,cattle,Con_5Raster) %>% as.data.frame(.,xy=T)

names(Global_birdEntropy_df) <- c('x','y','Large wading birds','Others','Seabirds','Shorebirds','Waterfowl','pop_log','poul_log','cattle_log','pop','poul','cattle',"Country")
Global_birdEntropy_df <- Global_birdEntropy_df%>%
  dplyr::select('x','y','Seabirds','Shorebirds','Waterfowl','Large wading birds','Others','pop_log','poul_log','cattle_log','pop','poul','cattle',"Country")
summary(Global_birdEntropy_df)
#fwrite(Global_birdEntropy_df,'/root/autodl-tmp/root/autodl-tmp/WAEdata_new_y/result/Global_birdEntropy_df.csv')

#CH
listFuncCH <- list.files('/root/autodl-tmp/root/autodl-tmp/WAEdata_new_y/result/WAE_data/', pattern = 'CH', full.names = T)
AE_FuncCH<-rast(listFuncCH)%>% resample(globalRaster)%>% mask(globalCountry) %>% crop(ext(globalRaster))

pop<-rast("/root/autodl-tmp/root/autodl-tmp/WAEdata_new_y/result/otherdata/pop.tif") %>% resample(globalRaster)%>% mask(globalCountry) %>% crop(ext(globalRaster))
poul<-rast("/root/autodl-tmp/root/autodl-tmp/WAEdata_new_y/result/otherdata/poultry.tif")%>% resample(globalRaster)%>% mask(globalCountry) %>% crop(ext(globalRaster))
cattle<-rast("/root/autodl-tmp/root/autodl-tmp/WAEdata_new_y/result/otherdata/cattle.tif")%>% resample(globalRaster)%>% mask(globalCountry) %>% crop(ext(globalRaster))

pop_log<- rast("/root/autodl-tmp/root/autodl-tmp/WAEdata_new_y/result/otherdata/pop_log.tif")%>% resample(globalRaster)%>% mask(globalCountry) %>% crop(ext(globalRaster))
poul_log<- rast("/root/autodl-tmp/root/autodl-tmp/WAEdata_new_y/result/otherdata/poul_log.tif")%>% resample(globalRaster)%>% mask(globalCountry) %>% crop(ext(globalRaster))
cattle_log<- rast("/root/autodl-tmp/root/autodl-tmp/WAEdata_new_y/result/otherdata/cattle_log.tif")%>% resample(globalRaster)%>% mask(globalCountry) %>% crop(ext(globalRaster))
values(pop_log)[is.infinite(values(pop_log))] <- NA
values(poul_log)[is.infinite(values(poul_log))] <- NA
values(cattle_log)[is.infinite(values(cattle_log))] <- NA
# globalSHP <- vect('/root/autodl-tmp/zyresult/Con_popentrpoul_sf_EU.shp')
# globalSHP2 <- terra::aggregate(globalSHP,'name_ec')
countrys <- vect ('/root/autodl-tmp/root/autodl-tmp/WAEdata_new_y/result/otherdata/Con_EU_dis.shp')
Con_5Raster <- subset(countrys, countrys$EU753 %in% c("United States", "EU", "China", "India"))%>% rasterize(.,globalRaster,field='EU753')

Global_CHEntropy_df <-c(AE_FuncCH,pop_log,poul_log,cattle_log,pop,poul,cattle,Con_5Raster) %>% as.data.frame(.,xy=T)

names(Global_CHEntropy_df) <- c('x','y','Large wading birds','Others','Seabirds','Shorebirds','Waterfowl','pop_log','poul_log','cattle_log','pop','poul','cattle',"Country")
Global_CHEntropy_df <- Global_CHEntropy_df%>%
  dplyr::select('x','y','Seabirds','Shorebirds','Waterfowl','Large wading birds','Others','pop_log','poul_log','cattle_log','pop','poul','cattle',"Country")
summary(Global_CHEntropy_df)
#fwrite(Global_CHEntropy_df,'/root/autodl-tmp/root/autodl-tmp/WAEdata_new_y/result/Global_CHEntropy_df.csv')



#linear regression---------
Contry5c <- c("United States", "EU", "China", "India")

Global_birdEntropy_df5_long <- gather(Global_birdEntropy_df, key = "Variable", value = "Value",
                                      -`x`,-`y`,-`pop_log`, -`poul_log`, -`cattle_log`,-`pop`, -`poul`, -`cattle`, -`Country`) %>%na.omit()
Global_birdEntropy_df5_long <- subset(Global_birdEntropy_df5_long, Value>0)
summary(Global_birdEntropy_df5_long)#ALL

Global_CHEntropy_df5_long <- gather(Global_CHEntropy_df, key = "Variable", value = "Value",
                                    -`x`,-`y`,-`pop_log`, -`poul_log`, -`cattle_log`,-`pop`, -`poul`, -`cattle`, -`Country`) %>%na.omit()
Global_CHEntropy_df5_long <- subset(Global_CHEntropy_df5_long, Value>0)
summary(Global_CHEntropy_df5_long)#CH

All5_poppoulcatresults_df <- data.frame(Variable = character(),  b = numeric(),  SE = numeric(), b_SE = character(), z = numeric(),  P = numeric(),  R2 = numeric(), 
                                        Country = character(), Type = character(), Part = character()  ,stringsAsFactors = FALSE)

Contry5c <- c("United States", "EU", "China", "India")
parts <- c('pop','poul','cattle')
types <- c('All', 'CH')
variables <- c('Seabirds', 'Shorebirds', 'Waterfowl', 'Large wading birds', 'Others')

# 
for(con in Contry5c){#
  
  for(typ in types){
    if(typ == 'All'){#
      
      for(part in parts){#
        
        for (var in variables) {#
          
          con_birdEntropy_df <- subset(Global_birdEntropy_df5_long,Global_birdEntropy_df5_long$Country==con)
          vardata = subset(con_birdEntropy_df,con_birdEntropy_df$Variable==var)
          
          if(part =='pop'){#
            model <- lm(`pop_log` ~ Value, vardata)
          }
          if(part =='poul'){#
            model <- lm(`poul_log` ~ Value, vardata)
          }
          if(part =='cattle'){#cattle
            model <- lm(`cattle_log` ~ Value, vardata)
            
          }
          model_summary <- summary(model)
          
          # 
          b <- model_summary$coefficients[2, 1]
          SE <- model_summary$coefficients[2, 2]
          b_SE <- sprintf("%.2f ± %.2f", b, SE)
          z <- b / SE
          P <- model_summary$coefficients[2, 4]
          R2 <- model_summary$r.squared
          
          # 
          All5_poppoulcatresults_df <- rbind(All5_poppoulcatresults_df, data.frame(Variable = var, b = b, SE = SE, b_SE = b_SE, z = z, P = P, R2 = R2, 
                                                                                   Country = con, Type = typ, Part = part))
        }
      }
    }
    
    if(typ == 'CH'){#CH
      
      for(part in parts){#
        
        for (var in variables) {#
          
          con_CHEntropy_df <- subset(Global_CHEntropy_df5_long,Global_CHEntropy_df5_long$Country==con)
          vardata = subset(con_CHEntropy_df,con_CHEntropy_df$Variable==var)
          
          if(part =='pop'){#
            model <- lm(`pop_log` ~ Value, vardata)
          }
          if(part =='poul'){#
            model <- lm(`poul_log` ~ Value, vardata)
          }
          if(part =='cattle'){#cattle
            model <- lm(`cattle_log` ~ Value, vardata)
            
          }
          model_summary <- summary(model)
          
          # 
          b <- model_summary$coefficients[2, 1]
          SE <- model_summary$coefficients[2, 2]
          b_SE <- sprintf("%.2f ± %.2f", b, SE)
          z <- b / SE
          P <- model_summary$coefficients[2, 4]
          R2 <- model_summary$r.squared
          
          # 
          All5_poppoulcatresults_df <- rbind(All5_poppoulcatresults_df, data.frame(Variable = var, b = b, SE = SE, b_SE = b_SE, z = z, P = P, R2 = R2, 
                                                                                   Country = con, Type = typ, Part = part))
        }
        
      }
      
    }
    
  }
}


All5_poppoulcatresults_df


#fig----------
All5_poppoulcatresults_df<-fread('/root/autodl-tmp/root/autodl-tmp/WAEdata_new_y/result/All5_poppoulcatresults_df.csv')
lineardatalong <- All5_poppoulcatresults_df %>%data.table()
lineardatalong$r <- sqrt(lineardatalong$R2)
lineardatalong<-lineardatalong %>%
  mutate(Pd = cut(P, breaks = c(-Inf, 0.01, 0.05, Inf),
                  labels = c("< 0.01", "0.01 - 0.05", "> 0.05")))%>%
  # mutate(R2d = cut(R2, breaks = c(-Inf, 0.3, 0.4,0.5,0.6, Inf),
  #                  labels = c("< 0.3", "0.3 - 0.4", "0.4 - 0.5", "0.5 - 0.6", "> 0.6"))
  mutate(rd = cut(r, breaks = c(-Inf, 0.4,0.5,0.6,0.7, Inf),
                  labels = c("< 0.4", "0.4 - 0.5", "0.5 - 0.6", "0.6 - 0.7", "> 0.7")))
head(lineardatalong)

get_max_r2_rows <- function(dt, type) {
  max_r2 <- dt[Type == type, .(Max_r = max(r)), by = .(Country, Part)]    # 计算每个Country和Part组合的最大R2值
  setkey(dt, Country, Part, r)                # 将最大r值合并回原始数据集
  setkey(max_r2, Country, Part, Max_r)
  max_rows <- max_r2[dt, nomatch=0]      # 筛选出最大R2值对应的行
  #max_rows[, Max_r := NULL]    # Max_r
  return(max_rows)
}


max_r2_rows_ALL <- get_max_r2_rows(lineardatalong, "All") 
max_r2_rows_CH <- get_max_r2_rows(lineardatalong, "CH")
lineardata_max <- rbind(max_r2_rows_ALL, max_r2_rows_CH)
print(lineardata_max)

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


#colors_All <- c("< 0.4"="#F5EAD0", "0.4 - 0.5" = "#E3CB8F", "0.5 - 0.6" = "#BF8D44", "0.6 - 0.7" = "#8E581C", "> 0.7" = "#59310D")
colors_All <- c("< 0.4"="#FAE3D7", "0.4 - 0.5" = "#F2B396", "0.5 - 0.6" = "#D9705A", "0.6 - 0.7" = "#B42D34", "> 0.7" = "#6E0D20")
#colors_CH <- c("< 0.3"="#E7F2F1", "0.3 - 0.4" = "#AFDCD5", "0.4 - 0.5" = "#6DB2AD", "0.5 - 0.6" = "#347C78", "> 0.6" = "#174A41")
colors_CH <- c("< 0.4"="#FAE3D7", "0.4 - 0.5" = "#F2B396", "0.5 - 0.6" = "#D9705A", "0.6 - 0.7" = "#B42D34", "> 0.7" = "#6E0D20")


data_All <- lineardata_max[lineardata_max$Type == "All", ]
data_CH <- lineardata_max[lineardata_max$Type == "CH", ]


#colors_gradientALL <- c("#F5EAD0", "#E3CB8F", "#BF8D44", "#8E581C", "#59310D")
colors_gradientALL <-c("#FAE3D7", "#F2B396", "#D9705A", "#B42D34", "#6E0D20")
colors_gradientCH <- c("#FAE3D7", "#F2B396", "#D9705A", "#B42D34", "#6E0D20")
library(scales)
# 创建一个连续的图例
ggplot(data_All, aes(x=Country, y=Max_r, fill=Max_r)) +
  geom_bar(stat="identity", position=position_dodge()) +
  scale_fill_gradientn(colors=colors_gradientALL, 
                       values=rescale(c(0.4, 0.5, 0.6, 0.7, 0.8)),
                       name="r (All)",
                       limits=c(0.4, 0.8),
                       breaks=c(0.4, 0.5, 0.6, 0.7, 0.8),
                       labels=c("< 0.4", "0.5", "0.6", "0.7", "> 0.8")) +
  theme(text = element_text(size=10),
        legend.position = "bottom",
        legend.background = element_rect(fill="transparent"),
        legend.title = element_blank())

ggplot(data_CH, aes(x=Country, y=Max_r, fill=Max_r)) +
  geom_bar(stat="identity", position=position_dodge()) +
  scale_fill_gradientn(colors=colors_gradientCH, 
                       values=rescale(c(0.4, 0.5, 0.6, 0.7, 0.8)),
                       name="r (All)",
                       limits=c(0.4, 0.8),
                       breaks=c(0.4, 0.5, 0.6, 0.7, 0.8),
                       labels=c("< 0.4", "0.5", "0.6", "0.7", "> 0.8")) +
  theme(text = element_text(size=10),
        legend.position = "bottom",
        legend.background = element_rect(fill="transparent"),
        legend.title = element_blank())



#world map--------

#Con_popentrpoul_sf <- st_read("/root/autodl-tmp/zyresult/Con_popentrpoul_sf_EU.shp");head(Con_popentrpoul_sf)
Con_popentrpoul_sf <- vect ('/root/autodl-tmp/root/autodl-tmp/WAEdata_new_y/result/otherdata/Con_EU_dis.shp');head(Con_popentrpoul_sf)

Con_China<-subset(Con_popentrpoul_sf,Con_popentrpoul_sf$EU753=="China")
Con_India<-subset(Con_popentrpoul_sf,Con_popentrpoul_sf$EU753=="India") 
Con_EU<-subset(Con_popentrpoul_sf,Con_popentrpoul_sf$EU753=="EU")
Con_US<-subset(Con_popentrpoul_sf,Con_popentrpoul_sf$EU753=="United States") 

ggplot() +
  geom_sf(data = world.map, fill = "lightgrey", color = NA) + 
  geom_sf(data = Con_China, fill = "#baa1a1", color = NA) +    #lightcoral
  geom_sf(data = Con_India, fill = "#baa1a1", color = NA) + 
  geom_sf(data = Con_EU, fill = "#baa1a1", color = NA) + 
  #geom_sf(data = Con_Ethiopia, fill = "#baa1a1", color = NA) + 
  geom_sf(data = Con_US, fill = "#baa1a1", color = NA) + 
  coord_sf(crs = "+proj=longlat +datum=WGS84", xlim = c(-160, 165), ylim = c(-56, 90)) +
  theme_bw() + # 
  theme(legend.position = "none",
        axis.text = element_text(size=12),
        panel.grid.major = element_blank() ) # )

#little map--------
piedata_China <- lineardata_max[lineardata_max$Country == "China", ]
piedata_China$Value<-c(1,1,1,1,1,1)
color_China<-c("pop_All"="#B42D34","pop_CH"="#D9705A","poul_All"="#6E0D20","poul_CH"="#B42D34","cattle_All"="#FAE3D7","cattle_CH"="#FAE3D7")
#color_China<-c("pop_All"="lightgrey","pop_CH"="grey","poul_All"="lightgrey","poul_CH"="grey","cattle_All"="lightgrey","cattle_CH"="grey")

# 
desired_order <- c("pop_CH","cattle_All","cattle_CH", "poul_All", "poul_CH","pop_All")
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
    panel.background = element_rect(fill = "transparent"), # 
    plot.background = element_rect(fill = "transparent", color = NA), # 
    legend.position = "none"
  )
p
ggsave("/root/autodl-tmp/root/autodl-tmp/WAEdata_new_y/result/piedata_China.png", p, bg = "transparent", width = 10, height = 8, units = "in")







