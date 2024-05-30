
library(data.table)
library(dplyr)

basePath<-"/root/autodl-tmp/humPoulResult/data/"

world.map <- rnaturalearth::ne_countries(returnclass = "sf") |> dplyr::filter(continent != "Antarctica")
globalCountry <- vect(world.map) 
globalRaster <- rast(vals=1:259200,nrows=360, ncols=720,xmin=-180, xmax=180,ymin=-90, ymax=90,crs=crs)
coast <- ne_coastline(scale = "small", returnclass = "sf")
crs <- '+proj=longlat +datum=WGS84'


######a, b------------

Global_birdEntropy_df<-fread(paste0(basePath,'Global_birdEntropy_df.csv'))
Global_CHEntropy_df<-fread(paste0(basePath,'Global_CHEntropy_df.csv'))


Contry5c <-c('China','India','EU','Nigeria','USA')

Global_birdEntropy_df5_long <- gather(Global_birdEntropy_df, key = "Variable", value = "Value",
                                      -`Population(Log)`, -`Poltry(Log)`, -`pop`, -`poul`, -`Country`) %>%na.omit()
Global_birdEntropy_df5_long <- subset(Global_birdEntropy_df5_long, Value>0)
summary(Global_birdEntropy_df5_long)#ALL

Global_CHEntropy_df5_long <- gather(Global_CHEntropy_df, key = "Variable", value = "Value",
                                    -`Population(Log)`, -`Poltry(Log)`, -`pop`, -`poul`, -`Country`) %>%na.omit()
Global_CHEntropy_df5_long <- subset(Global_CHEntropy_df5_long, Value>0)
summary(Global_CHEntropy_df5_long)#CH



######corelation
All5_poppoulresults_df <- data.frame(Variable = character(),  b = numeric(),  SE = numeric(), b_SE = character(), z = numeric(),  P = numeric(),  R2 = numeric(), 
                                     Country = character(), Type = character(), Part = character()  ,stringsAsFactors = FALSE)

Contry5c <-c('China','India','EU','Nigeria','USA')
parts <- c('pop','poul')
types <- c('All', 'CH')
variables <- c('Seabirds', 'Shorebirds', 'Waterfowls', 'Large wading birds', 'Others')

# 
for(con in Contry5c){
  
  for(typ in types){
    if(typ == 'All'){
      
      for(part in parts){
        
        for (var in variables) {
          
          con_birdEntropy_df <- subset(Global_birdEntropy_df5_long,Global_birdEntropy_df5_long$Country==con)
          vardata = subset(con_birdEntropy_df,con_birdEntropy_df$Variable==var)
          
          if(part =='pop'){
            model <- lm(`Population(Log)` ~ Value, vardata)
          }
          else{
            model <- lm(`Poltry(Log)` ~ Value, vardata)
            
          }
          model_summary <- summary(model)
          
         
          b <- model_summary$coefficients[2, 1]
          SE <- model_summary$coefficients[2, 2]
          b_SE <- sprintf("%.2f ± %.2f", b, SE)
          z <- b / SE
          P <- model_summary$coefficients[2, 4]
          R2 <- model_summary$r.squared
          
          
          All5_poppoulresults_df <- rbind(All5_poppoulresults_df, data.frame(Variable = var, b = b, SE = SE, b_SE = b_SE, z = z, P = P, R2 = R2, 
                                                                             Country = con, Type = typ, Part = part))
        }
      }
    }
    
    if(typ == 'CH'){
      
      for(part in parts){
        
        for (var in variables) {
          
          con_CHEntropy_df <- subset(Global_CHEntropy_df5_long,Global_CHEntropy_df5_long$Country==con)
          vardata = subset(con_CHEntropy_df,con_CHEntropy_df$Variable==var)
          
          if(part =='pop'){
            model <- lm(`Population(Log)` ~ Value, vardata)
          }
          else{
            model <- lm(`Poltry(Log)` ~ Value, vardata)
            
          }
          model_summary <- summary(model)
          
          
          b <- model_summary$coefficients[2, 1]
          SE <- model_summary$coefficients[2, 2]
          b_SE <- sprintf("%.2f ± %.2f", b, SE)
          z <- b / SE
          P <- model_summary$coefficients[2, 4]
          R2 <- model_summary$r.squared
          
          
          All5_poppoulresults_df <- rbind(All5_poppoulresults_df, data.frame(Variable = var, b = b, SE = SE, b_SE = b_SE, z = z, P = P, R2 = R2, 
                                                                             Country = con, Type = typ, Part = part))
        }
        
      }
      
    }
    
  }
}


