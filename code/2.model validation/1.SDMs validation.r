#nowater
nowater <- fread('/root/autodl-tmp/result/auc_scoreLsit_nowater_xgb.txt',header = F)
nowater2 <- data.table(
  spName=str_extract(nowater$V1,'(?<=\')\\D+(?=maxent|xgb|rf)'),
  model=str_extract(nowater$V1,'maxent|xgb|rf'),
  month=str_extract(nowater$V1,'(?<=maxent|xgb|rf)\\d+(?=\')'),
  auc=str_extract(nowater$V1,'(?<=: ).*')%>% as.numeric(),
  label='nowater'
)
# nowater2$index <- paste0(nowater2$spName,'-',nowater2$month,'-',nowater2$model)
#有水体
water_xgb <- fread('/root/autodl-tmp/result/auc_scoreLsit_water_xgb.txt',header = F)
water_xgb2 <- data.table(
  spName=str_extract(water_xgb$V1,'(?<=\')\\D+(?=maxent|xgb|rf)'),
  model=str_extract(water_xgb$V1,'maxent|xgb|rf'),
  month=str_extract(water_xgb$V1,'(?<=maxent|xgb|rf)\\d+(?=\')'),
  auc=str_extract(water_xgb$V1,'(?<=: ).*') %>% as.numeric(),
  label='water_xgb'
)

# water_xgb2$index <- paste0(water_xgb2$spName,'-',water_xgb2$month,'-',water_xgb2$model)
allData <- rbind(nowater2,water_xgb2)
#拼接属性
allDf2 <- fread('/root/autodl-tmp/result/allDf2.csv')
allData2 <- left_join(allData,allDf2,by=c('spName'='latName')) %>% subset(is.na(order)==F)


#画图
allData3 <- allData2[auc>0,.(auc=mean(auc)),.(spName,label.x)]
#整体对比
ggplot(data=allData3,aes(label.x,auc))+
  geom_boxplot(aes(fill=label.x),outlier.alpha = 0.1,outlier.fill = 'grey')+
  stat_compare_means(
    label="p.signif",
    show.legend = F)+
  theme_bw()
#分模型对比
allData4 <- subset(allData2,allData2$auc>0)
  # allData2[auc>0,.(auc=mean(auc)),.(spName,label.x,model)]
ggplot(data=allData4,aes(model,auc))+
  geom_boxplot(aes(fill=label.x),outlier.alpha = 0.1,outlier.fill = 'grey')+
  stat_compare_means(aes(group=label.x),
                     label="p.signif",
                     show.legend = F)+
  theme_bw()

#分目对比
allData5 <- subset(allData2,allData2$auc>0)
  # allData2[auc>0,.(auc=mean(auc)),.(order,label.x)]
ggplot(data=allData5,aes(x=order,y=auc,fill=label.x))+
  geom_boxplot()+
  # geom_col(position='dodge2')+
  # coord_cartesian(ylim = c(0.6, 0.9))+
  stat_compare_means(aes(group=label.x),method = 't.test',
                     label="p.signif",
                     show.legend = F)+
  theme_bw()+
  theme(
    axis.text.x = element_text(angle=60,hjust = 0.9)
  )

