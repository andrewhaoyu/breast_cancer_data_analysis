setwd("/data/zhangh24/breast_cancer_data_analysis/")
data <- read.csv("./poland/result/known_snps/known_snps_result.csv",header=T)
idx <- which(data$standard_p<=0.001)
data.zoom <- data[-idx,]

data.plot <- data.zoom[,c(1,18,17,19)]
data.plot$log_standard <- -log10(data.plot$standard_p)
data.plot$log_mix_two_stage <- -log10(data.plot$mixed.model.global.test.for.association.baseline.and.ER.fixed.)
data.plot$heter <- factor(ifelse(data.plot$mixed.model.global.test.for.heterogeneity.baseline.and.ER.fixed.<=0.05,1,0))

library(ggplot2)
library(ggthemes)

png("./poland/result/known_snps/figure/Mixed two-stage model versus standard analysis.png",height=6,width = 8,units="in",res=400)
ggplot(data.plot,aes(x=log_standard,y=log_mix_two_stage,color=heter))+
  geom_point()+
  #theme_minimal()+
  xlab("standard analysis -log10(p-value)")+
  ylab("mixed two-stage model -log10(p-value)")+
  ggtitle("mixed two-stage model versus standard analysis")+
  scale_fill_Publication()+
  scale_colour_manual(name="Global heterogeneity test",values=c("#386cb0","#fdb462"),
          labels=c("p-value>=0.05","p-value<0.05"))+
  theme_Publication()
dev.off()
  
  
  
bp <- ggplot(data=PlantGrowth, aes(x=group, y=weight, fill=group)) + geom_boxplot()
bp + scale_fill_discrete(name="Experimental\nCondition")
