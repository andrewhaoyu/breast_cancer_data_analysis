setwd('/Users/zhangh24/GoogleDrive/breast_cancer_data_analysis/simulation/power/result')
data <- read.csv("power_result.csv")
colnames(data)
library(ggplot2)
idx <- which(data[,2]== "No heterogneity")
data.new <- data[idx,]
data.new[,3] <- as.factor(data.new[,3])

ggplot(data.new,aes(x=Sample.Size,y=Power,color=Method,
                shape = Scenarior ))+
  geom_line()+
  scale_x_discrete(c(5000,50000,100000))

+
  ylim(0,1)
#plot(data.new[,3],data.new[,4])
