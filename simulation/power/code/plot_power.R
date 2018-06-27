setwd('/Users/zhangh24/GoogleDrive/breast_cancer_data_analysis/simulation/power/result')
data <- read.csv("power_result.csv")
colnames(data)
library(ggplot2)
ggplot(data,aes(x = "Sample.Size",y="Poewr"))+
  geom_bar()
