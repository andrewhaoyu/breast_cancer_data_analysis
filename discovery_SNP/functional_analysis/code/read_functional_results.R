#Goal: read functional analysis results
setwd("/Users/zhangh24/GoogleDrive/breast_cancer_data_analysis/")
data <- read.csv("./data/Inquisit_result.csv",header=T,
                 stringsAsFactors = F)
length(table(data$Target.gene3))
length(table(data$Region2))
library(dplyr)
data.1 = data %>% filter(FINAL.INQ.SCORE.LEVEL5==1)
length(table(data.1$Target.gene3))
length(table(data.1$Region2))

data.2 = data %>% filter(FINAL.INQ.SCORE.LEVEL5==1
                         &Analysis1=="overall-analysis")
length(table(data.2$Target.gene3))
length(table(data.2$Region2))
