data <- read.csv("/Users/zhangh24/GoogleDrive/breast_cancer_data_analysis/discovery_SNP/result/discovery_snp_analysis_results.csv",header=T,stringsAsFactors = F)
library(tidyverse)
n.data <- nrow(data)
plac <- 2
for(j in c(11,13,15,17)){
  print(j)
  for(i in 1:n.data){
    print(i)
    temp <- data[i,j]
    temp2<- strsplit(temp,"-")
    odds.high <- as.numeric(gsub(")","",temp2[[1]][2]))
    temp3 <- strsplit(temp2[[1]][1],"\\(")
    odds <- as.numeric(temp3[[1]][1])
    odds.low <- as.numeric(temp3[[1]][2])
    odds <- round(odds,plac)
    odds.low <- round(odds.low,plac)
    odds.high <- round(odds.high,plac)
    new <- paste0(odds,"(",odds.low,"-",odds.high,")")
    data[i,j] <- new
  }
}
write.csv(data,file="/Users/zhangh24/GoogleDrive/breast_cancer_data_analysis/discovery_SNP/result/discovery_snp_analysis_results_Rclean.csv")
