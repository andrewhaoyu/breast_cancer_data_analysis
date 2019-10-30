data <- read.csv("/data/zhangh24/breast_cancer_data_analysis/discovery_SNP/conditional_analysis/result/condition.all.csv",header=T,stringsAsFactors = F)
data <- data[,-1]

library(tidyverse)
data.con <- data %>% filter(mark!="known snp")
Alleles <- data.con$Alleles
n.str <- length(Alleles)
Alleles.str <- strsplit(Alleles,":")
rs.new <- rep("c",n.str)
Alleles.new <- rep(NA,n.str)
for(i in 1:n.str){
    print(i)
    if(length(Alleles.str[[i]]==4)){
    rs.new[i] <- Alleles.str[[i]][1]
    Alleles.new[i] <- paste0(Alleles.str[[i]][3],"/",Alleles.str[[i]][4])  
  }else{
    rs.new[i] <- Alleles.str[[i]][1]
  }
}
idx.con <- which(data$mark!="known snp")
data$rs.id[idx.con] <- rs.new
data$Alleles[idx.con] <- Alleles.new

n.d <- nrow(data)

plac <- 2
for(j in c(7,9,11,13,15)){
  print(j)
  for(i in 1:n.d){
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
write.csv(data,file="/data/zhangh24/breast_cancer_data_analysis/discovery_SNP/conditional_analysis/result/condition.all.Rclean.csv")
