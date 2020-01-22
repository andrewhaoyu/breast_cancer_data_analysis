#-------------------------------------------------------------------
# Update Date: 11/20/2018
# Create Date: 11/20/2018
# Goal: prepare base dataset for PRSice2 (decided not using this software)
# Author: Haoyu Zhang
#-------------------------------------------------------------------
setwd('/data/zhangh24/breast_cancer_data_analysis/risk_prediction/')
load("./EB_whole_genome/result/whole_gonome.rdata")

SNP <- whole_genome$rs_id

n <- nrow(whole_genome)
A1 <- rep(0,n)
A2 <- rep(0,n)
temp <- strsplit(whole_genome$ var_name,"_")
for(i in 1:n){
  print(i)
  
  A1[i] <- temp[[i]][4]
  A2[i] <- temp[[i]][3]
}
whole_genome <- cbind(whole_genome,A1,A2)
colnames(whole_genome)[c(56:57)] <- c("effect_allele","referece_allele")
save(whole_genome,file = "./EB_whole_genome/result/whole_gonome.rdata")
##save the effect allele to the ER data
##
library(pegas)
C <- C[]
