#-------------------------------------------------------------------
# Update Date: 11/26/2018
# Create Date: 11/26/2018
# Goal: global heterogeneity test for whole genome
# Author: Haoyu Zhang
#-------------------------------------------------------------------

arg <- commandArgs(trailingOnly=T)
i1 <- as.numeric(arg[[1]])
library(bc2)
library(bcutility)
library(dplyr)
setwd('/spin1/users/zhangh24/breast_cancer_data_analysis/risk_prediction/')
load("./EB_whole_genome/result/whole_gonome.rdata")
n <- nrow(whole_genome)
size <- 1000
start.end <- startend(n,size,i1)
start <- start.end[1]
end <- start.end[2]
subtypes <- c("Luminial_A",
              "Luminal_B",
              "Luminal_B_HER2Neg",
              "HER2_Enriched",
              "TN")
log.odds <- whole_genome %>% select(subtypes)
covar <- whole_genome %>% select(24:48)

second.num <- 5

size <- 1000
start.end <- startend(n,size,i1)
start <- start.end[1]
end <- start.end[2]
p.heter.sub <- matrix(0,end-start+1,second.num)
temp = 1
for(j in start:end){
  
  print(j)
  log.odds.temp <-  as.numeric(log.odds[j,])
  covar.temp <-   matrix(
    as.numeric(covar[j,]),
    second.num,second.num)
  result.temp <- DisplaySecondStageTestResult(
    log.odds.temp,
    covar.temp,
    self.design = T)
  p.heter.sub[temp] <- result.temp[length(result.temp)]
  temp = temp+1
}

save(p.heter.sub,file=paste0("./risk_prediction/heter_whole_genome/result/p_heter_sub",i1,".Rdata"))

