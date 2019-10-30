args = commandArgs(trailingOnly = T)
i1 = as.numeric(args[[1]])
setwd('/data/zhangh24/breast_cancer_data_analysis/')
load(paste0("./discovery_SNP/CIMBA_BCAC_meta_analysis/result/CIMBA.BCAC.combine.Rdata"))
library(bc2)
n <- nrow(CIMBA.BCAC.combine)
start.end <- startend(n,1000,i1)
start <- start.end[1]
end <- start.end[2]
total <- end-start+1
############## Meta analysis between BCAC and CIMBA
############## CIMBA only get the triple negative results
############## Triple negative is the 5th subtypes in BCAC
MetaCIMBABCAC <- function(BCAC.vector,CIMBA.vector){
  BCAC.log.odds <- BCAC.vector[1:5]
  BCAC.sigma <- matrix(BCAC.vector[6:30],5,5)
  CIMBA.log.odds <- CIMBA.vector[1]
  CIMBA.sigma <- CIMBA.vector[2]
  meta.log.odds <- BCAC.log.odds
  meta.log.odds[5] <- (BCAC.sigma[5,5]^-1+CIMBA.sigma^-1)^-1*
    (BCAC.sigma[5,5]^-1*BCAC.log.odds[5]+CIMBA.sigma^-1*CIMBA.log.odds)
  meta.sigma <- BCAC.sigma
  meta.sigma[5,5] = (BCAC.sigma[5,5]^-1+CIMBA.sigma^-1)^-1
  meta.sigma[1:4,5] <- meta.sigma[5,5]*
    (BCAC.sigma[5,5]^-1)*BCAC.sigma[1:4,5]
  meta.sigma[5,1:4] <- meta.sigma[1:4,5] 
  meta.vector <- c(meta.log.odds,as.vector(meta.sigma))
  return(meta.vector)
}
meta.result.sub <- matrix(0,total,30)
temp = 1
for(i in start:end){
  print(i)
  BCAC.vector <- as.numeric(CIMBA.BCAC.combine[i,21:50])
  CIMBA.vector <- as.numeric(c(CIMBA.BCAC.combine[i,5],as.numeric(CIMBA.BCAC.combine[i,6]^2)))
  meta.result.sub[temp,] <- MetaCIMBABCAC(BCAC.vector,CIMBA.vector)
  temp = temp+1
}
save(meta.result.sub,file = paste0("./discovery_SNP/CIMBA_BCAC_meta_analysis/result/meta.result.sub",i1,".Rdata"))
