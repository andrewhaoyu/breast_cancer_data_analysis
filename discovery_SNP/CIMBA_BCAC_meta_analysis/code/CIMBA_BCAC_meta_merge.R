##########merge the CIMBA and BCAC meta analysis result together
setwd('/spin1/users/zhangh24/breast_cancer_data_analysis/')
load(paste0("./discovery_SNP/CIMBA_BCAC_meta_analysis/result/CIMBA.BCAC.combine.Rdata"))
total <- nrow(CIMBA.BCAC.combine)
meta.result <- matrix(0,total,30)
total <- 0
for(i1 in 1:1000){
  print(i1)
  load(paste0("./discovery_SNP/CIMBA_BCAC_meta_analysis/result/meta.result.sub",i1,".Rdata"))
  temp <- nrow(meta.result.sub)
  meta.result[total+(1:temp),] <- meta.result.sub
  total <- total+temp
}

CIMBA.BCAC.combine[,c(22:51)] <- meta.result
CIMBA.BCAC.meta.result <- CIMBA.BCAC.combine

save(CIMBA.BCAC.meta.result,file=paste0("./discovery_SNP/CIMBA_BCAC_meta_analysis/result/CIMBA.BCAC.meta.result.Rdata"))


############organize the data into genetic correlation analysis format
