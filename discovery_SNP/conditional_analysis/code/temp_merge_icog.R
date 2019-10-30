i1 = 2
load(paste0("/data/zhangh24/breast_cancer_data_analysis/discovery_SNP/conditional_analysis/result/conditional.snp.list.icog",i1,".Rdata"))
n.sub <- nrow(conditional.snp.list.icog[[2]])
load("/data/zhangh24/breast_cancer_data_analysis/discovery_SNP/conditional_analysis/result/all.conditional.snps.Rdata")

library(bc2)

extract.num <- nrow(all.conditional.snps)
snpid.result <- rep("c",extract.num)
library(bigmemory)
snpvalue.result <- matrix(0,n.sub,extract.num)

total <- 0
for(i1 in 1:564){
  print(i1)
  load(paste0("/data/zhangh24/breast_cancer_data_analysis/discovery_SNP/conditional_analysis/result/conditional.snp.list.icog",i1,".Rdata"))
  if(is.null(conditional.snp.list.icog)==T){
    temp <- 0
  }else{
    temp <- length(conditional.snp.list.icog[[1]])
    snpid.result[total+(1:temp)] <- conditional.snp.list.icog[[1]]
    snpvalue.result[,total+(1:temp)] <- conditional.snp.list.icog[[2]]
    total <- temp+total  
  }
}

snpid.result <- snpid.result[1:total]
snpvalue.result <- snpvalue.result[,1:total]
conditional.snp.list.icog <- list(snpid.result,snpvalue.result)
save(conditional.snp.list.icog,file=paste0("/data/zhangh24/breast_cancer_data_analysis/discovery_SNP/conditional_analysis/result/conditional.snp.list.icog.Rdata"))
