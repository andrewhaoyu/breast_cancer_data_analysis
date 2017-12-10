i1 = 2
load(paste0("/spin1/users/zhangh24/breast_cancer_data_analysis/discovery_SNP/conditional_analysis/result/conditional.snp.list.icog",i1,".Rdata"))
n.sub <- nrow(conditional.snp.list.icog[[2]])
load("/spin1/users/zhangh24/breast_cancer_data_analysis/discovery_SNP/conditional_analysis/result/all.conditional.snps.Rdata")

library(bc2)

extract.num <- nrow(all.conditional.snps)
snpid.result <- rep("c",extract.num)
library(bigmemory)
snpvalue.result <- big.matrix(n.sub,extract.num,init=0)

total <- 0
for(i1 in 1:564){
  print(i1)
  load(paste0("/spin1/users/zhangh24/breast_cancer_data_analysis/discovery_SNP/conditional_analysis/result/conditional.snp.list.icog",i1,".Rdata"))
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
snpvalue.result <- sub.big.matrix(snpvalue.result,firstCol=1,lastCol=total)
conditional.snp.list.icog <- list(snpid.result,snpvalue.result)
save(conditional.snp.list.icog,file=paste0("/spin1/users/zhangh24/breast_cancer_data_analysis/discovery_SNP/conditional_analysis/result/conditional.snp.list.icog.Rdata"))