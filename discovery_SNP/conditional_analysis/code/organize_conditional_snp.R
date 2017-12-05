load("/spin1/users/zhangh24/breast_cancer_data_analysis/discovery_SNP/conditional_analysis/result/conditional.snp.list.icog.Rdata")
load("/spin1/users/zhangh24/breast_cancer_data_analysis/discovery_SNP/conditional_analysis/result/conditional.snp.list.onco.Rdata")
load("/spin1/users/zhangh24/breast_cancer_data_analysis/discovery_SNP/conditional_analysis/result/all.conditional.snps.Rdata")

n.icog <- length(conditional.snp.list.icog[[1]])
n.onco <- length(conditional.snp.list.onco[[1]])

idx.match.icog <- match(all.conditional.snps$SNP.ICOGS,conditional.snp.list.icog[[1]])

test.icog <- conditional.snp.list.icog[[1]][idx.match.icog]
test.icog.value <- conditional.snp.list.icog[[2]][,idx.match.icog]

conditional.snp.list.icog.clean <- list(test.icog,
                                        test.icog.value)
idx.match.onco <- match(all.conditional.snps$SNP.ONCO,conditional.snp.list.onco[[1]])

test.onco <- conditional.snp.list.onco[[1]][idx.match.onco]
test.onco.value <- conditional.snp.list.onco[[2]][,idx.match.onco]
conditional.snp.list.onco.clean <- list(test.onco,
                                        test.onco.value)
save(conditional.snp.list.icog.clean,file="/spin1/users/zhangh24/breast_cancer_data_analysis/discovery_SNP/conditional_analysis/result/conditional.snp.list.icog.clean.Rdata")
save(conditional.snp.list.onco.clean,file="/spin1/users/zhangh24/breast_cancer_data_analysis/discovery_SNP/conditional_analysis/result/conditional.snp.list.onco.clean.Rdata")




load("/spin1/users/zhangh24/breast_cancer_data_analysis/discovery_SNP/conditional_analysis/result/conditional.snp.list.icog.clean.Rdata")
load("/spin1/users/zhangh24/breast_cancer_data_analysis/discovery_SNP/conditional_analysis/result/conditional.snp.list.onco.clean.Rdata")
load("/spin1/users/zhangh24/breast_cancer_data_analysis/discovery_SNP/conditional_analysis/result/all.conditional.snps.Rdata")



n.condition <- nrow(all.conditional.snps)
library(bc2)


for(i in 1:1000){
  print(i)
  start.end <- startend(n.condition,1000,i)
  start <- start.end[1]
  end <- start.end[2]
  conditional.snp.list.icog.clean.sub <- list(conditional.snp.list.icog.clean[[1]][start:end],conditional.snp.list.icog.clean[[2]][,start:end])
  conditional.snp.list.onco.clean.sub <- list(conditional.snp.list.onco.clean[[1]][start:end],
                                              conditional.snp.list.onco.clean[[2]][,start:end])
  save(conditional.snp.list.icog.clean.sub,file=paste0("/spin1/users/zhangh24/breast_cancer_data_analysis/discovery_SNP/conditional_analysis/result/conditional.snp.list.icog.clean.sub",i,".Rdata"))
  save(conditional.snp.list.onco.clean.sub,file=paste0("/spin1/users/zhangh24/breast_cancer_data_analysis/discovery_SNP/conditional_analysis/result/conditional.snp.list.onco.clean.sub",i,".Rdata"))
}