setwd('/data/zhangh24/breast_cancer_data_analysis/')
snp.name <- rep("c",20)
p.value <- rep(0,20)
for(i1 in 1:20){
  load(paste0("./discovery_SNP/conditional_analysis_novel/result/standard_conditional_result",i1,".Rdata"))
p.value[i1] <- result[[1]]
snp.name[i1] <- result[[2]]
}
cbind(snp.name,p.value)
write.csv(cbind(snp.name,p.value),file="./discovery_SNP/conditional_analysis_novel/result/conditional_analysis_standard_result.csv")
