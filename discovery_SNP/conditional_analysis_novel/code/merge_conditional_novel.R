setwd('/spin1/users/zhangh24/breast_cancer_data_analysis/')
p.value.twostage <- rep(0,11)
p.value.standard <- rep(0,11)
for(i1 in 1:11){
  load(paste0("./discovery_SNP/conditional_analysis_novel/novel_conditional_reuslt",i1,".Rdata"))
  p.value.twostage[i1] <- result[[1]]
  p.value.standard[i1] <- result[[2]]
}
