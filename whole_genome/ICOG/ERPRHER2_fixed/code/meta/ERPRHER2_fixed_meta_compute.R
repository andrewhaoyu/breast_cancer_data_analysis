args = commandArgs(trailingOnly = T)
i1 = as.numeric(args[[1]])

library(bc2)
load(paste0("/data/zhangh24/breast_cancer_data_analysis/whole_genome/ICOG/ERPRHER2_fixed/result/meta/meta_result_shared_1p_sub",i1,".Rdata"))
load(paste0("/data/zhangh24/breast_cancer_data_analysis/whole_genome/ICOG/ERPRHER2_fixed/result/meta/icog_onco_score_infor_sub",i1,".Rdata"))



second.num <- 4

n <- nrow(meta_result_shared_1p_sub)

pvalue <- rep(0,n)

for(i in 1:n){
  print(i)
  icog_onco_infor_oneline <- icog_onco_score_infor_sub[i,]
  pvalue[i] <- MetaPfunction(icog_onco_infor_oneline,second.num)
  
}

meta_result_shared_1p_sub <- cbind(meta_result_shared_1p_sub,pvalue)
save(meta_result_shared_1p_sub,
     file=paste0("/data/zhangh24/breast_cancer_data_analysis/whole_genome/ICOG/ERPRHER2_fixed/result/meta/meta_result_shared_1p_sub",i1,".Rdata"))






