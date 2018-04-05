setwd("/spin1/users/zhangh24/breast_cancer_data_analysis/")
p.random <- rep(0,178)
p.fixed <- rep(0,178)
p.random.adjust <- rep(0,178)
p.fixed.adjust <- rep(0,178)
for(i1 in 1:178){
  load(paste0("./known_SNPs/known_SNPs_analysis_G_revised/additive_model_third_order/result/heter_result_",i1,".Rdata"))
  p.random[i1] <- as.numeric(heter.result[[2]][1])
  p.fixed[i1] <- as.numeric(heter.result[[1]][14])
}
#temp <- cbind(p.value,p.fixed)
p.fixed.adjust <- p.adjust(p.fixed,method="BH")
p.random.adjust <- p.adjust(p.random,method="BH")

p.adjunt <- cbind(p.fixed,p.fixed.adjust,p.random,p.random.adjust)
write.csv(p.adjunt,file=paste0("./known_SNPs/known_SNPs_analysis_G_revised/additive_model_third_order/result/padjunt.csv"))
