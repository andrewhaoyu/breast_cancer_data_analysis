setwd("/spin1/users/zhangh24/breast_cancer_data_analysis/")
p.value <- rep(0,178)
p.fixed <- rep(0,178)
p.value.adjust <- rep(0,178)
for(i1 in 1:178){
  load(paste0("./known_SNPs/known_SNPs_analysis_G_revised/additive_model_third_order/result/heter_result_",i1,".Rdata"))
  p.value[i1] <- heter.result[[2]][1]
  p.fixed[i1] <- heter.result[[1]][14]
}
temp <- cbind(p.value,p.fixed)
p.value.adjust <- p.adjust(p.value,method="BH")
p.adjunt <- cbind(p.value,p.value.adjust)
write.csv(p.adjunt,file=paste0("./known_SNPs/known_SNPs_analysis_G_revised/additive_model_third_order/result/padjunt.csv"))
