setwd("/data/zhangh24/breast_cancer_data_analysis/")
result <- NULL
for(i1 in 1:178){
  load(paste0("./known_SNPs/known_SNPs_analysis_G_revised/batch_effect_test/result/heter_result_",i1,".Rdata"))
  result <- rbind(result,heter.result[[1]])
}
save(result,file ="./known_SNPs/known_SNPs_analysis_G_revised/batch_effect_test/result/permutate_genoytp_result.rdata")

