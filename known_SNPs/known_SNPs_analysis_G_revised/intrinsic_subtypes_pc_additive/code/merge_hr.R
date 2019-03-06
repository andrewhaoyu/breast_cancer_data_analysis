setwd("/spin1/users/zhangh24/breast_cancer_data_analysis/known_SNPs/known_SNPs_analysis_G_revised/intrinsic_subtypes_pc_additive/result/")
library(xlsx)
result <-  NULL
first.stage <- NULL

for(i in 1:178){
  print(i)
  load(paste0("heter_result_hr",i,".Rdata"))
  
  result <- rbind(result,DisplaySecondStageTestResult(heter.result[[1]],
                                                      heter.result[[2]]))
}
result <- result[,c(1:6)]


