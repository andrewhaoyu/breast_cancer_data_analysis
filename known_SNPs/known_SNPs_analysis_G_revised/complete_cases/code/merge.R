

setwd("/data/zhangh24/breast_cancer_data_analysis/")
library(readr)
library(devtools)
library(CompQuadForm)


result <- NULL
for(i1 in 1:178){
  print(i1)
  load(paste0("./known_SNPs/known_SNPs_analysis_G_revised/intrinsic_subtypes_pc_additive/result/heter_result_complete_",i1,".Rdata"))
  result <- rbind(result,DisplaySecondStageTestResult(heter.result[[1]],heter.result[[2]]))
  #first.stage <- rbind(first.stage,heter.result[[2]])
}

generate_self_design_second_stage_parameter_names = function(tumor_characteristics){
  result = NULL
  for(i in 1:length(tumor_characteristics)){
    result = c(result,paste0(tumor_characteristics[i]," odds ratio(95%CI)"),
               paste0(tumor_characteristics[i]," P_Value"))
  }
  result = c(result,"Wald global test p value",
             "Wald global heterogneity test p value")
  return(result)
}


tumor.characteristics <- c("Luminial A","Luminal B",
                           "Luminal B HER2Neg",
                           "HER2 Enriched",
                           "Triple Negative")
generate_self_design_second_stage_parameter_names(tumor.characteristics)

colnames(result) <- generate_self_design_second_stage_parameter_names(tumor.characteristics)[1:12]

result <- as.data.frame(result)
write.csv(result,file = "./known_SNPs/known_SNPs_analysis_G_revised/complete_cases/result/known_intrinsic_subtypes_complete_data.csv")



