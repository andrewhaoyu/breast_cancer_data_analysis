setwd("/data/zhangh24/breast_cancer_data_analysis/known_SNPs/known_SNPs_analysis_G_revised/intrinsic_subtypes/result")
library(xlsx)
# generate_second_stage_parameter_names = function(tumor_characteristics){
#   result = c("baseline effect (95%CI)",
#              "P_value for baseline effect")
#   for(i in 1:length(tumor_characteristics)){
#     result = c(result,paste0(tumor_characteristics[i]," main effect(95%CI)"),
#                paste0(tumor_characteristics[i]," main effect P_Value"))
#   }
#   result = c(result,"global test p value",
#              "global heterogneity test p value")
#   return(result)
# }

result <-  NULL
loglikelihood <- NULL
AIC <- NULL
o

for(i in 1:181){
  print(i)
  load(paste0("heter_result_",i,".Rdata"))
  result <- rbind(result,heter.result[[1]])
  loglikelihood <- c(loglikelihood,heter.result[[4]])
  AIC <- c(AIC,heter.result[[5]])
}



colnames(result)  <- c("Luminial A Odds","Luminial A P","Luminal B Odds",
                       "Luminal B P",
                                            "Luminal B HER2- Odds",
                       "Luminal B HER2- P",
                                            "HER2 Enriched Odds",
                       "HER2 Enriched P",
                                            "Triple Negative Odds",
                       "Triple Negative P",
                       "global association",
                       "global heterogeneity")

result <- data.frame(result,loglikelihood=loglikelihood,AIC=AIC)

write.xlsx(result,file="./meta.xlsx",sheetName="second_stage")
