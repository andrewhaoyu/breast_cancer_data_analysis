setwd("./known_SNPs_analysis_revised/additive_model/result/")
library(xlsx)
generate_second_stage_parameter_names = function(tumor_characteristics){
  result = c("baseline effect (95%CI)",
             "P_value for baseline effect")
  for(i in 1:length(tumor_characteristics)){
    result = c(result,paste0(tumor_characteristics[i]," main effect(95%CI)"),
               paste0(tumor_characteristics[i]," main effect P_Value"))
  }
  result = c(result,"global test p value",
             "global heterogneity test p value")
  return(result)
}

result <-  NULL


for(i in 1:181){
  print(i)
  load(paste0("heter_result_",i,".Rdata"))
  result <- rbind(result,heter.result[[1]])
}

tumor.characteristics <- c("PR","ER","HER2")
generate_second_stage_parameter_names(tumor.characteristics)

colnames(result) <- generate_second_stage_parameter_names(tumor.characteristics)

result <- as.data.frame(result)

write.xlsx(result,file="./meta.xlsx",sheetName="second_stage")
