
change_binary_to_negative_positive = function(x){
  if(x==0){
    return("-")
  }else{
    return("+")
  }
}
generate_first_stage_parameter_names = function(tumor_characteristics,z_standard){
  max.z_standard = apply(z_standard,2,max)
  idx.not.binary = which(max.z_standard!=1)
  idx.binary = which(max.z_standard==1)
  result= NULL
  for(i in 1:nrow(z_standard)){
    names_each_row = NULL
    for(j in 1:ncol(z_standard)){
      if(j%in%idx.binary){
        temp = paste0(tumor_characteristics[j],change_binary_to_negative_positive(z_standard[i,j]))
      }else{
        temp = paste0(tumor_characteristics[j],z_standard[i,j])
      }
      names_each_row = paste0(names_each_row,temp)
    }
    result = c(result,paste0(names_each_row," OR (95%CI)"),paste("P_value for OR of ",names_each_row))
  }
  return(result)
}

setwd("/spin1/users/zhangh24/breast_cancer_data_analysis/")
library(readr)
library(devtools)
library(CompQuadForm)
library(bc2)
i1 = 1
data1 <- fread("./data/iCOGS_euro_v10_10232017.csv",header=T)
data1 <- as.data.frame(data1)
y.pheno.mis1 <- cbind(data1$Behaviour1,data1$PR_status1,data1$ER_status1,data1$HER2_status1,data1$Grade1)
colnames(y.pheno.mis1) = c("Behavior","PR","ER","HER2","Grade")
# Grade1.fake <- data1$Grade1
# Grade1.fake[data1$Grade1==2|data1$Grade1==3] <- 1
# Grade1.fake[data1$Grade1==1] <- 0
#y.pheno.mis1 <- cbind(data1$Behaviour1,data1$PR_status1,data1$ER_status1,data1$HER2_status1,Grade1.fake)
# y.pheno.mis1 <- cbind(data1$Behaviour1,data1$PR_status1,data1$ER_status1,data1$HER2_status1)

x.test.all.mis1 <- data1[,c(27:203)]
###pc1-10 and age
x.covar.mis1 <- data1[,c(5:14)]



x.all.mis1 <- as.matrix(cbind(x.test.all.mis1[,i1],x.covar.mis1))
colnames(x.all.mis1)[1] <- "gene"



Heter.result.Icog = TwoStageModel(y.pheno.mis1,baselineonly = NULL,additive = x.all.mis1,pairwise.interaction = NULL,saturated = NULL,missingTumorIndicator = 888,missingDataAlgorithm = "EM")
z.standard <- Heter.result.Icog[[12]]
#z.standard <- Heter.result.Icog[[12]]





















setwd("./known_SNPs/known_SNPs_analysis_G_revised/ERPosvesNeg/result/")
library(xlsx)
generate_self_design_second_stage_parameter_names = function(tumor_characteristics){
  # result = NULL
  # for(i in 1:length(tumor_characteristics)){
  #   result = c(result,paste0(tumor_characteristics[i]," odds ratio(95%CI)"),
  #              paste0(tumor_characteristics[i]," P_Value"))
  # }
  result = c("PR+ OR (95% CI)",
             "PR+ pvalue",
             "PR- OR (95% CI)",
             "PR- pvalue",
             "ER+ OR (95% CI)",
             "ER+ pvalue",
             "ER- OR (95% CI)",
             "ER- pvalue",
             "HER2+ OR (95% CI)",
             "HER2+ pvalue",
             "HER2- OR (95% CI)",
             "HER2- pvalue",
             "Grade1 OR (95% CI)",
             "Grade1 pvalue",
             "Grade2 OR (95% CI)",
             "Grade2 pvalue",
             "Grade3 OR (95% CI)",
             "Grade3 pvalue"
             )
  return(result)
}

result <-  NULL
first.stage <- NULL

for(i in 1:178){
  print(i)
  load(paste0("heter_result_",i,".Rdata"))
  result.temp <- heter.result[[1]][c(1:4,7:10,13:16,19:24)]
  result <- rbind(result,result.temp)
  #first.stage <- rbind(first.stage,heter.result[[2]])
}



tumor.characteristics <- c("Luminial A","Luminal B",
                           "Luminal B HER2Neg - Luminal A",
                           "HER2 Enriched - Luminal A",
                           "Triple Negative - Luminal A")
generate_self_design_second_stage_parameter_names(tumor.characteristics)

colnames(result) <- generate_self_design_second_stage_parameter_names(tumor.characteristics)

result <- as.data.frame(result)










write.xlsx(result,file="./binary_model.xlsx",sheetName="all binary results")
