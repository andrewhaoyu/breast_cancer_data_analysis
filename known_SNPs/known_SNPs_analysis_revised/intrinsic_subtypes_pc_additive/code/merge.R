
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
i1 = 1
library(bc2)
setwd("/spin1/users/zhangh24/breast_cancer_data_analysis/")
data1 <- read.csv("./data/iCOGS_euro_v10_05242017.csv",header=T)
y.pheno.mis1 <- cbind(data1$Behaviour1,data1$PR_status1,data1$ER_status1,data1$HER2_status1)
colnames(y.pheno.mis1) = c("Behavior","PR","ER","HER2")
# Grade1.fake <- data1$Grade1
# Grade1.fake[data1$Grade1==2|data1$Grade1==3] <- 1
# Grade1.fake[data1$Grade1==1] <- 0
#y.pheno.mis1 <- cbind(data1$Behaviour1,data1$PR_status1,data1$ER_status1,data1$HER2_status1,Grade1.fake)
# y.pheno.mis1 <- cbind(data1$Behaviour1,data1$PR_status1,data1$ER_status1,data1$HER2_status1)

x.test.all.mis1 <- data1[,c(27:206)]

x.covar.mis1 <- data1[,5:14]
x.all.mis1 <- as.matrix(cbind(x.test.all.mis1[,i1],x.covar.mis1))
colnames(x.all.mis1)[1] <- "gene"

Heter.result.Icog = EMmvpoly(y.pheno.mis1,baselineonly = NULL,additive = x.all.mis1,pairwise.interaction = NULL,saturated = NULL,missingTumorIndicator = 888)
z.standard <- Heter.result.Icog[[12]]
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


















setwd("/spin1/users/zhangh24/breast_cancer_data_analysis/known_SNPs/known_SNPs_analysis_revised/intrinsic_subtypes_pc_additive/result")
library(xlsx)
generate_self_design_second_stage_parameter_names = function(tumor_characteristics){
  result = NULL
  for(i in 1:length(tumor_characteristics)){
    result = c(result,paste0(tumor_characteristics[i]," odds ratio(95%CI)"),
               paste0(tumor_characteristics[i]," P_Value"))
  }
  result = c(result,"Wald global test p value",
             "Wald global heterogneity test p value",
             "Score global test p value",
             "Mixed Model global test p value",
             "loglikelihood",
             "AIC")
  return(result)
}

result <-  NULL
first.stage <- NULL


for(i in 1:181){
  print(i)
  load(paste0("heter_result_",i,".Rdata"))
  result <- rbind(result,heter.result[[1]])
  first.stage <- rbind(first.stage,heter.result[[2]])
}

tumor.characteristics <- c("Luminal A","Luminal B","HER2 Enriched","Triple Neg")
generate_self_design_second_stage_parameter_names(tumor.characteristics)

colnames(result) <- generate_self_design_second_stage_parameter_names(tumor.characteristics)

result <- as.data.frame(result)
tumor.characteristics <- c("PR","ER","HER2")

colnames(first.stage) = generate_first_stage_parameter_names(tumor.characteristics,z.standard)

write.xlsx(result,file="./intrinsic_subtypes_model.xlsx",sheetName="intrinsic_subtypes_2nd_stage")
write.xlsx(first.stage,file="./intrinsic_subtypes_model.xlsx",sheetName="
           intrinsic_subtypes_1st_stage",append=T)
