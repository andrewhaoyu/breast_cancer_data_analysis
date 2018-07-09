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

# p.adjunt <- cbind(p.fixed,p.fixed.adjust,p.random,p.random.adjust)
# write.csv(p.adjunt,file=paste0("./known_SNPs/known_SNPs_analysis_G_revised/additive_model_third_order/result/padjunt.csv"))











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
library(data.table)
i1 = 1
data1 <- fread("./data/iCOGS_euro_v10_10232017.csv",header=T)
data1 <- as.data.frame(data1)
y.pheno.mis1 <- cbind(data1$Behaviour1,data1$ER_status1,data1$PR_status1,data1$HER2_status1,data1$Grade1)
colnames(y.pheno.mis1) = c("Behavior","ER","PR","HER2","Grade")
z.standard <- GenerateZstandard(y.pheno.mis1)
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
library(xlsx)
generate_second_stage_parameter_names = function(tumor_characteristics){
  result = c("baseline effect (95%CI)",
             "P_value for baseline effect")
  for(i in 1:length(tumor_characteristics)){
    result = c(result,paste0(tumor_characteristics[i]," main effect(95%CI)"),
               paste0(tumor_characteristics[i]," main effect P_Value"))
  }
  result = c(result,
             "Triple-negative odds ratio",
             "Triple-negative p value",
             "Wald global test p value",
             "Wald global heterogneity test p value",
             #"Score global test p value"
             "Mixed Model global test p value",
             "Mixed Model global heterogeneity test p value"
             )
            
  return(result)
}

result <-  NULL
p_mglobal <- NULL
p_mheter <- NULL
#first.stage <- NULL

for(i1 in 1:178){
  print(i1)
  load(paste0("./known_SNPs/known_SNPs_analysis_G_revised/additive_model_third_order/result/heter_result_",i1,".Rdata"))
  result <- rbind(result,heter.result[[1]])
  p_mglobal <- c(p_mglobal,heter.result[[2]][1])
  p_mheter <- c(p_mheter,heter.result[[3]][1])
                    
  #first.stage <- rbind(first.stage,heter.result[[2]])
}
result <- cbind(result,p_mglobal,p_mheter)
tumor.characteristics <- c("ER","PR","HER2","Grade")
generate_second_stage_parameter_names(tumor.characteristics)
colnames(result) <- generate_second_stage_parameter_names(tumor.characteristics)



p.wald.assoc <- as.vector(result[,13])
p.wald.assoc[is.na(p.wald.assoc)] <- 0
p.wald.assoc.adjust <- p.adjust(p.wald.assoc,method="BH")
p.wald.heter <- result[,14]
p.wald.heter.adjust <- p.adjust(p.wald.heter,method="BH")
p.mixed.assoc.baselineER.fixed <- result[,15]
p.mixed.assoc.baselineER.fixed.adjust <- 
  p.adjust(p.mixed.assoc.baselineER.fixed,
           method="BH")
p.mixed.heter.baselineER.fixed <- result[,16]
p.mixed.heter.baselineER.fixed.adjust <- 
  p.adjust(p.mixed.heter.baselineER.fixed,
           method="BH")


pvalue = data.frame(p.wald.assoc,
                    p.wald.assoc.adjust,
                    p.wald.heter,
                    p.wald.heter.adjust,
                    p.mixed.assoc.baselineER.fixed,
                    p.mixed.assoc.baselineER.fixed.adjust,
                    p.mixed.heter.baselineER.fixed,
                    p.mixed.heter.baselineER.fixed.adjust)

colnames(pvalue) = c("Wald global test p value",
                     "Wald global test p value (BH adjust)",
                     "Wald global heterogneity test p value",
                     "Wald global heterogneity test p value (BH adjust)",
                     
                     "Mixed Model global test p value",
                     "Mixed Model global test p value",
                     "Mixed Model global heterogeneity test p value",
                     "Mixed Model global heterogeneity test p value")

result <- result[,-c(13:16)]

names.temp <- colnames(result)
p.value.names <- colnames(pvalue)
full.names <- c(names.temp,p.value.names)
result <- data.frame(result,pvalue)

colnames(result) <- full.names
write.csv(result,file = paste0("./known_SNPs/known_SNPs_analysis_G_revised/additive_model_third_order/result/known_third_order.csv"))


