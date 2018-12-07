#-------------------------------------------------------------------
# Update Date: 12/07/2018
# Create Date: 12/07/2018
# Goal: merge the additive two-stage model results
# Author: Haoyu Zhang
#-------------------------------------------------------------------
#---------------------------------------#---------------------------------------
library(bc2)
places <- 2
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
# Grade1.fake <- data1$Grade1
# Grade1.fake[data1$Grade1==2|data1$Grade1==3] <- 1
# Grade1.fake[data1$Grade1==1] <- 0
#y.pheno.mis1 <- cbind(data1$Behaviour1,data1$PR_status1,data1$ER_status1,data1$HER2_status1,Grade1.fake)
# y.pheno.mis1 <- cbind(data1$Behaviour1,data1$PR_status1,data1$ER_status1,data1$HER2_status1)

x.test.all.mis1 <- data1[,c(27:203)]
###pc1-10 and age
x.covar.mis1 <- data1[,c(5:14)]



x.all.mis1 <- as.matrix(cbind(x.test.all.mis1[,i1],x.covar.mis1))




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














discovery_snp <- read.csv("/spin1/users/zhangh24/breast_cancer_data_analysis/data/discovery_snp_summary_new.csv",header=T)







setwd("/spin1/users/zhangh24/breast_cancer_data_analysis/discovery_SNP/additive_model/result")
library(xlsx)
generate_second_stage_parameter_names = function(tumor_characteristics){
  result = c("baseline effect (95%CI)",
             "P_value for baseline effect")
  for(i in 1:length(tumor_characteristics)){
    result = c(result,paste0(tumor_characteristics[i]," main effect(95%CI)"),
               paste0(tumor_characteristics[i]," main effect P_Value"))
  }
  result = c(result,"Wald global test p value",
             "Wald global heterogneity test p value",
             "Score global test p value",
             "Mixed Model global test p value (baseline fixed)",
             "Mixed Model global heterogeneity test p value (baseline fixed)",
             "Mixed Model global test p value (baseline+ER fixed)",
             "Mixed Model global heterogeneity test p value (baseline+ER fixed)",
             "loglikelihood",
             "AIC")
  return(result)
}

result <-  NULL
first.stage <- NULL

for(i in 1:35){
  print(i)
  load(paste0("heter_result_",i,".Rdata"))
  result <- rbind(result,heter.result[[1]])
  first.stage <- rbind(first.stage,heter.result[[2]])
}

tumor.characteristics <- c("ER","PR","HER2","Grade")
generate_second_stage_parameter_names(tumor.characteristics)
colnames(result) <- generate_second_stage_parameter_names(tumor.characteristics)
#take out the 5 SNPs after we changed the threshold to 5E-08
new_discovery_snp <- read.csv("/spin1/users/zhangh24/breast_cancer_data_analysis/data/discovery_snp_summary_new_Dec07.csv")
library(dplyr)
new.chr.pos = paste0(new_discovery_snp$CHR,":",new_discovery_snp$Pos)
old.chr.pos = as.character(discovery_snp$chr.pos)
idx <- which(old.chr.pos%in%
               new.chr.pos)
snp.names <- discovery_snp[idx,"SNP.ICOGS"]
result <- result[idx,]

p.wald.assoc <- result[,11]
p.wald.assoc[is.na(p.wald.assoc)] <- 0
p.wald.assoc.adjust <- p.adjust(p.wald.assoc,method="BH")
p.wald.heter <- result[,12]
p.wald.heter.adjust <- p.adjust(p.wald.heter,method="BH")
# p.score.assoc <- result[,13]
# p.score.assoc.adjust <- p.adjust(p.score.assoc,method="BH")
# p.mixed.assoc.baseline.fixed <- result[,14]
# p.mixed.assoc.baseline.fixed.adjust <- 
#   p.adjust(p.mixed.assoc.baseline.fixed,
#            method="BH")
# p.mixed.heter.baseline.fixed <- result[,15]
# p.mixed.heter.baseline.fixed.adjust <- 
#   p.adjust(p.mixed.heter.baseline.fixed,
#            method="BH")
p.mixed.assoc.baselineER.fixed <- result[,13]
p.mixed.assoc.baselineER.fixed.adjust <- 
  p.adjust(p.mixed.assoc.baselineER.fixed,
           method="BH")
p.mixed.heter.baselineER.fixed <- result[,14]
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
                     "Mixed Model global test p value (baseline+ER fixed)",
                     "Mixed Model global test p value (baseline+ER fixed, BH adjust)",
                     "Mixed Model global heterogeneity test p value (baseline+ER fixed)",
                     "Mixed Model global heterogeneity test p value (baseline+ER fixed BH adjust)")


result <- result[,-c(11:16)]
result <- data.frame(result,pvalue)

# result <- result[,-c(11:17)]
# 
# names.temp <- colnames(result)
# p.value.names <- colnames(pvalue)
# full.names <- c(names.temp[1:10],p.value.names,names.temp[11:12])
# result <- data.frame(result[1:10],pvalue,result[,11:12])
# 
# colnames(result) <- full.names


#result <- result[,-c(12,14,15:16,17:20,22,24)]









#colnames(first.stage) <- generate_first_stage_parameter_names(tumor.characteristics,z.standard)

#result <- data.frame(result)
row.names(result) <- snp.names
write.xlsx(result,file="./additive_model_G_age_2des.xlsx",sheetName="additive_model_2nd_stage_age_minor_allele")
#write.xlsx(first.stage,file="./additive_model.xlsx",sheetName="additive_model_1st_stage",append=T)

