
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




Heter.result.Icog = TwoStageModel(y.pheno.mis1,baselineonly = NULL,additive = x.all.mis1,pairwise.interaction = NULL,saturated = NULL,missingTumorIndicator = 888,
                                  missingDataAlgorithm = "EM")

Heter.result.Icog = TwoStageModel(y.pheno.mis1,baselineonly = NULL,additive = x.all.mis1,pairwise.interaction = NULL,saturated = NULL,missingTumorIndicator = 888,
                                  missingDataAlgorithm = "OneStepMLE")

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






















setwd("/spin1/users/zhangh24/breast_cancer_data_analysis/known_SNPs/G_onemle/additive/result")
library(xlsx)
# generate_second_stage_parameter_names = function(tumor_characteristics){
#   result = c("baseline effect (95%CI)",
#              "P_value for baseline effect")
#   for(i in 1:length(tumor_characteristics)){
#     result = c(result,paste0(tumor_characteristics[i]," main effect(95%CI)"),
#                paste0(tumor_characteristics[i]," main effect P_Value"))
#   }
#   result = c(result,"Wald global test p value",
#              "Wald global heterogneity test p value",
#              "Score global test p value",
#              "Mixed Model global test p value (baseline fixed)",
#              "Mixed Model global heterogeneity test p value (baseline fixed)",
#              "Mixed Model global test p value (baseline+ER fixed)",
#              "Mixed Model global heterogeneity test p value (baseline+ER fixed)",
#              "loglikelihood",
#              "AIC")
#   return(result)
# }

generate_second_stage_parameter_names = function(tumor_characteristics){
  result = c("baseline effect (95%CI)",
             "P_value for baseline effect")
  for(i in 1:length(tumor_characteristics)){
    result = c(result,paste0(tumor_characteristics[i]," main effect(95%CI)"),
               paste0(tumor_characteristics[i]," main effect P_Value"))
  }
  result = c(result,"Wald global test p value",
             "Wald global heterogneity test p value")
  return(result)
}

result <-  NULL
#first.stage <- NULL

for(i in 1:178){
  print(i)
  load(paste0("heter_result_",i,".Rdata"))
  result <- rbind(result,heter.result[[1]])
  #first.stage <- rbind(first.stage,heter.result[[2]])
}

tumor.characteristics <- c("ER","PR","HER2","Grade")
generate_second_stage_parameter_names(tumor.characteristics)
colnames(result) <- generate_second_stage_parameter_names(tumor.characteristics)




result.temp <- result

p.value.onestepmle <- result.temp[,11]

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




Heter.result.Icog = TwoStageModel(y.pheno.mis1,baselineonly = NULL,additive = x.all.mis1,pairwise.interaction = NULL,saturated = NULL,missingTumorIndicator = 888,
                                  missingDataAlgorithm = "EM")
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






















setwd("/spin1/users/zhangh24/breast_cancer_data_analysis/known_SNPs/known_SNPs_analysis_G_revised/additive_model/result")
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

for(i in 1:178){
  print(i)
  load(paste0("heter_result_",i,".Rdata"))
  result <- rbind(result,heter.result[[1]])
  first.stage <- rbind(first.stage,heter.result[[2]])
}

tumor.characteristics <- c("ER","PR","HER2","Grade")
generate_second_stage_parameter_names(tumor.characteristics)
colnames(result) <- generate_second_stage_parameter_names(tumor.characteristics)


p.value.em <- result[,11]

temp.result <- cbind(p.value.em,p.value.onestepmle)
p.value.em <- p.value.em[-c(108:109)]
p.value.onestepmle <- p.value.onestepmle[-c(108:109)]
temp.result <- temp.result[-c(108:109),]
log10.p.value.em <- -log10(p.value.em)
log10.p.value.onestepmle <- -log10(p.value.onestepmle)
cbind(log10.p.value.em,log10.p.value.onestepmle)
log
plot(log10.p.value.em,log10.p.value.onestepmle,xlab="EM algorithm -log10(p)",
     ylab="one step mle -log10(p)",main="compare onestepmle to EM",xlim=c(0,10),ylim=c(0,10))
abline(a=0,b=1,color="red")

model1 <- lm(log10.p.value.onestepmle~log10.p.value.em)
summary(model1)
