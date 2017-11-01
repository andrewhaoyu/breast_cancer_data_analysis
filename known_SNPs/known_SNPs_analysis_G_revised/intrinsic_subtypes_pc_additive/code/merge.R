
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
y.pheno.mis1 <- cbind(data1$Behaviour1,data1$PR_status1,data1$ER_status1,data1$HER2_status1,data1$Grade1)
colnames(y.pheno.mis1) = c("Behavior","PR","ER","HER2","Grade")
# Grade1.fake <- data1$Grade1
# Grade1.fake[data1$Grade1==2|data1$Grade1==3] <- 1
# Grade1.fake[data1$Grade1==1] <- 0
#y.pheno.mis1 <- cbind(data1$Behaviour1,data1$PR_status1,data1$ER_status1,data1$HER2_status1,Grade1.fake)
# y.pheno.mis1 <- cbind(data1$Behaviour1,data1$PR_status1,data1$ER_status1,data1$HER2_status1)

x.test.all.mis1 <- data1[,c(27:203)]
###pc1-10 and age
x.covar.mis1 <- data1[,c(5:14,204)]

age <- data1[,204]
idx.complete <- which(age!=888)
x.all.mis1 <- as.matrix(cbind(x.test.all.mis1[,i1],x.covar.mis1))
colnames(x.all.mis1)[1] <- "gene"
y.pheno.mis1 <- y.pheno.mis1[idx.complete,]
x.all.mis1 <- x.all.mis1[idx.complete,]

#z.standard <- Heter.result.Icog[[12]]

Heter.result.Icog = EMmvpoly(y.pheno.mis1,baselineonly = NULL,additive = x.all.mis1,pairwise.interaction = NULL,saturated = NULL,missingTumorIndicator = 888)
z.standard <- Heter.result.Icog[[12]]



















setwd("/spin1/users/zhangh24/breast_cancer_data_analysis/known_SNPs/known_SNPs_analysis_G_revised/intrinsic_subtypes_pc_additive/result/")
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
             "Mixed Model heterogeneity test p value",
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



tumor.characteristics <- c("Luminial A","Luminal B",
                           "Luminal B HER2Neg - Luminal A",
                           "HER2 Enriched - Luminal A",
                           "Triple Negative - Luminal A")
generate_self_design_second_stage_parameter_names(tumor.characteristics)

colnames(result) <- generate_self_design_second_stage_parameter_names(tumor.characteristics)

result <- as.data.frame(result)


p.wald.assoc <- result[,11]
p.wald.assoc[which(is.na(p.wald.assoc))] <- 0
p.wald.assoc.adjust <- p.adjust(p.wald.assoc,method="BH")
p.wald.heter <- result[,12]
p.wald.heter.adjust <- p.adjust(p.wald.heter,method="BH")
p.score.assoc <- result[,13]
p.score.assoc.adjust <- p.adjust(p.score.assoc,method="BH")
p.mixed.assoc <- result[,14]
p.mixed.assoc.adjust <- 
  p.adjust(p.mixed.assoc,
           method="BH")
p.mixed.heter <- result[,15]
p.mixed.heter.adjust <- 
  p.adjust(p.mixed.heter,
           method="BH")

pvalue = data.frame(p.wald.assoc,
                    p.wald.assoc.adjust,
                    p.wald.heter,
                    p.wald.heter.adjust,
                    p.score.assoc,
                    p.score.assoc.adjust,
                    p.mixed.assoc,
                    p.mixed.assoc.adjust,
                    p.mixed.heter,
                    p.mixed.heter.adjust
)

colnames(pvalue) = c("Wald global test p value",
                     "Wald global test p value (BH adjust)",
                     "Wald global heterogneity test p value",
                     "Wald global heterogneity test p value (BH adjust)",
                     "Score global test p value",
                     "Score global test p value (BH adjust)",
                     "Mixed Model global test p value ",
                     "Mixed Model global test p value (BH adjust)",
                     "Mixed Model global heterogeneity test p value",
                     "Mixed Model global heterogeneity test p value (BH adjust)")

result <- result[,-c(11:15)]

result <- cbind(result[1:10],pvalue,result[,11:12])





generate_self_design_second_stage_parameter_names_2 = function(tumor_characteristics){
  result = NULL
  for(i in 1:length(tumor_characteristics)){
    result = c(result,paste0(tumor_characteristics[i]," odds ratio(95%CI)"),
               paste0(tumor_characteristics[i]," P_Value"))
  }
  return(result)
}

tumor.characteristics <- c("Luminal A","Luminal B","Luminal B HER2 Neg","HER2 Enriched","Triple Neg")
generate_self_design_second_stage_parameter_names_2(tumor.characteristics)

colnames(first.stage) = generate_self_design_second_stage_parameter_names_2(tumor.characteristics)







write.xlsx(result,file="./intrinsic_subtypes_model_G.xlsx",sheetName="intrinsic_subtypes_model_G_2nd_stage")
write.xlsx(first.stage,file="./intrinsic_subtypes_model_G.xlsx",sheetName="intrinsic_subtypes_model_G_1st_stage",append=T)


