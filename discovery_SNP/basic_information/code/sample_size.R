
setwd("/spin1/users/zhangh24/breast_cancer_data_analysis/")

library(readr)
library(devtools)
library(CompQuadForm)
library(bc2)
library(data.table)
data1 <- fread("./data/iCOGS_euro_v10_10232017.csv",header=T)
data1 <- as.data.frame(data1)
y.pheno.mis1 <- cbind(data1$Behaviour1,data1$ER_status1,data1$PR_status1,data1$HER2_status1,data1$Grade1)
colnames(y.pheno.mis1) = c("Behavior","ER","PR","HER2","Grade")
data2 <- fread("./data/Onco_euro_v10_10232017.csv",header=T)
data2 <- as.data.frame(data2)
y.pheno.mis2 <- cbind(data2$Behaviour1,data2$ER_status1,data2$PR_status1,data2$HER2_status1,data2$Grade1)
icog_tumor <- cbind(table(data1$study,data1$Behaviour1),
      table(data1$study,data1$ER_status1),
      table(data1$study,data1$PR_status1),
      table(data1$study,data1$HER2_status1),
      table(data1$study,data1$Grade1))
rowSums(icog_tumor)
onco_tumor <- cbind(table(data2$study,data2$Behaviour1)
  ,table(data2$study,data2$ER_status1),
      table(data2$study,data2$PR_status1),
      table(data2$study,data2$HER2_status1),
      table(data2$study,data2$Grade1))
rowSums(onco_tumor)
tumor_sample_size <- rbind(icog_tumor,
                           rowSums(icog_tumor),
                           onco_tumor,
                           rowSums(onco_tumor),
                           rowSums(icog_tumor)+rowSums(onco_tumor)
      )
write.csv(tumor_sample_size,file = "")

#y.pheno.mis2 <- cbind(data2$Behaviour1,data2$PR_status1,data2$ER_status1,data2$HER2_status1)
colnames(y.pheno.mis2) = c("Behaviour","ER",
                           "PR","HER2","Grade")
