
SubtypesTrans <- function(casecon,ER,PR,HER2,grade){
  n <- length(casecon)
  result <- rep("unknown",n)
  idx.con <- which(casecon==0)
  result[idx.con] <- "control"
  idx.LA <- which(HER2==0&(ER==1|PR==1)&(grade==1|grade==2))
  idx.LB <- which(HER2==1&(ER==1|PR==1))
  idx.LUBHER2 <- which(HER2==0&(ER==1|PR==1)&grade==3)
  idx.HER2 <- which(HER2==1&ER==0&PR==0)
  idx.TN <- which(HER2==0&ER==0&PR==0)
  #idx.mis <- which(HER2==888|ER==888|PR==888|Grade==888)
  result[idx.LA] <- "Luminal_A"
  result[idx.LB] <- "Luminal_B"
  result[idx.LUBHER2] <- "Luminal_B_HER2Neg"
  result[idx.HER2] <- "HER2Enriched"
  result[idx.TN] <- "TripleNeg"
  result <- factor(result,levels=c("control",
                                   "Luminal_A",
                                   "Luminal_B",
                                   "Luminal_B_HER2Neg",
                                   "HER2Enriched",
                                   "TripleNeg",
                                   "unknown"))
  
  return(result)
  
}


setwd("/data/zhangh24/breast_cancer_data_analysis/")

library(readr)
library(devtools)
library(CompQuadForm)
library(bc2)
library(data.table)
data1 <- fread("./data/iCOGS_euro_v10_10232017.csv",header=T)
data1 <- as.data.frame(data1)
y.pheno.mis1 <- cbind(data1$Behaviour1,data1$ER_status1,data1$PR_status1,data1$HER2_status1,data1$Grade1)
colnames(y.pheno.mis1) = c("Behavior","ER","PR","HER2","Grade")
intrinsic_icog <- SubtypesTrans(y.pheno.mis1[,1],
                                y.pheno.mis1[,2],
                                y.pheno.mis1[,3],
                                y.pheno.mis1[,4],
                                y.pheno.mis1[,5])
data2 <- fread("./data/Onco_euro_v10_10232017.csv",header=T)
data2 <- as.data.frame(data2)
y.pheno.mis2 <- cbind(data2$Behaviour1,data2$ER_status1,data2$PR_status1,data2$HER2_status1,data2$Grade1)
intrinsic_onco <- SubtypesTrans(y.pheno.mis2[,1],
                                y.pheno.mis2[,2],
                                y.pheno.mis2[,3],
                                y.pheno.mis2[,4],
                                y.pheno.mis2[,5])

icog_tumor <- cbind(table(data1$study,data1$Behaviour1),
      table(data1$study,data1$ER_status1),
      table(data1$study,data1$PR_status1),
      table(data1$study,data1$HER2_status1),
      table(data1$study,data1$Grade1))
icog_intrinsic_subtype_sample_size <- cbind(table(data1$study,data1$Behaviour1),
  table(data1$study,intrinsic_icog))

onco_intrinsic_subtype_sample_size <- cbind(table(data2$study,data2$Behaviour1),
                                            table(data2$study,intrinsic_onco))

intrinsic_sample_size <- rbind(icog_intrinsic_subtype_sample_size,
                               colSums(icog_intrinsic_subtype_sample_size),
                               onco_intrinsic_subtype_sample_size,
                               colSums(onco_intrinsic_subtype_sample_size),
                               colSums(icog_intrinsic_subtype_sample_size)  + colSums(onco_intrinsic_subtype_sample_size)                       
  
)
write.csv(intrinsic_sample_size,file = "./discovery_SNP/basic_information/result/intrinsic_sample_size.csv")

icog_intrinsic <- 
rowSums(icog_tumor)
onco_tumor <- cbind(table(data2$study,data2$Behaviour1)
  ,table(data2$study,data2$ER_status1),
      table(data2$study,data2$PR_status1),
      table(data2$study,data2$HER2_status1),
      table(data2$study,data2$Grade1))
rowSums(onco_tumor)
tumor_sample_size <- rbind(icog_tumor,
                           colSums(icog_tumor),
                           onco_tumor,
                           colSums(onco_tumor),
                           colSums(icog_tumor)+colSums(onco_tumor)
      )
write.csv(tumor_sample_size,file = "./discovery_SNP/basic_information/result/tumor_sample_size.csv")

#y.pheno.mis2 <- cbind(data2$Behaviour1,data2$PR_status1,data2$ER_status1,data2$HER2_status1)
colnames(y.pheno.mis2) = c("Behaviour","ER",
                           "PR","HER2","Grade")
