SubtypesTrans <- function(casecon,ER,PR,HER2,grade){
  n <- length(casecon)
  result <- rep("unknown",n)
  idx.con <- which(casecon==0)
  result[idx.con] <- "control"
  idx.LA <- which((ER==1|PR==1)&HER2==0)
  idx.LB <- which((ER==1|PR==1)&HER2==1&grade!=3)
  idx.LB.HERneg <- which((ER==1|PR==1)&HER2==1&grade==3)
  idx.TN <- which(ER==0&PR==0&HER2==0)
  idx.HER2 <- which(ER==0&PR==0&HER2==1)
  result[idx.LA] <- "Luminal_A"
  result[idx.LB] <- "Luminal_B"
  result[idx.LB.HERneg] <- "Luminal_B_HER2neg"
  result[idx.HER2] <- "HER2_enriched"
  result[idx.TN] <- "triple_neg"
  return(result)
  
}
library(R.utils)
library(data.table)
setwd("/data/zhangh24/breast_cancer_data_analysis/")
data1 <- fread("./data/iCOGS_euro_v10_10232017.csv",header=T)
data1 <- as.data.frame(data1)
y.pheno.mis1 <- cbind(data1$Behaviour1,data1$ER_status1,data1$PR_status1,data1$HER2_status1,data1$Grade1)
data2 <- fread("./data/Onco_euro_v10_10232017.csv",header=T)
data2 <- as.data.frame(data2)
y.pheno.mis2 <- cbind(data2$Behaviour1,data2$ER_status1,data2$PR_status1,data2$HER2_status1,data2$Grade1)

y.pheno.mis <- rbind(y.pheno.mis1,
                     y.pheno.mis2)
idx <- which(y.pheno.mis[,1]==1&
               y.pheno.mis[,2]!=888&
               y.pheno.mis[,3]!=888&
               y.pheno.mis[,4]!=888&
               y.pheno.mis[,5]!=888)
y.pheno.case.complete <- y.pheno.mis[idx,]
result <- GenerateFreqTable(y.pheno.case.complete)
result[,5] <- result[,5]/length(idx)
write.csv(result,file = "/data/zhangh24/breast_cancer_data_analysis/data/BCAC_freq_result.csv")
