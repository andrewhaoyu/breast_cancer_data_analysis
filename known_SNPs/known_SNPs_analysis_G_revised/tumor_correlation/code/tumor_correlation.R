#Goal: estimate the correlation between the tumor characteristics

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
  idx.complete1 <- which(y.pheno.mis1[,1]==1&
                           y.pheno.mis1[,2]!=888&
                           y.pheno.mis1[,3]!=888&
                           y.pheno.mis1[,4]!=888&
                           y.pheno.mis1[,5]!=888)
  y.pheno.complete1 <- y.pheno.mis1[idx.complete1,]
  table(y.pheno.complete1[,1])
  table(y.pheno.complete1[,2])
  table(y.pheno.complete1[,3])
  table(y.pheno.complete1[,4])
  table(y.pheno.complete1[,5])
  data2 <- fread("./data/Onco_euro_v10_10232017.csv",header=T)
  data2 <- as.data.frame(data2)
  y.pheno.mis2 <- cbind(data2$Behaviour1,data2$ER_status1,data2$PR_status1,data2$HER2_status1,data2$Grade1)
  idx.complete2 <- which(y.pheno.mis2[,1]==1&
                           y.pheno.mis2[,2]!=888&
                           y.pheno.mis2[,3]!=888&
                           y.pheno.mis2[,4]!=888&
                           y.pheno.mis2[,5]!=888)
  y.pheno.complete2 <- y.pheno.mis2[idx.complete2,]
  table(y.pheno.complete2[,1])
  table(y.pheno.complete2[,2])
  table(y.pheno.complete2[,3])
  table(y.pheno.complete2[,4])
  table(y.pheno.complete2[,5])
  y.pheno.complete.all <- rbind(y.pheno.complete1,
                                y.pheno.complete2)  
  correlation.result <- cor(y.pheno.complete.all[,2:5])  
  