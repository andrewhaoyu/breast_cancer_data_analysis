#Goal: permutate phenotypes to verify batch effect
#update date: 01172020
#install_github("andrewhaoyu/bc2", args = c('--library="/home/zhangh24/R/x86_64-pc-linux-gnu-library/3.5"'))
###1 represent Icog
###2 represent Onco

rm(list=ls())
#commandarg <- commandArgs(trailingOnly=F)
#myarg <- commandarg[length(commandarg)]
#myarg <- sub("-","",myarg)
#i1 <- as.numeric(myarg)
args = commandArgs(trailingOnly = T)
i1 = as.numeric(args[[1]])

setwd("/data/zhangh24/breast_cancer_data_analysis/")
library(readr)
library(devtools)
library(CompQuadForm)
library(bc2, lib.loc ="/home/zhangh24/R/x86_64-pc-linux-gnu-library/3.6/")
library(data.table)
set.seed(i1)
if(i1<=177){
  ##analysis for Icog
  data1 <- fread("./data/iCOGS_euro_v10_10232017.csv",header=T)
  data1 <- as.data.frame(data1)
  y.pheno.mis1 <- cbind(data1$Behaviour1,data1$ER_status1,data1$PR_status1,data1$HER2_status1,data1$Grade1)
  colnames(y.pheno.mis1) = c("Behavior","ER","PR","HER2","Grade")
  n <- nrow(y.pheno.mis1)
  idx.case <- which(y.pheno.mis1[,1]==1)
  idx <- sample(idx.case,length(idx.case),replace = FALSE)
  y.pheno.mis1[idx.case,c(2:5)] <- y.pheno.mis1[idx,c(2:5)]
  # Grade1.fake <- data1$Grade1
  # Grade1.fake[data1$Grade1==2|data1$Grade1==3] <- 1
  # Grade1.fake[data1$Grade1==1] <- 0
  #y.pheno.mis1 <- cbind(data1$Behaviour1,data1$PR_status1,data1$ER_status1,data1$HER2_status1,Grade1.fake)
  # y.pheno.mis1 <- cbind(data1$Behaviour1,data1$PR_status1,data1$ER_status1,data1$HER2_status1)
  
  
  
  
  
  
  x.test.all.mis1 <- data1[,c(27:203)]
  x.test.all.mis1 <- x.test.all.mis1
  ###pc1-10 and age
  x.covar.mis1 <- data1[,c(5:14)]
  
  idx.control <- which(y.pheno.mis1[,1]==0)
  maf <- sum(x.test.all.mis1[idx.control,i1])/(2*length(idx.control))
  if(maf>=0.5){
    x.test.all.mis1[,i1] < 2 - x.test.all.mis1[,i1]
  }
 
  
  x.all.mis1 <- as.matrix(cbind(x.test.all.mis1[,i1],x.covar.mis1))
  colnames(x.all.mis1)[1] <- "gene"
  
  Heter.result.Icog = TwoStageModel(y.pheno.mis1,baselineonly = NULL,additive = x.all.mis1[,1:11],pairwise.interaction = NULL,saturated = NULL,missingTumorIndicator = 888)
  z.standard <- Heter.result.Icog[[12]]
  z.additive.design <- as.matrix(cbind(1,z.standard))
  M <- nrow(z.standard)
  number.of.tumor <- ncol(z.standard)
  log.odds.icog <- Heter.result.Icog[[1]][(M+1):(M+1+number.of.tumor)]
  
  sigma.log.odds.icog <- Heter.result.Icog[[2]][(M+1):(M+1+number.of.tumor),(M+1):(M+1+number.of.tumor)]
  beta.icog <- z.additive.design%*%log.odds.icog
  beta.sigma.icog <- z.additive.design%*%sigma.log.odds.icog%*%t(z.additive.design)
  loglikelihood.icog <- Heter.result.Icog[[8]]
  rm(Heter.result.Icog)
  
  #analysis for Onco Array
  #data2 <- read.csv("./V10/Onco_euro_v10_05242017.csv",header=T)
  data2 <- fread("./data/Onco_euro_v10_10232017.csv",header=T)
  data2 <- as.data.frame(data2)
  y.pheno.mis2 <- cbind(data2$Behaviour1,data2$ER_status1,data2$PR_status1,data2$HER2_status1,data2$Grade1)
  #y.pheno.mis2 <- cbind(data2$Behaviour1,data2$PR_status1,data2$ER_status1,data2$HER2_status1)
  colnames(y.pheno.mis2) = c("Behaviour","ER",
                             "PR","HER2","Grade")
  n <- nrow(y.pheno.mis2)
  idx.case <- which(y.pheno.mis2[,1]==1)
  idx <- sample(idx.case,length(idx.case),replace = FALSE)
  y.pheno.mis2[idx.case,c(2:5)] <- y.pheno.mis2[idx,c(2:5)]
  x.test.all.mis2 <- data2[,c(27:203)]
  #x.test.all.mis2 <- 2-x.test.all.mis2
  x.covar.mis2 <- data2[,c(5:14)]
  
  if(maf>=0.5){
    x.test.all.mis2[,i1] <- 2-x.test.all.mis2[,i1]
  }
  
  
  x.all.mis2 <- as.matrix(cbind(x.test.all.mis2[,i1],x.covar.mis2))
  colnames(x.all.mis2)[1] = "gene"
  
  
  
  Heter.result.Onco = TwoStageModel(y.pheno.mis2,baselineonly = NULL,additive = x.all.mis2,pairwise.interaction = NULL,saturated = NULL,missingTumorIndicator = 888)
  z.standard <- Heter.result.Onco[[12]]
  M <- nrow(z.standard)
  number.of.tumor <- ncol(z.standard)
  log.odds.onco <- Heter.result.Onco[[1]][(M+1):(M+1+number.of.tumor)]
  nparm <- length(Heter.result.Onco)
  sigma.log.odds.onco <- Heter.result.Onco[[2]][(M+1):(M+1+number.of.tumor),(M+1):(M+1+number.of.tumor)]
  loglikelihood.onco <- Heter.result.Onco[[8]]
  #rm(Heter.result.Onco)  
  
  
  meta.result <- LogoddsMetaAnalysis(log.odds.icog,
                                     sigma.log.odds.icog,
                                     log.odds.onco,
                                     sigma.log.odds.onco)
  
  second.stage.logodds.meta <- meta.result[[1]]
  second.stage.sigma.meta <- meta.result[[2]]
  
  
  
  test.result.second.wald <- DisplaySecondStageTestResult(second.stage.logodds.meta,second.stage.sigma.meta)
  
  heter.result <- list(data.frame(test.result.second.wald))
  save(heter.result,file=paste0("./known_SNPs/known_SNPs_analysis_G_revised/batch_effect_test/result/heter_result_permu_pheno_",i1,".Rdata"))
  
  
  
}else{
  data2 <- fread("./data/Onco_euro_v10_10232017.csv",header=T)
  data2 <- as.data.frame(data2)
  names2 = colnames(data2)
  y.pheno.mis2 <- cbind(data2$Behaviour1,data2$ER_status1,data2$PR_status1,data2$HER2_status1,data2$Grade1)
  #y.pheno.mis2 <- cbind(data2$Behaviour1,data2$PR_status1,data2$ER_status1,data2$HER2_status1)
  colnames(y.pheno.mis2) = c("Behaviour","ER",
                             "PR","HER2","Grade")
  idxi1 = which(names2=="rs554219")
  
  x.test.all.mis2 <- data2
  x.covar.mis2 <- data2[,c(5:14)]
  
  n <- nrow(y.pheno.mis2)
  idx.case <- which(y.pheno.mis2[,1]==1)
  idx <- sample(idx.case,length(idx.case),replace = FALSE)
  y.pheno.mis2[idx.case,c(2:5)] <- y.pheno.mis2[idx,c(2:5)]
  
  x.all.mis2 <- as.matrix(cbind(x.test.all.mis2[,idxi1],x.covar.mis2))
  colnames(x.all.mis2)[1] = "gene"
  
  
  
  
  Heter.result.Onco = TwoStageModel(y.pheno.mis2,baselineonly = NULL,additive = x.all.mis2,pairwise.interaction = NULL,saturated = NULL,missingTumorIndicator = 888)
  z.standard <- Heter.result.Onco[[12]]
  z.additive.design <- as.matrix(cbind(1,z.standard))
  M <- nrow(z.standard)
  number.of.tumor <- ncol(z.standard)
  log.odds.onco <- Heter.result.Onco[[1]][(M+1):(M+1+number.of.tumor)]
  nparm = length(Heter.result.Onco[[1]])
  sigma.log.odds.onco <- Heter.result.Onco[[2]][(M+1):(M+1+number.of.tumor),(M+1):(M+1+number.of.tumor)]
  beta.onco <- z.additive.design%*%log.odds.onco
  beta.sigma.onco <- z.additive.design%*%sigma.log.odds.onco%*%t(z.additive.design)
  loglikelihood.onco <- Heter.result.Onco[[8]]
  
  
  M <- nrow(z.standard)
  number.of.tumor <- ncol(z.standard)
  log.odds.onco <- Heter.result.Onco[[1]][(M+1):(M+1+number.of.tumor)]
  sigma.log.odds.onco <- Heter.result.Onco[[2]][(M+1):(M+1+number.of.tumor),(M+1):(M+1+number.of.tumor)]
  
  second.stage.logodds.meta <-log.odds.onco
  second.stage.sigma.meta <- sigma.log.odds.onco
  
  
  
  
  
  test.result.second.wald <- DisplaySecondStageTestResult(second.stage.logodds.meta,second.stage.sigma.meta)
  
  
  heter.result <- list(data.frame(test.result.second.wald))
  
  save(heter.result,file=paste0("./known_SNPs/known_SNPs_analysis_G_revised/batch_effect_test/result/heter_result_permu_pheno_",i1,".Rdata"))
  
  
  
  
}



