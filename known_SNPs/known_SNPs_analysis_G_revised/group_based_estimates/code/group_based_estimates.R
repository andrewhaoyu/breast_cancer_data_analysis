#-------------------------------------------------------------------
# Update Date: 01/16/2020
# Goal: estimate log odds ratio and var of 178 SNPs for 5 intrinic subtypes for complete data
# Author: Haoyu Zhang
#-------------------------------------------------------------------
#---------------------------------------
result <- NULL
#install_github("andrewhaoyu/bc2", args = c('--library="/home/zhangh24/R/x86_64-pc-linux-gnu-library/3.4"'))
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
library(bcutility,lib.loc = "/spin1/home/linux/zhangh24/R/x86_64-pc-linux-gnu-library/3.6/")
library(data.table)

z.design <- matrix(c(
  c(0,1,1,1,0,0,0,0,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0),
  c(0,0,0,0,0,1,1,0,0,0,0,0,1,1,1,0,0,0,0,0,1,1,1),
  c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,0,0,0,0),
  c(0,0,0,0,1,0,0,0,0,0,0,1,0,0,0,0,0,0,0,1,0,0,0),
  c(1,0,0,0,0,0,0,1,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0)
),ncol=5)
z.design.score.baseline <- z.design[,1,drop=F]
z.design.score.casecase <- z.design[,2:ncol(z.design)]


colnames(z.design) <- c("Luminial A","Luminal B",
                        "Luminal B HER2Neg ",
                        "HER2 Enriched",
                        "Triple Negative")


if(i1<=177){
  data1 <- fread("./data/iCOGS_euro_v10_10232017.csv",header=T)
  data1 <- as.data.frame(data1)
  y.pheno.mis1 <- cbind(data1$Behaviour1,data1$ER_status1,data1$PR_status1,data1$HER2_status1,data1$Grade1)
  colnames(y.pheno.mis1) = c("Behavior","ER","PR","HER2","Grade")
  z.standard <- GenerateZstandard(y.pheno.mis1)
  # Grade1.fake <- data1$Grade1
  # Grade1.fake[data1$Grade1==2|data1$Grade1==3] <- 1
  # Grade1.fake[data1$Grade1==1] <- 0
  #y.pheno.mis1 <- cbind(data1$Behaviour1,data1$PR_status1,data1$ER_status1,data1$HER2_status1,Grade1.fake)
  # y.pheno.mis1 <- cbind(data1$Behaviour1,data1$PR_status1,data1$ER_status1,data1$HER2_status1)
  
  x.test.all.mis1 <- data1[,c(27:203)]
  #x.test.all.mis1 <- 2-x.test.all.mis1
  ###pc1-10 and age
  x.covar.mis1 <- data1[,c(5:14)]
  
  idx.control <- which(y.pheno.mis1[,1]==0)
  maf <- sum(x.test.all.mis1[idx.control,i1])/(2*length(idx.control))
  if(maf>=0.5){
    x.test.all.mis1[,i1] < 2 - x.test.all.mis1[,i1]
  }
  idx <- which(y.pheno.mis1[,l]!=888|
                 y.pheno.mis1[,1]==0)
  y.pheno.mis1 = y.pheno.mis1[idx,]
  x.covar.mis1 = x.covar.mis1[idx,]
  x.test.all.mis1 = x.test.all.mis1[idx,]
  x.all.mis1 <- as.matrix(cbind(x.test.all.mis1[,i1],x.covar.mis1))
  colnames(x.all.mis1)[1] <- "gene"
subtype.icog <- y.pheno.mis1[,l]
subtype.icog[is.na(subtype.icog)] <- "control"
subtype.icog[subtype.icog==0] <- "negative"
subtype.icog[subtype.icog==1] <- "positive"  
subtype.icog = as.factor(subtype.icog)
  library(nnet)
snpvalue <-  x.test.all.mis1[,1] 
model1 <- multinom(subtype.icog~ snpvalue+
                       as.matrix(x.covar.mis1), maxit= 500)
  coef.1 <- coef(model1)
  covar.1 <- vcov(model1)
  jdx <- grep("snpvalue",colnames(covar.1))
  
  
  log.odds.icog <- coef(model1)[,2]
  
  sigma.log.odds.icog <- covar.1[jdx,jdx]
  
  
  
  
  Heter.result.Icog = EMmvpolySelfDesign(y.pheno.mis1,x.self.design = x.all.mis1[,1,drop=F],z.design=z.design,baselineonly = NULL,additive = x.all.mis1[,2:ncol(x.all.mis1)],pairwise.interaction = NULL,saturated = NULL,missingTumorIndicator = 888)
  z.standard <- Heter.result.Icog[[12]]
  z.additive.design <- as.matrix(cbind(1,z.standard))
  M <- nrow(z.standard)
  number.of.tumor <- ncol(z.standard)
  log.odds.icog <- Heter.result.Icog[[1]][(M+1):(M+1+number.of.tumor)]
  nparm <- length(Heter.result.Icog[[1]])  
  sigma.log.odds.icog <- Heter.result.Icog[[2]][(M+1):(M+1+number.of.tumor),(M+1):(M+1+number.of.tumor)]
  
  
  
  names1 = colnames(data1)[27:206]
  
  
  
  #analysis for Onco Array
  #data2 <- read.csv("./V10/Onco_euro_v10_05242017.csv",header=T)
  data2 <- fread("./data/Onco_euro_v10_10232017.csv",header=T)
  data2 <- as.data.frame(data2)
  y.pheno.mis2 <- cbind(data2$Behaviour1,data2$PR_status1,data2$ER_status1,data2$HER2_status1,data2$Grade1)
  #y.pheno.mis2 <- cbind(data2$Behaviour1,data2$PR_status1,data2$ER_status1,data2$HER2_status1)
  colnames(y.pheno.mis2) = c("Behaviour","PR",
                             "ER","HER2","Grade")
  
  x.test.all.mis2 <- data2[,c(27:203)]
  x.test.all.mis2 <- x.test.all.mis2
  x.covar.mis2 <- data2[,c(5:14)]
  
  if(maf>=0.5){
    x.test.all.mis2[,i1] <- 2-x.test.all.mis2[,i1]
  }
  subtypes.onco <- GenerateIntrinsicmis(y.pheno.mis2[,2],
                                        y.pheno.mis2[,3],
                                        y.pheno.mis2[,4],
                                        y.pheno.mis2[,5])
  idx <- which(subtypes.onco!="mis")
  y.pheno.mis2 = y.pheno.mis2[idx,]
  x.covar.mis2 = x.covar.mis2[idx,]
  x.test.all.mis2 = x.test.all.mis2[idx,]
  
  
  
  x.all.mis2 <- as.matrix(cbind(x.test.all.mis2[,i1],x.covar.mis2))
  colnames(x.all.mis2)[1] = "gene"
  idx.control <- which(y.pheno.mis2[,1]==0)
  freq = sum(x.test.all.mis2[idx.control,i1])/(2*length(idx.control))
  
  
  Heter.result.Onco = EMmvpolySelfDesign(y.pheno.mis2,x.self.design = x.all.mis2[,1,drop=F],z.design = z.design,baselineonly = NULL,additive = x.all.mis2[,2:ncol(x.all.mis2)],pairwise.interaction = NULL,saturated = NULL,missingTumorIndicator = 888)
  z.standard <- Heter.result.Onco[[12]]
  M <- nrow(z.standard)
  number.of.tumor <- ncol(z.standard)
  log.odds.onco <- Heter.result.Onco[[1]][(M+1):(M+1+number.of.tumor)]
  nparm <- length(Heter.result.Onco[[1]])
  sigma.log.odds.onco <- Heter.result.Onco[[2]][(M+1):(M+1+number.of.tumor),(M+1):(M+1+number.of.tumor)]
  
  
  meta.result <- LogoddsMetaAnalysis(log.odds.icog,
                                     sigma.log.odds.icog,
                                     log.odds.onco,
                                     sigma.log.odds.onco)
  
  second.stage.logodds.meta <- meta.result[[1]]
  second.stage.sigma.meta <- meta.result[[2]]
  
  
  heter.result <- list(second.stage.logodds.meta,second.stage.sigma.meta,
                       freq)
  save(heter.result,file=paste0("./known_SNPs/known_SNPs_analysis_G_revised/intrinsic_subtypes_pc_additive/result/heter_result_complete_",i1,".Rdata"))
  
  
  
}else{
  data2 <- fread("./data/Onco_euro_v10_10232017.csv",header=T)
  data2 <- as.data.frame(data2)
  names2 = colnames(data2)
  y.pheno.mis2 <- cbind(data2$Behaviour1,data2$PR_status1,data2$ER_status1,data2$HER2_status1,data2$Grade1)
  #y.pheno.mis2 <- cbind(data2$Behaviour1,data2$PR_status1,data2$ER_status1,data2$HER2_status1)
  colnames(y.pheno.mis2) = c("Behaviour","PR",
                             "ER","HER2","Grade")
  idxi1 = which(names2=="rs554219")
  
  x.test.all.mis2 <- data2
  x.covar.mis2 <- data2[,c(5:14)]
  subtypes.onco <- GenerateIntrinsicmis(y.pheno.mis2[,2],
                                        y.pheno.mis2[,3],
                                        y.pheno.mis2[,4],
                                        y.pheno.mis2[,5])
  idx <- which(subtypes.onco!="mis")
  y.pheno.mis2 = y.pheno.mis2[idx,]
  x.covar.mis2 = x.covar.mis2[idx,]
  x.test.all.mis2 = x.test.all.mis2[idx,]
  
  
  x.all.mis2 <- as.matrix(cbind(x.test.all.mis2[,idxi1],x.covar.mis2))
  colnames(x.all.mis2)[1] = "gene"
  
  idx.control <- which(y.pheno.mis2[,1]==0)
  freq = sum(x.test.all.mis2[idx.control,idxi1])/(2*length(idx.control))
  
  
  
  
  
  Heter.result.Onco = EMmvpolySelfDesign(y.pheno.mis2,x.self.design = x.all.mis2[,1,drop=F],z.design = z.design,baselineonly = NULL,additive = x.all.mis2[,2:ncol(x.all.mis2)],pairwise.interaction = NULL,saturated = NULL,missingTumorIndicator = 888)
  z.standard <- Heter.result.Onco[[12]]
  M <- nrow(z.standard)
  number.of.tumor <- ncol(z.standard)
  log.odds.onco <- Heter.result.Onco[[1]][(M+1):(M+1+number.of.tumor)]
  nparm = length(Heter.result.Onco[[1]])
  sigma.log.odds.onco <- Heter.result.Onco[[2]][(M+1):(M+1+number.of.tumor),(M+1):(M+1+number.of.tumor)]
  
  heter.result <- list(log.odds.onco,sigma.log.odds.onco,
                       freq)
  
  
  save(heter.result,file=paste0("./known_SNPs/known_SNPs_analysis_G_revised/intrinsic_subtypes_pc_additive/result/heter_result_complete_",i1,".Rdata"))
  
  
  
  
  
}










