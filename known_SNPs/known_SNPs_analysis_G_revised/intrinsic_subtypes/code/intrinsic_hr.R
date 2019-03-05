#install_github("andrewhaoyu/bc2", args = c('--library="/home/zhangh24/R/x86_64-pc-linux-gnu-library/3.4"'))
#install_github("andrewhaoyu/bcutility",args = c('--library="/home/zhangh24/R/x86_64-pc-linux-gnu-library/3.4"'))
#install_github("andrewhaoyu/bc2",branch="version 0.0.2")
###1 represent Icog
###2 represent Onco
###load_all("/Users/zhangh24/GoogleDrive/bc2")

rm(list=ls())
#commandarg <- commandArgs(trailingOnly=F)
#myarg <- commandarg[length(commandarg)]
#myarg <- sub("-","",myarg)
#i1 <- as.numeric(myarg)
args = commandArgs(trailingOnly = T)
i1 = as.numeric(args[[1]])

setwd("/data/zhangh24/breast_cancer_data_analysis")
library(readr)
library(devtools)
#library(bc2,lib.loc ='/Users/zhangh24/Library/R/3.4/library')
library(bc2)
library(bcutility)
z.design <- matrix(c(
  c(0,1,1,1,0,0,0,0,1,1,1,0,0,0,0,0,1,1,0,0,0,0,0),
  c(0,0,0,0,0,1,1,0,0,0,0,0,1,1,1,0,0,0,0,0,1,1,1),
  c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0),
  c(0,0,0,0,1,0,0,0,0,0,0,1,0,0,0,0,0,0,0,1,0,0,0),
  c(1,0,0,0,0,0,0,1,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0)
),ncol=5)
z.design[,1] <- z.design[,1]+z.design[,2]+z.design[,3]
z.design <- z.design[,c(1,4,5)]
rowSums(z.design)
colnames(z.design) <- c("HR+",
                        "HER2 Enriched",
                        "Triple Negative")





if(i1<=177){
  ##analysis for Icog
  data1 <- fread("./data/iCOGS_euro_v10_10232017.csv",header=T)
  y.pheno.mis1 <- cbind(data1$Behaviour1,data1$PR_status1,data1$ER_status1,data1$HER2_status1,data1$Grade1)
  # y.pheno.mis1 <- cbind(data1$Behaviour1,data1$PR_status1,
  #data1$ER_status1,data1$HER2_status1)
  #y.pheno.mis1 <- cbind(data1$Behaviour1,data1$PR_status1,data1$ER_status1,data1$HER2_status1)
  colnames(y.pheno.mis1) <- c("Behaviour","PR","ER","HER2","Grade")
  # colnames(y.pheno.mis1) <- c("Behaviour","PR","ER","HER2")
  # Grade1.fake <- data1$Grade1
  # Grade1.fake[data1$Grade1==2|data1$Grade1==3] <- 1
  # Grade1.fake[data1$Grade1==1] <- 0
  #y.pheno.mis1 <- cbind(data1$Behaviour1,data1$PR_status1,data1$ER_status1,data1$HER2_status1,Grade1.fake)
  # y.pheno.mis1 <- cbind(data1$Behaviour1,data1$PR_status1,data1$ER_status1,data1$HER2_status1)
  
  x.test.all.mis1 <- data1[,c(27:206)]
  
  x.covar.mis1 <- data1[,5:14]
  x.all.mis1 <- as.matrix(cbind(x.test.all.mis1[,i1],x.covar.mis1))
  colnames(x.all.mis1)[1] <- "gene"
  
  Heter.result.G.Icog = EMmvpolySelfDesign(y.pheno.mis1, x= x.all.mis1 ,z.design = z.design,missingTumorIndicator = 888,z.all = NULL)
  
  M <- 23
  number.of.tumor <- 4
  log.odds.icog <- Heter.result.G.Icog[[1]][(M+1):(M+1+number.of.tumor)]
  sigma.icog <- Heter.result.G.Icog[[2]]
  sigma.log.odds.icog <- sigma.icog[(M+1):(M+1+number.of.tumor),(M+1):(M+1+number.of.tumor)]
  loglikelihood.icog <- Heter.result.G.Icog[[8]]
  
  
  
  #analysis for Onco Array
  #data2 <- read.csv("./V10/Onco_euro_v10_05242017.csv",header=T)
  data2 <- fread("./data/Onco_euro_v10_10232017.csv",header=T)
  names1 = colnames(data1)[27:206]
  rm(data1)
  names2 = colnames(data2)[27:212]
  
  idxi1 = which(names2==names1[i1])
  y.pheno.mis2 <- cbind(data2$Behaviour1,data2$PR_status1,data2$ER_status1,data2$HER2_status1,data2$Grade1)
  #y.pheno.mis2 <- cbind(data2$Behaviour1,data2$PR_status1,data2$ER_status1,data2$HER2_status1)
  
  x.test.all.mis2 <- data2[,c(27:212)]
  x.covar.mis2 <- data2[,5:14]
  x.all.mis2 <- as.matrix(cbind(x.test.all.mis2[,idxi1],x.covar.mis2))
  
  
  
  Heter.result.G.Onco = EMmvpolySelfDesign(y.pheno.mis2,x= x.all.mis2 ,z.design = z.design,missingTumorIndicator = 888,z.all = NULL)
  
  number.of.tumor <- 4
  log.odds.onco <- Heter.result.G.Onco[[1]][(M+1):(M+1+number.of.tumor)]
  
  sigma.onco <- Heter.result.G.Onco[[2]]
  sigma.log.odds.onco <- sigma.onco[(M+1):(M+1+number.of.tumor),(M+1):(M+1+number.of.tumor)]
  loglikelihood.onco <- Heter.result.G.Onco[[8]]
  
  
  meta.result <- LogoddsMetaAnalysis(log.odds.icog,
                                     sigma.log.odds.icog,
                                     log.odds.onco,
                                     sigma.log.odds.onco)
  
  logodds <- meta.result[[1]]
  sigma <- meta.result[[2]]
  loglikelihood.meta <- loglikelihood.icog+loglikelihood.onco
  AIC <- 2*length(Heter.result.G.Onco[[1]])-2*loglikelihood.meta
  test.result <- DisplayTestResult(logodds,sigma)
  
  heter.result <- list(test.result = test.result,
                       logodds = logodds,
                       sigma= sigma,
                       loglikelihood.meta=loglikelihood.meta,
                       AIC=AIC)
  save(heter.result,file=paste0("./known_SNPs/known_SNPs_analysis_G_revised/intrinsic_subtypes/result/heter_result_hr",i1,".Rdata"))
  
  
  
}else{
  M <- 23
  data2 <- fread("./data/Onco_euro_v10_10232017.csv",header=T)
  #data2 <- read.csv("/data/zhangh20/breast_cancer/known_SNPs_anlysis/Onco_euro_v10_rs554219.csv",1)
  
  # names1 = colnames(data1)[27:206]
  #rm(data1)
  names2 = colnames(data2)
  
  idxi1 = which(names2=="rs554219")
  y.pheno.mis2 <- cbind(data2$Behaviour1,data2$PR_status1,data2$ER_status1,data2$HER2_status1,data2$Grade1)
  x.test.all.mis2 <- data2
  x.covar.mis2 <- data2[,5:14]
  x.all.mis2 <- as.matrix(cbind(x.test.all.mis2[,idxi1],x.covar.mis2))
  
  
  
  Heter.result.G.Onco = EMmvpolySelfDesign(y.pheno.mis2,x= x.all.mis2 ,z.design = z.design,missingTumorIndicator = 888,z.all = NULL)
  
  number.of.tumor <- 4
  log.odds.onco <- Heter.result.G.Onco[[1]][(M+1):(M+1+number.of.tumor)]
  
  sigma.onco <- Heter.result.G.Onco[[2]]
  sigma.log.odds.onco <- sigma.onco[(M+1):(M+1+number.of.tumor),(M+1):(M+1+number.of.tumor)]
  loglikelihood.onco <- Heter.result.G.Onco[[8]]
  
  
  
  logodds <- log.odds.onco
  sigma <- sigma.log.odds.onco
  
  loglikelihood.meta <- loglikelihood.onco
  AIC <- 2*length(Heter.result.G.Onco[[1]])-2*loglikelihood.meta
  test.result <- DisplayTestResult(logodds,sigma)
  
  heter.result <- list(test.result = test.result,
                       logodds = logodds,
                       sigma= sigma,
                       loglikelihood.meta=loglikelihood.meta,
                       AIC=AIC)
  save(heter.result,file=paste0("./known_SNPs/known_SNPs_analysis_G_revised/intrinsic_subtypes/result/heter_result_hr",i1,".Rdata"))
  
  
  
  
  
}



