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
library(data.table)
#library(bc2,lib.loc ='/Users/zhangh24/Library/R/3.4/library')
#install_github("andrewhaoyu/bc2")
#install_github("andrewhaoyu/bcutility")
#library(bcutility,lib.loc = "/spin1/home/linux/zhangh24/R/x86_64-pc-linux-gnu-library/3.6/")
#library(bc2, lib.loc ="/home/zhangh24/R/x86_64-pc-linux-gnu-library/3.6/")
library(bc2)
library(bcutility)
# z.design <- matrix(c(
#   c(0,1,1,1,0,0,0,0,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0),
#   c(0,0,0,0,0,1,1,0,0,0,0,0,1,1,1,0,0,0,0,0,1,1,1),
#   c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,0,0,0,0),
#   c(0,0,0,0,1,0,0,0,0,0,0,1,0,0,0,0,0,0,0,1,0,0,0),
#   c(1,0,0,0,0,0,0,1,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0)
# ),ncol=5)
# colnames(z.design) <- c("Luminial A-like","Luminal B-like",
#                     "Luminal B HER2 negative-like",
#                     "HER2 Enriched-like",
#                     "Triple Negative")




if(i1<=177){
  ##analysis for Icog
  data1 <- as.data.frame(fread("./data/iCOGS_euro_v10_10232017.csv",header=T))
  #prepare the phenotypes data for iCOGs 
  y.pheno.mis1 <- cbind(data1$Behaviour1,data1$ER_status1,data1$PR_status1,data1$HER2_status1,data1$Grade1)
  #ER, PR, HER2 is binary with negative as 0 and positive as 1
  #Grade is ordinal with 1, 2, 3
  #controls don't have tumor characteristics data (all NA)
  #cases with missing tumor characteristics marked as 888
  colnames(y.pheno.mis1) <- c("Behaviour","ER","PR","HER2","Grade")
  #generate the z standard matrix
  #z standard matrix is link to link subtypes with tumor characteristics
  z.standard <- GenerateZstandard(y.pheno.mis1)
  #each row of z.standard represent a subtype
  #each column represent a tumor marker
  #e.g. first row is ER-PR-HER2-Grade 1
  #total number of subtypes
  #subtypes with less than 10 cases are automatically removed
  M <- nrow(z.standard)
  #construct the z design matrix for intrinsic subtypes
  #intrinsic subtypes are defined as follows
  #Luminal-A like: ER or PR +, HER2-, grade is 1 or 2
  #Luminal-B like: ER or PR +, HER2+
  #Luminal B HER2 negative-like: ER or PR+, HER2-, grade 3
  #HER2 enriched-like: both ER and PR-, HER2+
  #Triple negative: ER, PR, HER2-
  #Prepare z design matrix
  z.design <- matrix(0,M,5)
  #define Luminal-A like
  idx.1 <- which((z.standard[,1]==1|z.standard[,2]==1)
                 &z.standard[,3]==0
                 &(z.standard[,4]==1|z.standard[,4]==2))
  z.design[idx.1,1] <- 1
  #define Luminal-B like
  idx.2 <- which((z.standard[,1]==1|z.standard[,2]==1)
                 &z.standard[,3]==1)
  z.design[idx.2,2] <- 1
  #for Luminal B HER2 negative-like
  idx.3 <- which((z.standard[,1]==1|z.standard[,2]==1)
                 &z.standard[,3]==0
                 &z.standard[,4]==3)
  z.design[idx.3,3] <- 1
  #for HER2 enriched-like
  idx.4 <- which(z.standard[,1]==0&z.standard[,2]==0
                 &z.standard[,3]==1)
  z.design[idx.4,4] <- 1
  #for Triple negative
  idx.5 <- which(z.standard[,1]==0&z.standard[,2]==0
                 &z.standard[,3]==0)
  z.design[idx.5,5] <- 1
  
  #genotype data for 178 known SNPs (178th SNP only exit on Oncoarray)
  x.test.all.mis1 <- data1[,c(27:203)]
  #prepare covariates table: PC1-10
  #we only adjusted PC1-10 in known SNPs analyses
  #in genome-wide analyses, we will need to adjust age and PC1-10
  x.covar.mis1 <- data1[,5:14]
  
  #################################
  #this section is not important
  #it's mainly used to convert the effect-size to minor allele
  #no need for this in analyses
  idx.control <- which(y.pheno.mis1[,1]==0)
  maf <- sum(x.test.all.mis1[idx.control,i1])/(2*length(idx.control))
  if(maf>=0.5){
    x.test.all.mis1[,i1] < 2 - x.test.all.mis1[,i1]
  }
  #################################
  
  
  x.all.mis1 <- as.matrix(cbind(x.test.all.mis1[,i1],x.covar.mis1))
  colnames(x.all.mis1)[1] <- "gene"
  
  Heter.result.G.Icog = EMmvpolySelfDesign(y.pheno.mis1, 
                                           x= x.all.mis1 ,
                                           z.design = z.design,
                                           missingTumorIndicator = 888,
                                           z.all = NULL)
  
  M <- 23
  number.of.tumor <- 4
  log.odds.icog <- Heter.result.G.Icog[[1]][(M+1):(M+1+number.of.tumor)]
  sigma.icog <- Heter.result.G.Icog[[2]]
  sigma.log.odds.icog <- sigma.icog[(M+1):(M+1+number.of.tumor),(M+1):(M+1+number.of.tumor)]
  loglikelihood.icog <- Heter.result.G.Icog[[8]]
  
  
  
  #analysis for Onco Array
  #data2 <- read.csv("./V10/Onco_euro_v10_05242017.csv",header=T)
  data2 <- as.data.frame(fread("./data/Onco_euro_v10_10232017.csv",header=T))
  names1 = colnames(data1)[c(27:203)]
  rm(data1)
  names2 = colnames(data2)[27:205]
  
  idxi1 = which(names2==names1[i1])
  y.pheno.mis2 <- cbind(data2$Behaviour1,data2$ER_status1,data2$PR_status1,data2$HER2_status1,data2$Grade1)
  #y.pheno.mis2 <- cbind(data2$Behaviour1,data2$PR_status1,data2$ER_status1,data2$HER2_status1)
  
  x.test.all.mis2 <- data2[,c(27:205)]
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
  save(heter.result,file=paste0("./known_SNPs/known_SNPs_analysis_G_revised/intrinsic_subtypes/result/heter_result_",i1,".Rdata"))
  
  
  
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
  save(heter.result,file=paste0("./known_SNPs/known_SNPs_analysis_G_revised/intrinsic_subtypes/result/heter_result_",i1,".Rdata"))
  
  
  
  
  
}



