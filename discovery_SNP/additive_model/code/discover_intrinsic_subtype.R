# Goal: estimate log odds ratio and var of discovery SNPs
#install_github("andrewhaoyu/bc2",ref='development', args = c('--library="/home/zhangh24/R/x86_64-pc-linux-gnu-library/3.4"'))
#install_github("andrewhaoyu/bc2", args = c('--library="/home/zhangh24/R/x86_64-pc-linux-gnu-library/3.4"'))
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

setwd("/spin1/users/zhangh24/breast_cancer_data_analysis/")

library(readr)
library(devtools)
library(CompQuadForm)
library(bc2)
library(data.table)
discovery.snp.icog <- as.data.frame(fread("/spin1/users/zhangh24/breast_cancer_data_analysis/discovery_SNP/result/discovery_icog_data.csv",header=T))
colnames(discovery.snp.icog)
#onco.julie <- fread("/spin1/users/zhangh24/breast_cancer_data_analysis/data/Julie_snp_onco.csv")
#onco.julie <- onco.julie[,-1]
x.test.all.mis1 <- discovery.snp.icog

#sum(x.test.all.mis2[,11])/(2*nrow(x.test.all.mis2))

discovery_snp <- read.csv("/spin1/users/zhangh24/breast_cancer_data_analysis/data/discovery_snp_summary_new.csv",header=T)

#x.test.all.mis1 <- as.data.frame(cbind(icog.julie,discovery.snp.icog))
#x.test.all.mis2 <- as.data.frame(cbind(onco.julie,discovery.snp.onco))

z.design <- matrix(c(
  c(0,1,1,1,0,0,0,0,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0),
  c(0,0,0,0,0,1,1,0,0,0,0,0,1,1,1,0,0,0,0,0,1,1,1),
  c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,0,0,0,0),
  c(0,0,0,0,1,0,0,0,0,0,0,1,0,0,0,0,0,0,0,1,0,0,0),
  c(1,0,0,0,0,0,0,1,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0)
),ncol=5)
rowSums(z.design)



colnames(z.design) <- c("Luminial A","Luminal B",
                        "Luminal B HER2Neg",
                        "HER2 Enriched",
                        "Triple Negative")




#load("./whole_genome_age/ONCO/ERPRHER2GRADE_fixed_baseline/result/delta0.onco.Rdata")
#load("./whole_genome_age/ONCO/ERPRHER2GRADE_fixed_baseline/result/z.standard.Rdata")

#z.design.support <- cbind(1,z.standard[,1])
#z.design.test <- z.standard[,2:4]
#if(i1<=28){
  ##analysis for Icog
  data1 <- fread("./data/iCOGS_euro_v10_10232017.csv",header=T)
  data1 <- as.data.frame(data1)
  y.pheno.mis1 <- cbind(data1$Behaviour1,data1$ER_status1,data1$PR_status1,data1$HER2_status1,data1$Grade1)
  colnames(y.pheno.mis1) = c("Behavior","ER","PR","HER2","Grade")
  x.covar.mis1 <- data1[,c(5:14,204)]
  if(discovery_snp$exp_freq_a1[i1]>0.5){
    x.test.all.mis1[,i1] = 2-x.test.all.mis1[,i1]
  }
  
  age <- data1[,204]
  snpvalue <- x.test.all.mis1[,i1]
  #x.all.mis1 <- as.matrix(cbind(x.test.all.mis1[,i1],x.covar.mis1))
  age <- data1[,204]
  idx.complete <- which(age!=888)
  snpvalue <- snpvalue[idx.complete]
  y.pheno.mis1 <- y.pheno.mis1[idx.complete,]
  #x.all.mis1 <- x.all.mis1[idx.complete,]
  x.covar.mis1 <- x.covar.mis1[idx.complete,]
  #colnames(x.all.mis1)[1] <- "gene"
  
  
  Heter.result.Icog = EMmvpolySelfDesign(y.pheno.mis1,x.self.design = snpvalue,z.design=z.design,baselineonly = NULL,additive = x.covar.mis1,pairwise.interaction = NULL,saturated = NULL,missingTumorIndicator = 888)
  z.standard <- Heter.result.Icog[[12]]
  M <- nrow(z.standard)
  number.of.tumor <- ncol(z.standard)
  log.odds.icog <- Heter.result.Icog[[1]][(M+1):(M+1+number.of.tumor)]
  nparm <- length(Heter.result.Icog[[1]])  
  sigma.log.odds.icog <- Heter.result.Icog[[2]][(M+1):(M+1+number.of.tumor),(M+1):(M+1+number.of.tumor)]  
  
  
  sum(snpvalue)/(2*length(snpvalue))
  
  
  

  
  #analysis for Onco Array
  #data2 <- read.csv("./V10/Onco_euro_v10_05242017.csv",header=T)
  data2 <- fread("./data/Onco_euro_v10_10232017.csv",header=T)
  data2 <- as.data.frame(data2)
  y.pheno.mis2 <- cbind(data2$Behaviour1,data2$ER_status1,data2$PR_status1,data2$HER2_status1,data2$Grade1)
  #y.pheno.mis2 <- cbind(data2$Behaviour1,data2$PR_status1,data2$ER_status1,data2$HER2_status1)
  colnames(y.pheno.mis2) = c("Behaviour","ER",
                             "PR","HER2","Grade")
  #x.test.all.mis2 <- data2[,c(27:203)]
  discovery.snp.onco <- as.data.frame(fread("/spin1/users/zhangh24/breast_cancer_data_analysis/discovery_SNP/result/discovery_onco_data.csv",header=T))
  
  x.test.all.mis2 <- discovery.snp.onco
  x.covar.mis2 <- data2[,c(5:14,204)]
  if(discovery_snp$exp_freq_a1[i1]>0.5){
    x.test.all.mis2[,i1] = 2-x.test.all.mis2[,i1]
  }
  sum(x.test.all.mis2[,i1])/(2*nrow(x.test.all.mis2))
 # sum(x.test.all.mis1[,i1])/(2*nrow(x.test.all.mis1))
  x.all.mis2 <- as.matrix(cbind(x.test.all.mis2[,i1],x.covar.mis2))
  ages <- data2[,204]
  idx.complete <- which(ages!=888)
  
  y.pheno.mis2 <- y.pheno.mis2[idx.complete,]
  x.all.mis2 <- x.all.mis2[idx.complete,]
  x.covar.mis2 <- x.covar.mis2[idx.complete,]
  colnames(x.all.mis2)[1] = "gene"
  
  snpvalue = x.all.mis2[,1,drop=F]

  
  
  
  
  Heter.result.Onco = EMmvpolySelfDesign(y.pheno.mis2,x.self.design = snpvalue,z.design = z.design,baselineonly = NULL,additive = x.covar.mis2,pairwise.interaction = NULL,saturated = NULL,missingTumorIndicator = 888)
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
  
  
  
  test.result.second.wald <- DisplaySecondStageTestResult(second.stage.logodds.meta,second.stage.sigma.meta)
  
  save( test.result.second.wald,file=paste0("./discovery_SNP/additive_model/result/intrinsic_subtype_",i1,".Rdata"))
  
  
  
#}