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
  gene_value = x.test.all.mis1[,i1,drop = F]
  #Fit the two-stage model
  #Two stage have several z design matrix structure
  #baseline only assume all subtypes have the same effect
  #additive assumes all higher order interactions across to be 0
  #pair-wise interaction allows main effect and second order interactions across tumor markers
  #saturated model allows all interactions
  #we can also self design the z matrix
  #in this setting, we self define the z matrix to make the definition align with intrinsic subtypes
  #we only do this self define z matrix for genetic variant
  #meanwhile, we keep the additive model for all other covariates
  Heter.result.Icog = EMmvpolySelfDesign(y.pheno.mis1,
                                         x.self.design = gene_value,
                                         z.design=z.design,
                                         baselineonly = NULL,
                                         additive = x.covar.mis1,
                                         pairwise.interaction = NULL,
                                         saturated = NULL,missingTumorIndicator = 888)
  n.param <- ncol(z.design)
  #the log-odds ratio for second-stage parameters are saved in the first elelment
  #first M parameters are for intercept
  #We don't make any assumptions regarding the intercept
  #Therefore, the intercept has the same degree of freedom panelty as the M subtypes
  #The next five parameters are for the genetic variants
  #They represent the log-odds ratio for the five intrinsic subtypes
  log.odds.icog <- Heter.result.Icog[[1]][(M+1):(M+n.param)]
  #nparm <- length(Heter.result.Icog[[1]])  
  #variance matrix for the log-odds-ratio are saved in the second component
  sigma.log.odds.icog <- Heter.result.Icog[[2]][(M+1):(M+n.param),(M+1):(M+n.param)]
  names1 = colnames(data1)[27:206]
  
  
  
  #analysis for Onco Array
  #data2 <- read.csv("./V10/Onco_euro_v10_05242017.csv",header=T)
  data2 <- fread("./data/Onco_euro_v10_10232017.csv",header=T)
  data2 <- as.data.frame(data2)
  y.pheno.mis2 <- cbind(data2$Behaviour1,data2$ER_status1,data2$PR_status1,data2$HER2_status1,data2$Grade1)
  #y.pheno.mis2 <- cbind(data2$Behaviour1,data2$PR_status1,data2$ER_status1,data2$HER2_status1)
  colnames(y.pheno.mis2) = c("Behaviour","ER",
                             "PR","HER2","Grade")
  
  x.test.all.mis2 <- data2[,c(27:203)]
  x.test.all.mis2 <- x.test.all.mis2
  x.covar.mis2 <- data2[,c(5:14)]
  gene_value = x.test.all.mis2[,i1,drop = F]
  
  Heter.result.Onco = EMmvpolySelfDesign(y.pheno.mis2,
                                         x.self.design = gene_value,
                                         z.design = z.design,
                                         baselineonly = NULL,
                                         additive = x.covar.mis2,
                                         pairwise.interaction = NULL,
                                         saturated = NULL,missingTumorIndicator = 888)
  
  log.odds.onco <- Heter.result.Onco[[1]][(M+1):(M+n.param)]
  sigma.log.odds.onco <- Heter.result.Onco[[2]][(M+1):(M+n.param),(M+1):(M+n.param)]
  
  meta.result <- LogoddsMetaAnalysis(log.odds.icog,
                                     sigma.log.odds.icog,
                                     log.odds.onco,
                                     sigma.log.odds.onco)
  
  second.stage.logodds.meta <- meta.result[[1]]
  second.stage.sigma.meta <- meta.result[[2]]
  
  #function to calculate z-statistics
  z = second.stage.logodds.meta/sqrt(diag(second.stage.sigma.meta))
  #function to calculate P-value
  p = 2*pnorm(-abs(z), lower.tail = T)
  
  heter.result <- list(second.stage.logodds.meta,second.stage.sigma.meta)
  save(heter.result,file=paste0("./known_SNPs/known_SNPs_analysis_G_revised/intrinsic_subtypes_pc_additive/result/heter_result_origin",i1,".Rdata"))
  
  
  
}else{
  #I didn't update code in this section
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
  beta.onco <- z.trans%*%log.odds.onco
  beta.sigma.onco <- z.trans%*%sigma.log.odds.onco%*%t(z.trans)
  loglikelihood.onco <- Heter.result.Onco[[8]]
  
  
  
  
  standard.result.onco <- glm(y.pheno.mis2[,1]~
                                cbind( x.all.mis2[,1,drop=F],
                                       x.all.mis2[,2:ncol(x.all.mis2)]),
                              family = binomial)
  model.summary.onco <- summary(standard.result.onco)
  log.odds.onco.standard <-  model.summary.onco$coefficient[2,1]
  var.onco.standard <- model.summary.onco$coefficient[2,2]^2
  
  heter.result <- list(log.odds.onco,sigma.log.odds.onco,
                       log.odds.onco.standard,
                       var.onco.standard,
                       freq)
  
  
  save(heter.result,file=paste0("./known_SNPs/known_SNPs_analysis_G_revised/intrinsic_subtypes_pc_additive/result/heter_result_origin",i1,".Rdata"))
  
  
  
  
  
}




