rm(list=ls())
args = commandArgs(trailingOnly = T)
i1 = as.numeric(args[[1]])
i2 = as.numeric(args[[2]])



setwd("/spin1/users/zhangh24/breast_cancer_data_analysis/")
#library(readr)
library(devtools)
library(CompQuadForm)
library(bc2)
library(data.table)
library(bigmemory)
load("/spin1/users/zhangh24/breast_cancer_data_analysis/discovery_SNP/conditional_analysis/result/support.matrix.Rdata")
z.standard <- support.matrix[[1]]
z.additive.design <- support.matrix[[2]]
M <- support.matrix[[3]]
number.of.tumor <- support.matrix[[4]]
z.design.score.baseline <- support.matrix[[5]]
z.design.score.casecase <- support.matrix[[6]]
z.design.score.baseline.ER <- support.matrix[[7]]
z.design.score.casecase.ER <- support.matrix[[8]]

data1 <- fread("./data/iCOGS_euro_v10_10232017.csv",header=T)
data1 <- as.data.frame(data1)
age1 <- as.vector(data1[,204])
idx.complete1 <- which(age1!=888)
age1 <- age1[idx.complete1]

if(i1!=19){
  y.pheno.mis1 <- cbind(data1$Behaviour1,data1$ER_status1,data1$PR_status1,data1$HER2_status1,data1$Grade1)
  colnames(y.pheno.mis1) = c("Behavior","ER","PR","HER2","Grade")
  
  x.covar.mis1 <- data1[,c(5:14)]
  
  y.pheno.mis1 <- y.pheno.mis1[idx.complete1,]
  x.covar.mis1 <- x.covar.mis1[idx.complete1,]
  x.covar.mis1 <- cbind(x.covar.mis1,age1)
  data1.com <- data1[idx.complete1,]
  
  data2 <- fread("./data/Onco_euro_v10_10232017.csv",header=T)
  data2 <- as.data.frame(data2)
  age2 <- data2[,204]
  idx.complete2 <- which(age2!=888)
  age2 <- age2[idx.complete2]
  y.pheno.mis2 <- cbind(data2$Behaviour1,data2$ER_status1,data2$PR_status1,data2$HER2_status1,data2$Grade1)
  
  #y.pheno.mis2 <- cbind(data2$Behaviour1,data2$PR_status1,data2$ER_status1,data2$HER2_status1)
  colnames(y.pheno.mis2) = c("Behaviour","ER",
                             "PR","HER2","Grade")
  
  y.pheno.mis2 <- y.pheno.mis2[idx.complete2,]
  
  x.covar.mis2 <- data2[,c(5:14)]
  x.covar.mis2 <- x.covar.mis2[idx.complete2,]
  x.covar.mis2 <- cbind(x.covar.mis2,age2)
  data2.com <- data2[idx.complete2,]
  discovery.snp.icog <- as.data.frame(fread("/spin1/users/zhangh24/breast_cancer_data_analysis/data/discovery_icog_data.csv",header=T))
  
  discovery.snp.icog.complete <- discovery.snp.icog[idx.complete1,]
  
  discovery.snp.onco <- as.data.frame(fread("/spin1/users/zhangh24/breast_cancer_data_analysis/data/discovery_onco_data.csv"))
  discovery.snp.onco.complete <- discovery.snp.onco[idx.complete2,]
  
  
  
  
  all.countries <- unique(c(data1$StudyCountry,data2$StudyCountry))
  n.coun <- length(all.countries)
  
  
  idx.icog <- which(data1.com$StudyCountry!=all.countries[i2])
  y.pheno.mis1.sub <- y.pheno.mis1[idx.icog,]
  x.covar.mis1.sub <- x.covar.mis1[idx.icog,]
  discovery.snp.icog.sub <- discovery.snp.icog.complete[idx.icog,i1]
  x.all.mis1.sub <- cbind(discovery.snp.icog.sub,x.covar.mis1.sub)
  # x.all.mis1 <- cbind(discovery.snp.icog.complete[,i1],x.covar.mis1)
  Heter.result.Icog = EMmvpoly(y.pheno.mis1.sub,baselineonly = NULL,additive = x.all.mis1.sub,pairwise.interaction = NULL,saturated = NULL,missingTumorIndicator = 888)
  z.standard <- Heter.result.Icog[[12]]
  M <- nrow(z.standard)
  log.odds.icog <- Heter.result.Icog[[1]][(M+1):(M+1+number.of.tumor)]
  
  sigma.log.odds.icog <- Heter.result.Icog[[2]][(M+1):(M+1+number.of.tumor),(M+1):(M+1+number.of.tumor)]
  
  
  
  
  
  
  
  
  
  
  idx.onco <- which(data2.com$StudyCountry!=all.countries[i2])
  print(length(idx.onco))
  y.pheno.mis2.sub <- y.pheno.mis2[idx.onco,]
  x.covar.mis2.sub <- x.covar.mis2[idx.onco,]
  discovery.snp.onco.sub <- discovery.snp.onco.complete[idx.onco,i1]
  x.all.mis2.sub <- cbind(discovery.snp.onco.sub,x.covar.mis2.sub)
  Heter.result.Onco = EMmvpoly(y.pheno.mis2.sub,baselineonly = NULL,additive = x.all.mis2.sub,pairwise.interaction = NULL,saturated = NULL,missingTumorIndicator = 888)
  z.standard <- Heter.result.Onco[[12]]
  M <- nrow(z.standard)
  number.of.tumor <- ncol(z.standard)
  log.odds.onco <- Heter.result.Onco[[1]][(M+1):(M+1+number.of.tumor)]
  nparm <- length(Heter.result.Onco)
  sigma.log.odds.onco <- Heter.result.Onco[[2]][(M+1):(M+1+number.of.tumor),(M+1):(M+1+number.of.tumor)] 
  
  
  meta.result <- LogoddsMetaAnalysis(log.odds.icog,
                                     sigma.log.odds.icog,
                                     log.odds.onco,
                                     sigma.log.odds.onco)
  
  second.stage.logodds.meta <- meta.result[[1]]
  second.stage.sigma.meta <- meta.result[[2]]
  
  
  
  test.result.second.wald <- DisplaySecondStageTestResult(second.stage.logodds.meta,second.stage.sigma.meta)
  
  result <- list(second.stage.logodds.meta,second.stage.sigma.meta,
                 test.result.second.wald)
  
  save(result,file=paste0("./discovery_SNP/sensitivity_analysis/sensitivity_analysis",i1,"_",i2,".Rdata"))
  
  
  
  
  
  
  
  
  
  
  
}else{
  y.pheno.mis1 <- cbind(data1$Behaviour1,data1$ER_status1,data1$PR_status1,data1$HER2_status1,data1$Grade1)
  colnames(y.pheno.mis1) = c("Behavior","ER","PR","HER2","Grade")
  
  x.covar.mis1 <- data1[,c(5:14)]
  
  y.pheno.mis1 <- y.pheno.mis1[idx.complete1,]
  x.covar.mis1 <- x.covar.mis1[idx.complete1,]
  x.covar.mis1 <- cbind(x.covar.mis1,age1)
  data1.com <- data1[idx.complete1,]
  
  data2 <- fread("./data/Onco_euro_v10_10232017.csv",header=T)
  data2 <- as.data.frame(data2)
  age2 <- data2[,204]
  idx.complete2 <- which(age2!=888)
  age2 <- age2[idx.complete2]
  y.pheno.mis2 <- cbind(data2$Behaviour1,data2$ER_status1,data2$PR_status1,data2$HER2_status1,data2$Grade1)
  
  #y.pheno.mis2 <- cbind(data2$Behaviour1,data2$PR_status1,data2$ER_status1,data2$HER2_status1)
  colnames(y.pheno.mis2) = c("Behaviour","ER",
                             "PR","HER2","Grade")
  
  y.pheno.mis2 <- y.pheno.mis2[idx.complete2,]
  
  x.covar.mis2 <- data2[,c(5:14)]
  x.covar.mis2 <- x.covar.mis2[idx.complete2,]
  x.covar.mis2 <- cbind(x.covar.mis2,age2)
  data2.com <- data2[idx.complete2,]
  discovery.snp.icog <- as.data.frame(fread("/spin1/users/zhangh24/breast_cancer_data_analysis/data/discovery_icog_data.csv",header=T))
  
  discovery.snp.icog.complete <- discovery.snp.icog[idx.complete1,]
  
  discovery.snp.onco <- as.data.frame(fread("/spin1/users/zhangh24/breast_cancer_data_analysis/data/discovery_onco_data.csv"))
  discovery.snp.onco.complete <- discovery.snp.onco[idx.complete2,]
  
  
  
  
  all.countries <- unique(c(data1$StudyCountry,data2$StudyCountry))
  n.coun <- length(all.countries)
  
  
  idx.icog <- which(data1.com$StudyCountry!=all.countries[i2])
  y.pheno.mis1.sub <- y.pheno.mis1[idx.icog,]
  x.covar.mis1.sub <- x.covar.mis1[idx.icog,]
  discovery.snp.icog.sub <- discovery.snp.icog.complete[idx.icog,i1]
  x.all.mis1.sub <- cbind(discovery.snp.icog.sub,x.covar.mis1.sub)
  # x.all.mis1 <- cbind(discovery.snp.icog.complete[,i1],x.covar.mis1)
  
  
  
  
  idx.onco <- which(data2.com$StudyCountry!=all.countries[i2])
  y.pheno.mis2.sub <- y.pheno.mis2[idx.onco,]
  x.covar.mis2.sub <- x.covar.mis2[idx.onco,]
  discovery.snp.onco.sub <- discovery.snp.onco.complete[idx.onco,i1]
  x.all.mis2.sub <- cbind(discovery.snp.onco.sub,x.covar.mis2.sub)
  Heter.result.Onco = EMmvpoly(y.pheno.mis2.sub,baselineonly = NULL,additive = x.all.mis2.sub,pairwise.interaction = NULL,saturated = NULL,missingTumorIndicator = 888)
  z.standard <- Heter.result.Onco[[12]]
  M <- nrow(z.standard)
  number.of.tumor <- ncol(z.standard)
  log.odds.onco <- Heter.result.Onco[[1]][(M+1):(M+1+number.of.tumor)]
  nparm <- length(Heter.result.Onco)
  sigma.log.odds.onco <- Heter.result.Onco[[2]][(M+1):(M+1+number.of.tumor),(M+1):(M+1+number.of.tumor)] 
  
  
  second.stage.logodds.meta <- log.odds.onco
  second.stage.sigma.meta <-  sigma.log.odds.onco
  
  
  
  test.result.second.wald <- DisplaySecondStageTestResult(second.stage.logodds.meta,second.stage.sigma.meta)
  
  result <- list(second.stage.logodds.meta,second.stage.sigma.meta,
                 test.result.second.wald)
  
  save(result,file=paste0("./discovery_SNP/sensitivity_analysis/sensitivity_analysis",i1,"_",i2,".Rdata"))
  
  
  
  
  
  
  
  
  
  
  
}
