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
library(bc2)




if(i1<=180){
  ##analysis for Icog
  data1 <- read.csv("./data/iCOGS_euro_v10_05242017.csv",header=T)
  y.pheno.mis1 <- cbind(data1$Behaviour1,data1$PR_status1,data1$ER_status1,data1$HER2_status1)
  colnames(y.pheno.mis1) = c("Behavior","PR","ER","HER2")
  # Grade1.fake <- data1$Grade1
  # Grade1.fake[data1$Grade1==2|data1$Grade1==3] <- 1
  # Grade1.fake[data1$Grade1==1] <- 0
  #y.pheno.mis1 <- cbind(data1$Behaviour1,data1$PR_status1,data1$ER_status1,data1$HER2_status1,Grade1.fake)
  # y.pheno.mis1 <- cbind(data1$Behaviour1,data1$PR_status1,data1$ER_status1,data1$HER2_status1)
  
  x.test.all.mis1 <- data1[,c(27:206)]
  
  x.covar.mis1 <- data1[,5:14]
  x.all.mis1 <- as.matrix(cbind(x.test.all.mis1[,i1],x.covar.mis1))
  colnames(x.all.mis1)[1] <- "gene"
  
  Heter.result.Icog = EMmvpoly(y.pheno.mis1,baselineonly = NULL,additive = x.all.mis1,pairwise.interaction = NULL,saturated = NULL,missingTumorIndicator = 888)
  z.standard <- Heter.result.Icog[[12]]
  z.additive.design <- as.matrix(cbind(1,z.standard))
  M <- nrow(z.standard)
  number.of.tumor <- ncol(z.standard)
  log.odds.icog <- Heter.result.Icog[[1]][(M+1):(M+1+number.of.tumor)]
  
  sigma.log.odds.icog <- Heter.result.Icog[[2]][(M+1):(M+1+number.of.tumor),(M+1):(M+1+number.of.tumor)]
  beta.icog <- z.additive.design%*%log.odds.icog
  beta.sigma.icog <- z.additive.design%*%sigma.log.odds.icog%*%t(z.additive.design)
  loglikelihood.icog <- Heter.result.Icog[[8]]
  
  score.test.support.icog <- ScoreTestSupportMixedModel(
    y.pheno.mis1,
    baselineonly = NULL,
    additive = x.all.mis1[,2:3],
    pairwise.interaction = NULL,
    saturated = NULL,
    missingTumorIndicator = 888
  )
  score.test<- ScoreTestMixedModel(y=y.pheno.mis1,
                         x=x.all.mis1[,1],
                         second.stage.structure="additive",
                         score.test.support=score.test.support.icog,
                         missingTumorIndicator=888)
  z.design.score.baseline <- matrix(rep(1,8),ncol=1)
  z.design.score.casecase <- matrix(c(
    c(0,1,0,1,0,1,0,1),
    c(0,0,1,1,0,0,1,1),
    c(0,0,0,0,1,1,1,1)
  ),ncol=3)
  
  score.icog <- score.test[[1]]
  infor.icog <- score.test[[2]]
  
  DisplayFixedScoreTestResult(score.icog,infor.icog)
  
  score.test.support.icog.baseline <- ScoreTestSupport(
    y.pheno.mis1,
    baselineonly = NULL,
    additive = x.all.mis1[,2:11],
    pairwise.interaction = NULL,
    saturated = NULL,
    missingTumorIndicator = 888
  )
  score.test.icog.baseline<- ScoreTestSelfDesign(y=y.pheno.mis1,
                                                 x=x.all.mis1[,1,drop=F],
                                                 z.design=z.design.score.baseline,
                                                 score.test.support=score.test.support.icog.baseline,
                                                 missingTumorIndicator=888)
  
  score.icog.baseline <- score.test.icog.baseline[[1]]
  infor.icog.baseline <- score.test.icog.baseline[[2]]
  
  
  
  
  score.test.support.icog.casecase <- ScoreTestSupport(
    y.pheno.mis1,
    baselineonly = x.all.mis1[,1,drop=F],
    additive = x.all.mis1[,2:11],
    pairwise.interaction = NULL,
    saturated = NULL,
    missingTumorIndicator = 888
  )
  score.test.icog.casecase<- ScoreTestSelfDesign(y=y.pheno.mis1,
                         x=x.all.mis1[,1,drop=F],
                        z.design=z.design.score.casecase,
                         score.test.support=score.test.support.icog.casecase,
                         missingTumorIndicator=888)
  
  score.icog.casecase <- score.test.icog.casecase[[1]]
  infor.icog.casecase <- score.test.icog.casecase[[2]]
  
  DisplayMixedScoreTestResult(score.icog.baseline,infor.icog.baseline,score.icog.casecase,infor.icog.casecase)
  
  
  
  
  
  #analysis for Onco Array
  #data2 <- read.csv("./V10/Onco_euro_v10_05242017.csv",header=T)
  data2 <- read.csv("./data/Onco_euro_v10_05242017.csv",header=T)
  names1 = colnames(data1)[27:206]
  rm(data1)
  names2 = colnames(data2)[27:212]
  
  idxi1 = which(names2==names1[i1])
  y.pheno.mis2 <- cbind(data2$Behaviour1,data2$PR_status1,data2$ER_status1,data2$HER2_status1)
  #y.pheno.mis2 <- cbind(data2$Behaviour1,data2$PR_status1,data2$ER_status1,data2$HER2_status1)
  colnames(y.pheno.mis2) = c("Behaviour","PR",
                             "ER","HER2")
  
  x.test.all.mis2 <- data2[,c(27:212)]
  x.covar.mis2 <- data2[,5:14]
  x.all.mis2 <- as.matrix(cbind(x.test.all.mis2[,idxi1],x.covar.mis2))
  colnames(x.all.mis2)[1] = "gene"
  
  
  
  Heter.result.Onco = EMmvpoly(y.pheno.mis2,baselineonly = NULL,additive = x.all.mis2,pairwise.interaction = NULL,saturated = NULL,missingTumorIndicator = 888)
  
  z.standard <- Heter.result.Onco[[12]]
  z.additive.design <- as.matrix(cbind(1,z.standard))
  M <- nrow(z.standard)
  number.of.tumor <- ncol(z.standard)
  log.odds.onco <- Heter.result.Onco[[1]][(M+1):(M+1+number.of.tumor)]
  sigma.log.odds.onco <- Heter.result.Onco[[2]][(M+1):(M+1+number.of.tumor),(M+1):(M+1+number.of.tumor)]
  beta.onco <- z.additive.design%*%log.odds.onco
  beta.sigma.onco <- z.additive.design%*%sigma.log.odds.onco%*%t(z.additive.design)
  loglikelihood.onco <- Heter.result.Onco[[8]]
  
  score.test.support.onco.baseline <- ScoreTestSupport(
    y.pheno.mis2,
    baselineonly = NULL,
    additive = x.all.mis2[,2:11],
    pairwise.interaction = NULL,
    saturated = NULL,
    missingTumorIndicator = 888
  )
  score.test.onco.baseline<- ScoreTestSelfDesign(y=y.pheno.mis2,
                                                 x=x.all.mis2[,1,drop=F],
                                                 z.design=z.design.score.baseline,
                                                 score.test.support=score.test.support.onco.baseline,
                                                 missingTumorIndicator=888)
  
  
  
  
  score.test.support.onco.casecase <- ScoreTestSupport(
    y.pheno.mis2,
    baselineonly = x.all.mis2[,1,drop=F],
    additive = x.all.mis2[,2:11],
    pairwise.interaction = NULL,
    saturated = NULL,
    missingTumorIndicator = 888
  )
  score.test.onco.casecase<- ScoreTestSelfDesign(y=y.pheno.mis2,
                                                 x=x.all.mis2[,1,drop=F],
                                                 z.design=z.design.score.casecase,
                                                 score.test.support=score.test.support.onco.casecase,
                                                 missingTumorIndicator=888)
  
  
  
  
  
  
  meta.result <- LogoddsMetaAnalysis(log.odds.icog,
                                     sigma.log.odds.icog,
                                     log.odds.onco,
                                     sigma.log.odds.onco)
  logodds <- meta.result[[1]]
  sigma <- meta.result[[2]]
  
  score.meta.baseline <- ScoreMetaAnalysis(score)
  
  
  
  test.result <- DisplayTestResult(logodds,sigma)
  
  meta.result.beta <- LogoddsMetaAnalysis(beta.icog,
                                          beta.sigma.icog,
                                          beta.onco,
                                          beta.sigma.onco)
  
  heter.result <- list(test.result = test.result,
                       logodds = logodds,
                       sigma= sigma)
  save(heter.result,file=paste0("./known_SNPs_analysis_revised/additive_model/result/heter_result_",i1,".Rdata"))
  
  
  
}else{
  data2 <- read.csv("./data/Onco_euro_v10_rs554219.csv",header = T)
  #data2 <- read.csv("/data/zhangh20/breast_cancer/known_SNPs_anlysis/Onco_euro_v10_rs554219.csv",1)
  
  # names1 = colnames(data1)[27:206]
  #rm(data1)
  names2 = colnames(data2)
  
  idxi1 = which(names2=="rs554219")
  y.pheno.mis2 <- cbind(data2$Behaviour1,data2$PR_status1,data2$ER_status1,data2$HER2_status1)
  x.test.all.mis2 <- data2
  
  x.covar.mis2 <- data2[,5:14]
  x.all.mis2 <- as.matrix(cbind(x.test.all.mis2[,idxi1],x.covar.mis2))
  
  
  
  Heter.result.Onco = EMmvpoly(y.pheno.mis2,baselineonly = NULL,additive = x.all.mis2,pairwise.interaction = NULL,saturated = NULL,missingTumorIndicator = 888)
  
  
  
  M <- Heter.result.Onco[[5]]
  number.of.tumor <- Heter.result.Onco[[6]]
  log.odds.onco <- Heter.result.Onco[[1]][(M+1):(M+1+number.of.tumor)]
  info.onco <- Heter.result.Onco[[2]]
  sigma.onco <- solve(info.onco)
  sigma.log.odds.onco <- sigma.onco[(M+1):(M+1+number.of.tumor),(M+1):(M+1+number.of.tumor)]
  
  logodds <- log.odds.onco
  sigma <- sigma.log.odds.onco
  
  test.result <- DisplayTestResult(logodds,sigma)
  
  heter.result <- list(test.result = test.result,
                       logodds = logodds,
                       sigma= sigma)
  save(heter.result,file=paste0("./known_SNPs_analysis_revised/additive_model/result/heter_result_",i1,".Rdata"))
  
}



