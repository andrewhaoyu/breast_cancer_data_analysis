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

z.design <- matrix(c(
  c(0,1,1,1,0,0,0,0,1,1,1,0,0,0,0,0,1,1,0,0,0,0,0),
  c(0,0,0,0,0,1,1,0,0,0,0,0,1,1,1,0,0,0,0,0,1,1,1),
  c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0),
  c(0,0,0,0,1,0,0,0,0,0,0,1,0,0,0,0,0,0,0,1,0,0,0),
  c(1,0,0,0,0,0,0,1,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0)
),ncol=5)
colnames(z.design) <- c("Luminial A","Luminal B",
                        "Luminal B HER2-",
                        "HER2 Enriched",
                        "Triple Negative")


if(i1<=180){
  ##analysis for Icog
  data1 <- read.csv("./data/iCOGS_euro_v10_05242017.csv",header=T)
  y.pheno.mis1 <- cbind(data1$Behaviour1,data1$PR_status1,data1$ER_status1,data1$HER2_status1,data1$Grade1)
  colnames(y.pheno.mis1) = c("Behavior","PR","ER","HER2","Grade")
  # Grade1.fake <- data1$Grade1
  # Grade1.fake[data1$Grade1==2|data1$Grade1==3] <- 1
  # Grade1.fake[data1$Grade1==1] <- 0
  #y.pheno.mis1 <- cbind(data1$Behaviour1,data1$PR_status1,data1$ER_status1,data1$HER2_status1,Grade1.fake)
  # y.pheno.mis1 <- cbind(data1$Behaviour1,data1$PR_status1,data1$ER_status1,data1$HER2_status1)
  
  x.test.all.mis1 <- data1[,c(27:206)]
  
  x.covar.mis1 <- data1[,5:14]
  x.all.mis1 <- as.matrix(cbind(x.test.all.mis1[,i1],x.covar.mis1))
  colnames(x.all.mis1)[1] <- "gene"
  
  Heter.result.Icog = EMmvpolySelfDesign(y.pheno.mis1,x.self.design = x.all.mis1[,1,drop=F],z.design=z.design,baselineonly = NULL,additive = x.all.mis1[,2:11],pairwise.interaction = NULL,saturated = NULL,missingTumorIndicator = 888)
  z.standard <- Heter.result.Icog[[12]]
  z.additive.design <- as.matrix(cbind(1,z.standard))
  M <- nrow(z.standard)
  number.of.tumor <- ncol(z.standard)
  log.odds.icog <- Heter.result.Icog[[1]][(M+1):(M+1+number.of.tumor)]
  nparm <- length(Heter.result.Icog[[1]])  
  sigma.log.odds.icog <- Heter.result.Icog[[2]][(M+1):(M+1+number.of.tumor),(M+1):(M+1+number.of.tumor)]
  beta.icog <- z.design%*%log.odds.icog
  beta.sigma.icog <- z.design%*%sigma.log.odds.icog%*%t(z.design)
  loglikelihood.icog <- Heter.result.Icog[[8]]
  rm(Heter.result.Icog)  
  score.test.support.icog <- ScoreTestSupport(
    y.pheno.mis1,
    baselineonly = NULL,
    additive = x.all.mis1[,2:11],
    pairwise.interaction = NULL,
    saturated = NULL,
    missingTumorIndicator = 888
  )
  
  
  
  score.test.icog<- ScoreTestSelfDesign(y=y.pheno.mis1,
                                        x=x.all.mis1[,1,drop=F],
                                        z.design=z.design,
                                        score.test.support=score.test.support.icog,
                                        missingTumorIndicator=888)
  z.design.score.baseline <- matrix(rep(1,M),ncol=1)
  z.design.score.casecase <-z.standard
  
  score.icog <- score.test.icog[[1]]
  infor.icog <- score.test.icog[[2]]
  rm(score.test.support.icog)
  
  
  #analysis for Onco Array
  #data2 <- read.csv("./V10/Onco_euro_v10_05242017.csv",header=T)
  data2 <- read.csv("./data/Onco_euro_v10_05242017.csv",header=T)
  names1 = colnames(data1)[27:206]
  
  names2 = colnames(data2)[27:212]
  
  idxi1 = which(names2==names1[i1])
  y.pheno.mis2 <- cbind(data2$Behaviour1,data2$PR_status1,data2$ER_status1,data2$HER2_status1,data2$Grade1)
  #y.pheno.mis2 <- cbind(data2$Behaviour1,data2$PR_status1,data2$ER_status1,data2$HER2_status1)
  colnames(y.pheno.mis2) = c("Behaviour","PR",
                             "ER","HER2","Grade")
  
  x.test.all.mis2 <- data2[,c(27:212)]
  x.covar.mis2 <- data2[,5:14]
  x.all.mis2 <- as.matrix(cbind(x.test.all.mis2[,idxi1],x.covar.mis2))
  colnames(x.all.mis2)[1] = "gene"
  
  
  Heter.result.Onco = EMmvpolySelfDesign(y.pheno.mis2,x.self.design = x.all.mis2[,1,drop=F],z.design = z.design,baselineonly = NULL,additive = x.all.mis2[,2:11],pairwise.interaction = NULL,saturated = NULL,missingTumorIndicator = 888)
  z.standard <- Heter.result.Onco[[12]]
  M <- nrow(z.standard)
  number.of.tumor <- ncol(z.standard)
  log.odds.onco <- Heter.result.Onco[[1]][(M+1):(M+1+number.of.tumor)]
  nparm <- length(Heter.result.Onco[[1]])
  sigma.log.odds.onco <- Heter.result.Onco[[2]][(M+1):(M+1+number.of.tumor),(M+1):(M+1+number.of.tumor)]
  beta.onco <- z.design%*%log.odds.onco
  beta.sigma.onco <- z.design%*%sigma.log.odds.onco%*%t(z.design)
  loglikelihood.onco <- Heter.result.Onco[[8]]
  
  
  score.test.support.onco <- ScoreTestSupport(
    y.pheno.mis2,
    baselineonly = NULL,
    additive = x.all.mis2[,2:11],
    pairwise.interaction = NULL,
    saturated = NULL,
    missingTumorIndicator = 888
  )
  
  score.test.onco<- ScoreTestSelfDesign(y=y.pheno.mis2,
                                        x=x.all.mis2[,1,drop=F],
                                        z.design= z.design,
                                        score.test.support=score.test.support.onco,
                                        missingTumorIndicator=888)
  
  score.onco <- score.test.onco[[1]]
  infor.onco <- score.test.onco[[2]]
  
  
  
  meta.result <- LogoddsMetaAnalysis(log.odds.icog,
                                     sigma.log.odds.icog,
                                     log.odds.onco,
                                     sigma.log.odds.onco)
  
  second.stage.logodds.meta <- meta.result[[1]]
  second.stage.sigma.meta <- meta.result[[2]]
  
  
  
  test.result.second.wald <- DisplaySecondStageTestResult(second.stage.logodds.meta,second.stage.sigma.meta,self.design=T)
  
  
  beta.meta <- z.design%*%second.stage.logodds.meta
  beta.sigma.meta <- z.design%*%second.stage.sigma.meta%*%t(z.design)
  
  test.result.first.wald <- DisplayFirstStageTestResult(beta.meta,beta.sigma.meta)
  
  meta.result.score <- ScoreMetaAnalysis(score.icog,infor.icog,
                                         score.onco,infor.onco)
  score.meta <- meta.result.score[[1]]
  infor.meta <- meta.result.score[[2]]
  
  test.result.second.score <- ScoreGlobalMixedTestForAssoc(score.meta,infor.meta)
  test.result.second.score <- t(test.result.second.score)
  
  
  colnames(test.result.second.score) <- c("fixed effect score test global test for association","random effect model global test for association")
  
  loglikelihood <- loglikelihood.icog+loglikelihood.onco
  AIC <- 2*nparm-2*loglikelihood
  
  heter.result <- list(data.frame(test.result.second.wald,test.result.second.score,loglikelihood = loglikelihood,AIC=AIC),
                       data.frame(test.result.first.wald))
  
  save(heter.result,file=paste0("./known_SNPs/known_SNPs_analysis_G_revised/intrinsic_subtypes_pc_additive/result/heter_result_",i1,".Rdata"))
  
  
  
}else{
  data2 <- read.csv("./data/Onco_euro_v10_rs554219.csv",header = T)
  #data2 <- read.csv("/data/zhangh20/breast_cancer/known_SNPs_anlysis/Onco_euro_v10_rs554219.csv",1)
  
  # names1 = colnames(data1)[27:206]
  #rm(data1)
  names2 = colnames(data2)
  
  idxi1 = which(names2=="rs554219")
  y.pheno.mis2 <- cbind(data2$Behaviour1,data2$PR_status1,data2$ER_status1,data2$HER2_status1,data2$Grade1)
  x.test.all.mis2 <- data2
  colnames(y.pheno.mis2) <- c("Behavior","PR","ER","HER2","Grade")
  
  x.covar.mis2 <- data2[,5:14]
  x.all.mis2 <- as.matrix(cbind(x.test.all.mis2[,idxi1],x.covar.mis2))
  
  
  
  Heter.result.Onco = EMmvpolySelfDesign(y.pheno.mis2,x.self.design = x.all.mis2[,1,drop=F],z.design = z.design,baselineonly = NULL,additive = x.all.mis2[,2:11],pairwise.interaction = NULL,saturated = NULL,missingTumorIndicator = 888)
  z.standard <- Heter.result.Onco[[12]]
  M <- nrow(z.standard)
  number.of.tumor <- ncol(z.standard)
  log.odds.onco <- Heter.result.Onco[[1]][(M+1):(M+1+number.of.tumor)]
  nparm = length(Heter.result.Onco[[1]])
  sigma.log.odds.onco <- Heter.result.Onco[[2]][(M+1):(M+1+number.of.tumor),(M+1):(M+1+number.of.tumor)]
  beta.onco <- z.design%*%log.odds.onco
  beta.sigma.onco <- z.design%*%sigma.log.odds.onco%*%t(z.design)
  loglikelihood.onco <- Heter.result.Onco[[8]]
  
  
  score.test.support.onco <- ScoreTestSupport(
    y.pheno.mis2,
    baselineonly = NULL,
    additive = x.all.mis2[,2:11],
    pairwise.interaction = NULL,
    saturated = NULL,
    missingTumorIndicator = 888
  )
  
  score.test.onco<- ScoreTestSelfDesign(y=y.pheno.mis2,
                                        x=x.all.mis2[,1,drop=F],
                                        z.design= z.design,
                                        score.test.support=score.test.support.onco,
                                        missingTumorIndicator=888)
  
  score.onco <- score.test.onco[[1]]
  infor.onco <- score.test.onco[[2]]
  
  
  
  
  
  second.stage.logodds.meta <-log.odds.onco
  second.stage.sigma.meta <- sigma.log.odds.onco
  
  
  
  
  
  test.result.second.wald <- DisplaySecondStageTestResult(second.stage.logodds.meta,second.stage.sigma.meta,self.design=T)
  
  
  
  beta.meta <- z.design%*%second.stage.logodds.meta
  beta.sigma.meta <- z.design%*%second.stage.sigma.meta%*%t(z.design)
  
  test.result.first.wald <- DisplayFirstStageTestResult(beta.meta,beta.sigma.meta)
  
  
  score.meta <- score.onco
  infor.meta <- infor.onco
  
  test.result.second.score <- ScoreGlobalMixedTestForAssoc(score.meta,infor.meta)
  test.result.second.score <- t(test.result.second.score)
  colnames(test.result.second.score) <- c("fixed effect score test global test for association","random effect model global test for association")
  
  
  loglikelihood <- loglikelihood.onco
  AIC <- 2*nparm-2*loglikelihood
  
  heter.result <- list(data.frame(test.result.second.wald,test.result.second.score,loglikelihood = loglikelihood,AIC=AIC),
                       data.frame(test.result.first.wald))
  
  
  save(heter.result,file=paste0("./known_SNPs/known_SNPs_analysis_G_revised/intrinsic_subtypes_pc_additive/result/heter_result_",i1,".Rdata"))
  
  
  
  
  
}



