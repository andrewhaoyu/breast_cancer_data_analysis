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

z.design <- matrix(c(
  c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1),
  c(0,0,0,0,0,1,1,0,0,0,0,0,1,1,1,0,0,0,0,0,1,1,1),
  c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0),
  c(0,0,0,0,1,0,0,0,0,0,0,1,0,0,0,0,0,0,0,1,0,0,0),
  c(1,0,0,0,0,0,0,1,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0)
),ncol=5)
z.design.score.baseline <- z.design[,1,drop=F]
z.design.score.casecase <- z.design[,2:ncol(z.design)]
z.trans <- matrix(c(
  c(1,1,1,1,1),
  c(0,1,0,0,0),
  c(0,0,1,0,0),
  c(0,0,0,1,0),
  c(0,0,0,0,1)
),ncol = 5
)


colnames(z.design) <- c("Luminial A","Luminal B",
                        "Luminal B HER2Neg - Luminal A",
                        "HER2 Enriched - Luminal A",
                        "Triple Negative - Luminal A")


if(i1<=177){
  data1 <- fread("./data/iCOGS_euro_v10_10232017.csv",header=T)
  data1 <- as.data.frame(data1)
  y.pheno.mis1 <- cbind(data1$Behaviour1,data1$PR_status1,data1$ER_status1,data1$HER2_status1,data1$Grade1)
  colnames(y.pheno.mis1) = c("Behavior","PR","ER","HER2","Grade")
  # Grade1.fake <- data1$Grade1
  # Grade1.fake[data1$Grade1==2|data1$Grade1==3] <- 1
  # Grade1.fake[data1$Grade1==1] <- 0
  #y.pheno.mis1 <- cbind(data1$Behaviour1,data1$PR_status1,data1$ER_status1,data1$HER2_status1,Grade1.fake)
  # y.pheno.mis1 <- cbind(data1$Behaviour1,data1$PR_status1,data1$ER_status1,data1$HER2_status1)
  
  x.test.all.mis1 <- data1[,c(27:203)]
  x.test.all.mis1 <- x.test.all.mis1
  ###pc1-10 and age
  x.covar.mis1 <- data1[,c(5:14)]
  
  
  
  x.all.mis1 <- as.matrix(cbind(x.test.all.mis1[,i1],x.covar.mis1))
  colnames(x.all.mis1)[1] <- "gene"
  
  Heter.result.Icog = EMmvpolySelfDesign(y.pheno.mis1,x.self.design = x.all.mis1[,1,drop=F],z.design=z.design,baselineonly = NULL,additive = x.all.mis1[,2:ncol(x.all.mis1)],pairwise.interaction = NULL,saturated = NULL,missingTumorIndicator = 888)
  z.standard <- Heter.result.Icog[[12]]
  M <- nrow(z.standard)
  z.design.PR <- matrix(0,M,2)
  idx.PR.pos <- which(z.standard[,1]==1)
  z.design.PR[idx.PR.pos,1] <- 1
  idx.PR.neg <- which(z.standard[,1]==0)
  z.design.PR[idx.PR.neg,2] <- 1
  colnames(z.design.PR) <- c("PRPos","PRNeg")
  
  z.design.ER <- matrix(0,M,2)
  idx.ER.pos <- which(z.standard[,2]==1)
  z.design.ER[idx.ER.pos,1] <- 1
  idx.ER.neg <- which(z.standard[,2]==0)
  z.design.ER[idx.ER.neg,2] <- 1
  colnames(z.design.ER) <- c("ERPos","ERNeg")
  
  z.design.HER <- matrix(0,M,2)
  idx.HER.pos <- which(z.standard[,3]==1)
  z.design.HER[idx.HER.pos,1] <- 1
  idx.HER.neg <- which(z.standard[,3]==0)
  z.design.HER[idx.HER.neg,2] <- 1
  colnames(z.design.HER) <- c("HERPos","HERNeg")
  
  z.design.Grade <- matrix(0,M,2)
  z.design.Grade[,1] <- 1
  idx.Grade1 <- which(z.standard[,4]==1)
  z.design.Grade[idx.Grade1,2] <- 1
  idx.Grade2 <- which(z.standard[,4]==2)
  z.design.Grade[idx.Grade2,2] <- 2
  idx.Grade3 <- which(z.standard[,4]==3)
  z.design.Grade[idx.Grade3,2] <- 3
  colnames(z.design.Grade) <- c("Grade1","Gradecasecase")
  
  
  z.additive.design <- as.matrix(cbind(1,z.standard))
  M <- nrow(z.standard)
  number.of.tumor <- ncol(z.standard)
  log.odds.icog <- Heter.result.Icog[[1]][(M+1):(M+1+number.of.tumor)]
  nparm <- length(Heter.result.Icog[[1]])  
  sigma.log.odds.icog <- Heter.result.Icog[[2]][(M+1):(M+1+number.of.tumor),(M+1):(M+1+number.of.tumor)]
  beta.icog <- z.trans%*%log.odds.icog
  beta.sigma.icog <- z.trans%*%sigma.log.odds.icog%*%t(z.trans)
  loglikelihood.icog <- Heter.result.Icog[[8]]
  
  Heter.result.PR.Icog = EMmvpolySelfDesign(y.pheno.mis1,x.self.design = x.all.mis1[,1,drop=F],z.design=z.design.PR,baselineonly = NULL,additive = x.all.mis1[,2:ncol(x.all.mis1)],pairwise.interaction = NULL,saturated = NULL,missingTumorIndicator = 888)
  log.odds.PR.Icog <-   Heter.result.PR.Icog[[1]][(M+1):(M+2)]
  sigma.log.odds.PR.Icog <- Heter.result.PR.Icog[[2]][(M+1):(M+2),(M+1):(M+2)]
  
  Heter.result.ER.Icog = EMmvpolySelfDesign(y.pheno.mis1,x.self.design = x.all.mis1[,1,drop=F],z.design=z.design.ER,baselineonly = NULL,additive = x.all.mis1[,2:ncol(x.all.mis1)],pairwise.interaction = NULL,saturated = NULL,missingTumorIndicator = 888)
  log.odds.ER.Icog <-   Heter.result.ER.Icog[[1]][(M+1):(M+2)]
  sigma.log.odds.ER.Icog <- Heter.result.ER.Icog[[2]][(M+1):(M+2),(M+1):(M+2)]
  
  Heter.result.HER.Icog = EMmvpolySelfDesign(y.pheno.mis1,x.self.design = x.all.mis1[,1,drop=F],z.design=z.design.HER,baselineonly = NULL,additive = x.all.mis1[,2:ncol(x.all.mis1)],pairwise.interaction = NULL,saturated = NULL,missingTumorIndicator = 888)
  log.odds.HER.Icog <-   Heter.result.HER.Icog[[1]][(M+1):(M+2)]
  sigma.log.odds.HER.Icog <- Heter.result.HER.Icog[[2]][(M+1):(M+2),(M+1):(M+2)]
  
  Heter.result.Grade.Icog = EMmvpolySelfDesign(y.pheno.mis1,x.self.design = x.all.mis1[,1,drop=F],z.design=z.design.Grade,baselineonly = NULL,additive = x.all.mis1[,2:ncol(x.all.mis1)],pairwise.interaction = NULL,saturated = NULL,missingTumorIndicator = 888)
  log.odds.Grade.Icog <-   Heter.result.Grade.Icog[[1]][(M+1):(M+2)]
  sigma.log.odds.Grade.Icog <- Heter.result.Grade.Icog[[2]][(M+1):(M+2),(M+1):(M+2)]
  
  
  beta.grade.icog <- z.design.Grade%*%log.odds.Grade.Icog
  beta.sigma.grade.icog <- z.design.Grade%*%sigma.log.odds.Grade.Icog%*%t(z.design.Grade)
  
  p.grade <- c(1,8,16)
  beta.grade.icog <- beta.grade.icog[p.grade]  
  beta.sigma.grade.icog <- beta.sigma.grade.icog[p.grade,p.grade]  
  
  
  
  
  
  
  
  
  
  
  names1 = colnames(data1)[27:206]
  rm(data1)
  gc()
  
  
  
  #analysis for Onco Array
  #data2 <- read.csv("./V10/Onco_euro_v10_05242017.csv",header=T)
  data2 <- fread("./data/Onco_euro_v10_10232017.csv",header=T)
  data2 <- as.data.frame(data2)
  y.pheno.mis2 <- cbind(data2$Behaviour1,data2$PR_status1,data2$ER_status1,data2$HER2_status1,data2$Grade1)
  #y.pheno.mis2 <- cbind(data2$Behaviour1,data2$PR_status1,data2$ER_status1,data2$HER2_status1)
  colnames(y.pheno.mis2) = c("Behaviour","PR",
                             "ER","HER2","Grade")
  
  x.test.all.mis2 <- data2[,c(27:203)]
  x.test.all.mis2 <- 2-x.test.all.mis2
  x.covar.mis2 <- data2[,c(5:14,204)]
  
  
  x.all.mis2 <- as.matrix(cbind(x.test.all.mis2[,i1],x.covar.mis2))
  colnames(x.all.mis2)[1] = "gene"
  
  
  
  
  Heter.result.Onco = EMmvpolySelfDesign(y.pheno.mis2,x.self.design = x.all.mis2[,1,drop=F],z.design = z.design,baselineonly = NULL,additive = x.all.mis2[,2:ncol(x.all.mis2)],pairwise.interaction = NULL,saturated = NULL,missingTumorIndicator = 888)
  z.standard <- Heter.result.Onco[[12]]
  M <- nrow(z.standard)
  number.of.tumor <- ncol(z.standard)
  log.odds.onco <- Heter.result.Onco[[1]][(M+1):(M+1+number.of.tumor)]
  nparm <- length(Heter.result.Onco[[1]])
  sigma.log.odds.onco <- Heter.result.Onco[[2]][(M+1):(M+1+number.of.tumor),(M+1):(M+1+number.of.tumor)]
  beta.onco <- z.trans%*%log.odds.onco
  beta.sigma.onco <- z.trans%*%sigma.log.odds.onco%*%t(z.trans)
  loglikelihood.onco <- Heter.result.Onco[[8]]
  
  
  Heter.result.PR.Onco = EMmvpolySelfDesign(y.pheno.mis2,x.self.design = x.all.mis2[,1,drop=F],z.design=z.design.PR,baselineonly = NULL,additive = x.all.mis2[,2:ncol(x.all.mis2)],pairwise.interaction = NULL,saturated = NULL,missingTumorIndicator = 888)
  log.odds.PR.Onco <-   Heter.result.PR.Onco[[1]][(M+1):(M+2)]
  sigma.log.odds.PR.Onco <- Heter.result.PR.Onco[[2]][(M+1):(M+2),(M+1):(M+2)]
  
  Heter.result.ER.Onco = EMmvpolySelfDesign(y.pheno.mis2,x.self.design = x.all.mis2[,1,drop=F],z.design=z.design.ER,baselineonly = NULL,additive = x.all.mis2[,2:ncol(x.all.mis2)],pairwise.interaction = NULL,saturated = NULL,missingTumorIndicator = 888)
  log.odds.ER.Onco <-   Heter.result.ER.Onco[[1]][(M+1):(M+2)]
  sigma.log.odds.ER.Onco <- Heter.result.ER.Onco[[2]][(M+1):(M+2),(M+1):(M+2)]
  
  Heter.result.HER.Onco = EMmvpolySelfDesign(y.pheno.mis2,x.self.design = x.all.mis2[,1,drop=F],z.design=z.design.HER,baselineonly = NULL,additive = x.all.mis2[,2:ncol(x.all.mis2)],pairwise.interaction = NULL,saturated = NULL,missingTumorIndicator = 888)
  log.odds.HER.Onco <-   Heter.result.HER.Onco[[1]][(M+1):(M+2)]
  sigma.log.odds.HER.Onco <- Heter.result.HER.Onco[[2]][(M+1):(M+2),(M+1):(M+2)]
  
  Heter.result.Grade.Onco = EMmvpolySelfDesign(y.pheno.mis2,x.self.design = x.all.mis2[,1,drop=F],z.design=z.design.Grade,baselineonly = NULL,additive = x.all.mis2[,2:ncol(x.all.mis2)],pairwise.interaction = NULL,saturated = NULL,missingTumorIndicator = 888)
  log.odds.Grade.Onco <-   Heter.result.Grade.Onco[[1]][(M+1):(M+3)]
  sigma.log.odds.Grade.Onco <- Heter.result.Grade.Onco[[2]][(M+1):(M+3),(M+1):(M+3)]
  beta.grade.onco <- z.design.Grade%*%log.odds.Grade.Onco
  beta.sigma.onco <- z.design.Grade%*%sigma.log.odds.Grade.Onco%*%t(z.design.Grade)
  
  p.grade <- c(1,8,16)
  beta.grade.onco <- beta.grade.onco[p.grade]  
  beta.sigma.grade.onco <- beta.sigma.grade.onco[p.grade,p.grade]  
  
  
  
  
  
  
  # 
  # meta.result <- LogoddsMetaAnalysis(log.odds.icog,
  #                                    sigma.log.odds.icog,
  #                                    log.odds.onco,
  #                                    sigma.log.odds.onco)
  # 
  # second.stage.logodds.meta <- meta.result[[1]]
  # second.stage.sigma.meta <- meta.result[[2]]
  # 
  # 
  # 
  # test.result.second.wald <- DisplaySecondStageTestResult(second.stage.logodds.meta,second.stage.sigma.meta)
  
  meta.result.PR <- LogoddsMetaAnalysis(log.odds.PR.Icog,
                                        sigma.log.odds.PR.Icog,
                                        log.odds.PR.Onco,
                                        sigma.log.odds.PR.Onco)
  
  second.stage.logodds.meta.PR <- meta.result.PR[[1]]
  second.stage.sigma.meta.PR <- meta.result.PR[[2]]
  
  test.result.second.wald.PR <- DisplaySecondStageTestResult(second.stage.logodds.meta.PR,second.stage.sigma.meta.PR)
  
  
  meta.result.ER <- LogoddsMetaAnalysis(log.odds.ER.Icog,
                                        sigma.log.odds.ER.Icog,
                                        log.odds.ER.Onco,
                                        sigma.log.odds.ER.Onco)
  
  second.stage.logodds.meta.ER <- meta.result.ER[[1]]
  second.stage.sigma.meta.ER <- meta.result.ER[[2]]
  
  test.result.second.wald.ER <- DisplaySecondStageTestResult(second.stage.logodds.meta.ER,second.stage.sigma.meta.ER)
  
  
  meta.result.HER <- LogoddsMetaAnalysis(log.odds.HER.Icog,
                                         sigma.log.odds.HER.Icog,
                                         log.odds.HER.Onco,
                                         sigma.log.odds.HER.Onco)
  
  second.stage.logodds.meta.HER <- meta.result.HER[[1]]
  second.stage.sigma.meta.HER <- meta.result.HER[[2]]
  
  test.result.second.wald.HER <- DisplaySecondStageTestResult(second.stage.logodds.meta.HER,second.stage.sigma.meta.HER)
  
  meta.result.Grade <- LogoddsMetaAnalysis(log.odds.Grade.Icog,
                                           sigma.log.odds.Grade.Icog,
                                           log.odds.Grade.Onco,
                                           sigma.log.odds.Grade.Onco)
  
  second.stage.logodds.meta.Grade <- meta.result.Grade[[1]]
  second.stage.sigma.meta.Grade <- meta.result.Grade[[2]]
  
  test.result.second.wald.Grade <- DisplaySecondStageTestResult(second.stage.logodds.meta.Grade,second.stage.sigma.meta.Grade)
  
  
  
  
  
  
  
  
  
  heter.result <- list(data.frame(test.result.second.wald.PR,test.result.second.wald.ER,test.result.second.wald.HER,test.result.second.wald.Grade))
  
  
  
  
  save(heter.result,file=paste0("./known_SNPs/known_SNPs_analysis_G_revised/ERPosvesNeg/result/heter_result_",i1,".Rdata"))
  
  
  
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
  
  
  x.all.mis2 <- as.matrix(cbind(x.test.all.mis2[,idxi1],x.covar.mis2))
  colnames(x.all.mis2)[1] = "gene"
  
  
  
  
  
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
  
  
  M <- nrow(z.standard)
  z.design.PR <- matrix(0,M,2)
  idx.PR.pos <- which(z.standard[,1]==1)
  z.design.PR[idx.PR.pos,1] <- 1
  idx.PR.neg <- which(z.standard[,1]==0)
  z.design.PR[idx.PR.neg,2] <- 1
  colnames(z.design.PR) <- c("PRPos","PRNeg")
  
  z.design.ER <- matrix(0,M,2)
  idx.ER.pos <- which(z.standard[,2]==1)
  z.design.ER[idx.ER.pos,1] <- 1
  idx.ER.neg <- which(z.standard[,2]==0)
  z.design.ER[idx.ER.neg,2] <- 1
  colnames(z.design.ER) <- c("ERPos","ERNeg")
  
  z.design.HER <- matrix(0,M,2)
  idx.HER.pos <- which(z.standard[,3]==1)
  z.design.HER[idx.HER.pos,1] <- 1
  idx.HER.neg <- which(z.standard[,3]==0)
  z.design.HER[idx.HER.neg,2] <- 1
  colnames(z.design.HER) <- c("HERPos","HERNeg")
  
  z.design.Grade <- matrix(0,M,2)
  z.design.Grade[,1] <- 1
  idx.Grade1 <- which(z.standard[,4]==1)
  z.design.Grade[idx.Grade1,2] <- 1
  idx.Grade2 <- which(z.standard[,4]==2)
  z.design.Grade[idx.Grade2,2] <- 2
  idx.Grade3 <- which(z.standard[,4]==3)
  z.design.Grade[idx.Grade3,2] <- 3
  colnames(z.design.Grade) <- c("Grade1","Gradecasecase")
  
  
  
  Heter.result.PR.Onco = EMmvpolySelfDesign(y.pheno.mis2,x.self.design = x.all.mis2[,1,drop=F],z.design=z.design.PR,baselineonly = NULL,additive = x.all.mis2[,2:ncol(x.all.mis2)],pairwise.interaction = NULL,saturated = NULL,missingTumorIndicator = 888)
  log.odds.PR.Onco <-   Heter.result.PR.Onco[[1]][(M+1):(M+2)]
  sigma.log.odds.PR.Onco <- Heter.result.PR.Onco[[2]][(M+1):(M+2),(M+1):(M+2)]
  
  Heter.result.ER.Onco = EMmvpolySelfDesign(y.pheno.mis2,x.self.design = x.all.mis2[,1,drop=F],z.design=z.design.ER,baselineonly = NULL,additive = x.all.mis2[,2:ncol(x.all.mis2)],pairwise.interaction = NULL,saturated = NULL,missingTumorIndicator = 888)
  log.odds.ER.Onco <-   Heter.result.ER.Onco[[1]][(M+1):(M+2)]
  sigma.log.odds.ER.Onco <- Heter.result.ER.Onco[[2]][(M+1):(M+2),(M+1):(M+2)]
  
  Heter.result.HER.Onco = EMmvpolySelfDesign(y.pheno.mis2,x.self.design = x.all.mis2[,1,drop=F],z.design=z.design.HER,baselineonly = NULL,additive = x.all.mis2[,2:ncol(x.all.mis2)],pairwise.interaction = NULL,saturated = NULL,missingTumorIndicator = 888)
  log.odds.HER.Onco <-   Heter.result.HER.Onco[[1]][(M+1):(M+2)]
  sigma.log.odds.HER.Onco <- Heter.result.HER.Onco[[2]][(M+1):(M+2),(M+1):(M+2)]
  
  Heter.result.Grade.Onco = EMmvpolySelfDesign(y.pheno.mis2,x.self.design = x.all.mis2[,1,drop=F],z.design=z.design.Grade,baselineonly = NULL,additive = x.all.mis2[,2:ncol(x.all.mis2)],pairwise.interaction = NULL,saturated = NULL,missingTumorIndicator = 888)
  log.odds.Grade.Onco <-   Heter.result.Grade.Onco[[1]][(M+1):(M+3)]
  sigma.log.odds.Grade.Onco <- Heter.result.Grade.Onco[[2]][(M+1):(M+3),(M+1):(M+3)]
  
  
  
  
  
  
  
  
  
  
  second.stage.logodds.meta.PR <- log.odds.PR.Onco
  second.stage.sigma.meta.PR <- sigma.log.odds.PR.Onco
  
  test.result.second.wald.PR <- DisplaySecondStageTestResult(second.stage.logodds.meta.PR,second.stage.sigma.meta.PR)
  
  
  
  
  second.stage.logodds.meta.ER <- log.odds.ER.Onco
  second.stage.sigma.meta.ER <- sigma.log.odds.ER.Onco
  
  test.result.second.wald.ER <- DisplaySecondStageTestResult(second.stage.logodds.meta.ER,second.stage.sigma.meta.ER)
  
  
  
  second.stage.logodds.meta.HER <-  log.odds.HER.Onco
  second.stage.sigma.meta.HER <-  sigma.log.odds.HER.Onco
  
  test.result.second.wald.HER <- DisplaySecondStageTestResult(second.stage.logodds.meta.HER,second.stage.sigma.meta.HER)
  
  
  
  second.stage.logodds.meta.Grade <-  log.odds.Grade.Onco
  second.stage.sigma.meta.Grade <-  sigma.log.odds.Grade.Onco
  
  test.result.second.wald.Grade <- DisplaySecondStageTestResult(second.stage.logodds.meta.Grade,second.stage.sigma.meta.Grade)
  
  
  
  
  
  
  
  
  
  
  heter.result <- list(data.frame(test.result.second.wald.PR,test.result.second.wald.ER,test.result.second.wald.HER,test.result.second.wald.Grade))
  
  
  save(heter.result,file=paste0("./known_SNPs/known_SNPs_analysis_G_revised/ERPosvesNeg/result/heter_result_",i1,".Rdata")) 
  
  
  
  
}
