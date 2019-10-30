#install_github("andrewhaoyu/bc2", args = c('--library="/home/zhangh24/R/x86_64-pc-linux-gnu-library/3.4"'))
###1 represent icog
###2 represent onco

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
  
  Heter.result.icog = EMmvpolySelfDesign(y.pheno.mis1,x.self.design = x.all.mis1[,1,drop=F],z.design=z.design,baselineonly = NULL,additive = x.all.mis1[,2:ncol(x.all.mis1)],pairwise.interaction = NULL,saturated = NULL,missingTumorIndicator = 888)
  z.standard <- Heter.result.icog[[12]]
  M <- nrow(z.standard)
  z.design.PR <- matrix(0,M,2)
  idx.PR.pos <- which(z.standard[,1]==1)
  z.design.PR[idx.PR.pos,1] <- 1
  idx.PR.neg <- which(z.standard[,1]==0)
  z.design.PR[,2] <- 1
  colnames(z.design.PR) <- c("PRcase","PRNeg")
  
  z.design.ER <- matrix(0,M,2)
  idx.ER.pos <- which(z.standard[,2]==1)
  z.design.ER[idx.ER.pos,1] <- 1
  idx.ER.neg <- which(z.standard[,2]==0)
  z.design.ER[,2] <- 1
  colnames(z.design.ER) <- c("ERcase","ERNeg")
  
  z.design.HER <- matrix(0,M,2)
  idx.HER.pos <- which(z.standard[,3]==1)
  z.design.HER[idx.HER.pos,1] <- 1
  idx.HER.neg <- which(z.standard[,3]==0)
  z.design.HER[,2] <- 1
  colnames(z.design.HER) <- c("HERcase","HERNeg")
  
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
  log.odds.icog <- Heter.result.icog[[1]][(M+1):(M+1+number.of.tumor)]
  nparm <- length(Heter.result.icog[[1]])  
  sigma.log.odds.icog <- Heter.result.icog[[2]][(M+1):(M+1+number.of.tumor),(M+1):(M+1+number.of.tumor)]
  beta.icog <- z.trans%*%log.odds.icog
  beta.sigma.icog <- z.trans%*%sigma.log.odds.icog%*%t(z.trans)
  loglikelihood.icog <- Heter.result.icog[[8]]
  
  Heter.result.PR.icog = EMmvpolySelfDesign(y.pheno.mis1,x.self.design = x.all.mis1[,1,drop=F],z.design=z.design.PR,baselineonly = NULL,additive = x.all.mis1[,2:ncol(x.all.mis1)],pairwise.interaction = NULL,saturated = NULL,missingTumorIndicator = 888)
  log.odds.PR.icog <-   Heter.result.PR.icog[[1]][(M+1):(M+2)]
  sigma.log.odds.PR.icog <- Heter.result.PR.icog[[2]][(M+1):(M+2),(M+1):(M+2)]
  beta.PR.icog <- z.design.PR%*%log.odds.PR.icog
  beta.sigma.PR.icog <- z.design.PR%*%sigma.log.odds.PR.icog%*%t(z.design.PR)
  
  p.PR <- c(2,1)
  beta.PR.icog <- beta.PR.icog[p.PR]  
  beta.sigma.PR.icog <- beta.sigma.PR.icog[p.PR,p.PR]  
  
  
  Heter.result.ER.icog = EMmvpolySelfDesign(y.pheno.mis1,x.self.design = x.all.mis1[,1,drop=F],z.design=z.design.ER,baselineonly = NULL,additive = x.all.mis1[,2:ncol(x.all.mis1)],pairwise.interaction = NULL,saturated = NULL,missingTumorIndicator = 888)
  log.odds.ER.icog <-   Heter.result.ER.icog[[1]][(M+1):(M+2)]
  sigma.log.odds.ER.icog <- Heter.result.ER.icog[[2]][(M+1):(M+2),(M+1):(M+2)]
  
  beta.ER.icog <- z.design.ER%*%log.odds.ER.icog
  beta.sigma.ER.icog <- z.design.ER%*%sigma.log.odds.ER.icog%*%t(z.design.ER)
  
  p.ER <- c(3,1)
  beta.ER.icog <- beta.ER.icog[p.ER]  
  beta.sigma.ER.icog <- beta.sigma.ER.icog[p.ER,p.ER]  
  
  
  Heter.result.HER.icog = EMmvpolySelfDesign(y.pheno.mis1,x.self.design = x.all.mis1[,1,drop=F],z.design=z.design.HER,baselineonly = NULL,additive = x.all.mis1[,2:ncol(x.all.mis1)],pairwise.interaction = NULL,saturated = NULL,missingTumorIndicator = 888)
  log.odds.HER.icog <-   Heter.result.HER.icog[[1]][(M+1):(M+2)]
  sigma.log.odds.HER.icog <- Heter.result.HER.icog[[2]][(M+1):(M+2),(M+1):(M+2)]
  beta.HER.icog <- z.design.HER%*%log.odds.HER.icog
  beta.sigma.HER.icog <- z.design.HER%*%sigma.log.odds.HER.icog%*%t(z.design.HER)
  
  p.HER <- c(12,1)
  beta.HER.icog <- beta.HER.icog[p.HER]  
  beta.sigma.HER.icog <- beta.sigma.HER.icog[p.HER,p.HER]  
  
  
  
  Heter.result.Grade.icog = EMmvpolySelfDesign(y.pheno.mis1,x.self.design = x.all.mis1[,1,drop=F],z.design=z.design.Grade,baselineonly = NULL,additive = x.all.mis1[,2:ncol(x.all.mis1)],pairwise.interaction = NULL,saturated = NULL,missingTumorIndicator = 888)
  log.odds.Grade.icog <-   Heter.result.Grade.icog[[1]][(M+1):(M+2)]
  sigma.log.odds.Grade.icog <- Heter.result.Grade.icog[[2]][(M+1):(M+2),(M+1):(M+2)]
  
  
  beta.grade.icog <- z.design.Grade%*%log.odds.Grade.icog
  beta.sigma.grade.icog <- z.design.Grade%*%sigma.log.odds.Grade.icog%*%t(z.design.Grade)
  
  p.grade <- c(1,8,16)
  beta.grade.icog <- beta.grade.icog[p.grade]  
  beta.sigma.grade.icog <- beta.sigma.grade.icog[p.grade,p.grade]  
  
  
  
  
  
  
  
  
  
  
  names1 = colnames(data1)[27:206]
  rm(data1)
  gc()
  
  
  
  #analysis for onco Array
  #data2 <- read.csv("./V10/onco_euro_v10_05242017.csv",header=T)
  data2 <- fread("./data/Onco_euro_v10_10232017.csv",header=T)
  data2 <- as.data.frame(data2)
  y.pheno.mis2 <- cbind(data2$Behaviour1,data2$PR_status1,data2$ER_status1,data2$HER2_status1,data2$Grade1)
  #y.pheno.mis2 <- cbind(data2$Behaviour1,data2$PR_status1,data2$ER_status1,data2$HER2_status1)
  colnames(y.pheno.mis2) = c("Behaviour","PR",
                             "ER","HER2","Grade")
  
  x.test.all.mis2 <- data2[,c(27:203)]
  x.test.all.mis2 <- x.test.all.mis2
  x.covar.mis2 <- data2[,c(5:14,204)]
  
  
  x.all.mis2 <- as.matrix(cbind(x.test.all.mis2[,i1],x.covar.mis2))
  colnames(x.all.mis2)[1] = "gene"
  
  
  
  
  Heter.result.onco = EMmvpolySelfDesign(y.pheno.mis2,x.self.design = x.all.mis2[,1,drop=F],z.design = z.design,baselineonly = NULL,additive = x.all.mis2[,2:ncol(x.all.mis2)],pairwise.interaction = NULL,saturated = NULL,missingTumorIndicator = 888)
  z.standard <- Heter.result.onco[[12]]
  M <- nrow(z.standard)
  number.of.tumor <- ncol(z.standard)
  log.odds.onco <- Heter.result.onco[[1]][(M+1):(M+1+number.of.tumor)]
  nparm <- length(Heter.result.onco[[1]])
  sigma.log.odds.onco <- Heter.result.onco[[2]][(M+1):(M+1+number.of.tumor),(M+1):(M+1+number.of.tumor)]
  beta.onco <- z.trans%*%log.odds.onco
  beta.sigma.onco <- z.trans%*%sigma.log.odds.onco%*%t(z.trans)
  loglikelihood.onco <- Heter.result.onco[[8]]
  
  
  Heter.result.PR.onco = EMmvpolySelfDesign(y.pheno.mis2,x.self.design = x.all.mis2[,1,drop=F],z.design=z.design.PR,baselineonly = NULL,additive = x.all.mis2[,2:ncol(x.all.mis2)],pairwise.interaction = NULL,saturated = NULL,missingTumorIndicator = 888)
  log.odds.PR.onco <-   Heter.result.PR.onco[[1]][(M+1):(M+2)]
  sigma.log.odds.PR.onco <- Heter.result.PR.onco[[2]][(M+1):(M+2),(M+1):(M+2)]
  beta.PR.onco <- z.design.PR%*%log.odds.PR.onco
  beta.sigma.PR.onco <- z.design.PR%*%sigma.log.odds.PR.onco%*%t(z.design.PR)
  
   p.PR <- c(2,1)
  beta.PR.onco <- beta.PR.onco[p.PR]  
  beta.sigma.PR.onco <- beta.sigma.PR.onco[p.PR,p.PR]  
  
  
  Heter.result.ER.onco = EMmvpolySelfDesign(y.pheno.mis2,x.self.design = x.all.mis2[,1,drop=F],z.design=z.design.ER,baselineonly = NULL,additive = x.all.mis2[,2:ncol(x.all.mis2)],pairwise.interaction = NULL,saturated = NULL,missingTumorIndicator = 888)
  log.odds.ER.onco <-   Heter.result.ER.onco[[1]][(M+1):(M+2)]
  sigma.log.odds.ER.onco <- Heter.result.ER.onco[[2]][(M+1):(M+2),(M+1):(M+2)]
  
  beta.ER.onco <- z.design.ER%*%log.odds.ER.onco
  beta.sigma.ER.onco <- z.design.ER%*%sigma.log.odds.ER.onco%*%t(z.design.ER)
  
  p.ER <- c(3,1)
  beta.ER.onco <- beta.ER.onco[p.ER]  
  beta.sigma.ER.onco <- beta.sigma.ER.onco[p.ER,p.ER]  
  
  
  Heter.result.HER.onco = EMmvpolySelfDesign(y.pheno.mis2,x.self.design = x.all.mis2[,1,drop=F],z.design=z.design.HER,baselineonly = NULL,additive = x.all.mis2[,2:ncol(x.all.mis2)],pairwise.interaction = NULL,saturated = NULL,missingTumorIndicator = 888)
  log.odds.HER.onco <-   Heter.result.HER.onco[[1]][(M+1):(M+2)]
  sigma.log.odds.HER.onco <- Heter.result.HER.onco[[2]][(M+1):(M+2),(M+1):(M+2)]
  beta.HER.onco <- z.design.HER%*%log.odds.HER.onco
  beta.sigma.HER.onco <- z.design.HER%*%sigma.log.odds.HER.onco%*%t(z.design.HER)
  
  p.HER <- c(12,1)
  beta.HER.onco <- beta.HER.onco[p.HER]  
  beta.sigma.HER.onco <- beta.sigma.HER.onco[p.HER,p.HER]  
  
  Heter.result.Grade.onco = EMmvpolySelfDesign(y.pheno.mis2,x.self.design = x.all.mis2[,1,drop=F],z.design=z.design.Grade,baselineonly = NULL,additive = x.all.mis2[,2:ncol(x.all.mis2)],pairwise.interaction = NULL,saturated = NULL,missingTumorIndicator = 888)
  log.odds.Grade.onco <-   Heter.result.Grade.onco[[1]][(M+1):(M+2)]
  sigma.log.odds.Grade.onco <- Heter.result.Grade.onco[[2]][(M+1):(M+2),(M+1):(M+2)]
  beta.grade.onco <- z.design.Grade%*%log.odds.Grade.onco
  beta.sigma.grade.onco <- z.design.Grade%*%sigma.log.odds.Grade.onco%*%t(z.design.Grade)
  
  p.grade <- c(1,8,16)
  beta.grade.onco <- beta.grade.onco[p.grade]  
  beta.sigma.grade.onco <- beta.sigma.grade.onco[p.grade,p.grade]  
  
  
  
  meta.result.PR <- LogoddsMetaAnalysis(log.odds.PR.icog,
                                        sigma.log.odds.PR.icog,
                                        log.odds.PR.onco,
                                        sigma.log.odds.PR.onco)
  
  second.stage.logodds.meta.PR <- meta.result.PR[[1]]
  second.stage.sigma.meta.PR <- meta.result.PR[[2]]
  
  meta.result.PR.first <- LogoddsMetaAnalysis(beta.PR.icog,
                                              beta.sigma.PR.icog,
                                              beta.PR.onco,
                                              beta.sigma.PR.onco)
  first.stage.logodds.meta.PR <- meta.result.PR.first[[1]]
  first.stage.sigma.meta.PR <- meta.result.PR.first[[2]]
                                              
  test.result.second.wald.PR <- DisplaySecondStageTestResult(second.stage.logodds.meta.PR,second.stage.sigma.meta.PR)
  test.result.first.wald.PR <- DisplayFirstStageTestResult(first.stage.logodds.meta.PR,
                                                           first.stage.sigma.meta.PR)
  
  
  meta.result.ER <- LogoddsMetaAnalysis(log.odds.ER.icog,
                                        sigma.log.odds.ER.icog,
                                        log.odds.ER.onco,
                                        sigma.log.odds.ER.onco)
  
  second.stage.logodds.meta.ER <- meta.result.ER[[1]]
  second.stage.sigma.meta.ER <- meta.result.ER[[2]]
  meta.result.ER.first <- LogoddsMetaAnalysis(beta.ER.icog,
                                              beta.sigma.ER.icog,
                                              beta.ER.onco,
                                              beta.sigma.ER.onco)
  first.stage.logodds.meta.ER <- meta.result.ER.first[[1]]
  first.stage.sigma.meta.ER <- meta.result.ER.first[[2]]
  
  test.result.second.wald.ER <- DisplaySecondStageTestResult(second.stage.logodds.meta.ER,second.stage.sigma.meta.ER)
  test.result.first.wald.ER <- DisplayFirstStageTestResult(first.stage.logodds.meta.ER,
                                                           first.stage.sigma.meta.ER)
  
  
  meta.result.HER <- LogoddsMetaAnalysis(log.odds.HER.icog,
                                         sigma.log.odds.HER.icog,
                                         log.odds.HER.onco,
                                         sigma.log.odds.HER.onco)
  
  second.stage.logodds.meta.HER <- meta.result.HER[[1]]
  second.stage.sigma.meta.HER <- meta.result.HER[[2]]
  
  meta.result.HER.first <- LogoddsMetaAnalysis(beta.HER.icog,
                                              beta.sigma.HER.icog,
                                              beta.HER.onco,
                                              beta.sigma.HER.onco)
  first.stage.logodds.meta.HER <- meta.result.HER.first[[1]]
  first.stage.sigma.meta.HER <- meta.result.HER.first[[2]]
  
  test.result.second.wald.HER <- DisplaySecondStageTestResult(second.stage.logodds.meta.HER,second.stage.sigma.meta.HER)
  test.result.first.wald.HER <- DisplayFirstStageTestResult(first.stage.logodds.meta.HER,
                                                           first.stage.sigma.meta.HER)
  
  meta.result.Grade <- LogoddsMetaAnalysis(log.odds.Grade.icog,
                                           sigma.log.odds.Grade.icog,
                                           log.odds.Grade.onco,
                                           sigma.log.odds.Grade.onco)
  
  second.stage.logodds.meta.Grade <- meta.result.Grade[[1]]
  second.stage.sigma.meta.Grade <- meta.result.Grade[[2]]
  
  
  first.stage.logodds.meta.Grade <- z.design.Grade%*%second.stage.logodds.meta.Grade
  first.stage.sigma.meta.Grade <- z.design.Grade%*%second.stage.sigma.meta.Grade%*%t(z.design.Grade)
  
  p.grade <- c(1,8,16)
  first.stage.logodds.meta.Grade <-  first.stage.logodds.meta.Grade[p.grade]  
  first.stage.sigma.meta.Grade<- first.stage.sigma.meta.Grade[p.grade,p.grade]  
  
  
  
  
  
  test.result.second.wald.Grade <- DisplaySecondStageTestResult(second.stage.logodds.meta.Grade,second.stage.sigma.meta.Grade)
  test.result.first.wald.Grade <- DisplayFirstStageTestResult(first.stage.logodds.meta.Grade,
                                                            first.stage.sigma.meta.Grade)
  
  
  
  
  
  
  
  
  
  heter.result <- list(data.frame(test.result.second.wald.PR,test.result.second.wald.ER,test.result.second.wald.HER,test.result.second.wald.Grade),
                       data.frame(test.result.first.wald.PR,test.result.first.wald.ER,test.result.first.wald.HER,test.result.first.wald.Grade))
  
  
  
  
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
  
  
  
  
  
  Heter.result.onco = EMmvpolySelfDesign(y.pheno.mis2,x.self.design = x.all.mis2[,1,drop=F],z.design = z.design,baselineonly = NULL,additive = x.all.mis2[,2:ncol(x.all.mis2)],pairwise.interaction = NULL,saturated = NULL,missingTumorIndicator = 888)
  z.standard <- Heter.result.onco[[12]]
  M <- nrow(z.standard)
  number.of.tumor <- ncol(z.standard)
  log.odds.onco <- Heter.result.onco[[1]][(M+1):(M+1+number.of.tumor)]
  nparm = length(Heter.result.onco[[1]])
  sigma.log.odds.onco <- Heter.result.onco[[2]][(M+1):(M+1+number.of.tumor),(M+1):(M+1+number.of.tumor)]
  beta.onco <- z.trans%*%log.odds.onco
  beta.sigma.onco <- z.trans%*%sigma.log.odds.onco%*%t(z.trans)
  loglikelihood.onco <- Heter.result.onco[[8]]
  
  
  M <- nrow(z.standard)
  z.design.PR <- matrix(0,M,2)
  idx.PR.pos <- which(z.standard[,1]==1)
  z.design.PR[idx.PR.pos,1] <- 1
  idx.PR.neg <- which(z.standard[,1]==0)
  z.design.PR[,2] <- 1
  colnames(z.design.PR) <- c("PRPos","PRNeg")
  
  z.design.ER <- matrix(0,M,2)
  idx.ER.pos <- which(z.standard[,2]==1)
  z.design.ER[idx.ER.pos,1] <- 1
  idx.ER.neg <- which(z.standard[,2]==0)
  z.design.ER[,2] <- 1
  colnames(z.design.ER) <- c("ERPos","ERNeg")
  
  z.design.HER <- matrix(0,M,2)
  idx.HER.pos <- which(z.standard[,3]==1)
  z.design.HER[idx.HER.pos,1] <- 1
  idx.HER.neg <- which(z.standard[,3]==0)
  z.design.HER[,2] <- 1
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
  
  
  
  Heter.result.PR.onco = EMmvpolySelfDesign(y.pheno.mis2,x.self.design = x.all.mis2[,1,drop=F],z.design=z.design.PR,baselineonly = NULL,additive = x.all.mis2[,2:ncol(x.all.mis2)],pairwise.interaction = NULL,saturated = NULL,missingTumorIndicator = 888)
  log.odds.PR.onco <-   Heter.result.PR.onco[[1]][(M+1):(M+2)]
  sigma.log.odds.PR.onco <- Heter.result.PR.onco[[2]][(M+1):(M+2),(M+1):(M+2)]
  beta.PR.onco <- z.design.PR%*%log.odds.PR.onco
  beta.sigma.PR.onco <- z.design.PR%*%sigma.log.odds.PR.onco%*%t(z.design.PR)
  
  p.PR <- c(2,1)
  beta.PR.onco <- beta.PR.onco[p.PR]  
  beta.sigma.PR.onco <- beta.sigma.PR.onco[p.PR,p.PR]  
  
  
  Heter.result.ER.onco = EMmvpolySelfDesign(y.pheno.mis2,x.self.design = x.all.mis2[,1,drop=F],z.design=z.design.ER,baselineonly = NULL,additive = x.all.mis2[,2:ncol(x.all.mis2)],pairwise.interaction = NULL,saturated = NULL,missingTumorIndicator = 888)
  log.odds.ER.onco <-   Heter.result.ER.onco[[1]][(M+1):(M+2)]
  sigma.log.odds.ER.onco <- Heter.result.ER.onco[[2]][(M+1):(M+2),(M+1):(M+2)]
  
  beta.ER.onco <- z.design.ER%*%log.odds.ER.onco
  beta.sigma.ER.onco <- z.design.ER%*%sigma.log.odds.ER.onco%*%t(z.design.ER)
  
  p.ER <- c(3,1)
  beta.ER.onco <- beta.ER.onco[p.ER]  
  beta.sigma.ER.onco <- beta.sigma.ER.onco[p.ER,p.ER]  
  
  
  Heter.result.HER.onco = EMmvpolySelfDesign(y.pheno.mis2,x.self.design = x.all.mis2[,1,drop=F],z.design=z.design.HER,baselineonly = NULL,additive = x.all.mis2[,2:ncol(x.all.mis2)],pairwise.interaction = NULL,saturated = NULL,missingTumorIndicator = 888)
  log.odds.HER.onco <-   Heter.result.HER.onco[[1]][(M+1):(M+2)]
  sigma.log.odds.HER.onco <- Heter.result.HER.onco[[2]][(M+1):(M+2),(M+1):(M+2)]
  beta.HER.onco <- z.design.HER%*%log.odds.HER.onco
  beta.sigma.HER.onco <- z.design.HER%*%sigma.log.odds.HER.onco%*%t(z.design.HER)
  
  p.HER <- c(12,1)
  beta.HER.onco <- beta.HER.onco[p.HER]  
  beta.sigma.HER.onco <- beta.sigma.HER.onco[p.HER,p.HER]  
  
  Heter.result.Grade.onco = EMmvpolySelfDesign(y.pheno.mis2,x.self.design = x.all.mis2[,1,drop=F],z.design=z.design.Grade,baselineonly = NULL,additive = x.all.mis2[,2:ncol(x.all.mis2)],pairwise.interaction = NULL,saturated = NULL,missingTumorIndicator = 888)
  log.odds.Grade.onco <-   Heter.result.Grade.onco[[1]][(M+1):(M+2)]
  sigma.log.odds.Grade.onco <- Heter.result.Grade.onco[[2]][(M+1):(M+2),(M+1):(M+2)]
  beta.grade.onco <- z.design.Grade%*%log.odds.Grade.onco
  beta.sigma.grade.onco <- z.design.Grade%*%sigma.log.odds.Grade.onco%*%t(z.design.Grade)
  
  p.grade <- c(1,8,16)
  beta.grade.onco <- beta.grade.onco[p.grade]  
  beta.sigma.grade.onco <- beta.sigma.grade.onco[p.grade,p.grade]  
  
  
  second.stage.logodds.meta.PR <- log.odds.PR.onco
  second.stage.sigma.meta.PR <- sigma.log.odds.PR.onco
  first.stage.logodds.meta.PR <- beta.PR.onco
  first.stage.sigma.meta.PR <- beta.sigma.PR.onco
  
  test.result.second.wald.PR <- DisplaySecondStageTestResult(second.stage.logodds.meta.PR,second.stage.sigma.meta.PR)
  test.result.first.wald.PR <- DisplayFirstStageTestResult(first.stage.logodds.meta.PR,
                                                           first.stage.sigma.meta.PR)
  
  
  second.stage.logodds.meta.ER <- log.odds.ER.onco
  second.stage.sigma.meta.ER <- sigma.log.odds.ER.onco
  first.stage.logodds.meta.ER <- beta.ER.onco
  first.stage.sigma.meta.ER <- beta.sigma.ER.onco
  
  test.result.second.wald.ER <- DisplaySecondStageTestResult(second.stage.logodds.meta.ER,second.stage.sigma.meta.ER)
  test.result.first.wald.ER <- DisplayFirstStageTestResult(first.stage.logodds.meta.ER,
                                                           first.stage.sigma.meta.ER)
  
  
  second.stage.logodds.meta.HER <- log.odds.HER.onco
  second.stage.sigma.meta.HER <- sigma.log.odds.HER.onco
  first.stage.logodds.meta.HER <- beta.HER.onco
  first.stage.sigma.meta.HER <- beta.sigma.HER.onco
  
  test.result.second.wald.HER <- DisplaySecondStageTestResult(second.stage.logodds.meta.HER,second.stage.sigma.meta.HER)
  test.result.first.wald.HER <- DisplayFirstStageTestResult(first.stage.logodds.meta.HER,
                                                           first.stage.sigma.meta.HER)
  
  
  second.stage.logodds.meta.Grade <- log.odds.Grade.onco
  second.stage.sigma.meta.Grade <- sigma.log.odds.Grade.onco
  
  
  first.stage.logodds.meta.Grade <- z.design.Grade%*%second.stage.logodds.meta.Grade
  first.stage.sigma.meta.Grade <- z.design.Grade%*%second.stage.sigma.meta.Grade%*%t(z.design.Grade)
  
  p.grade <- c(1,8,16)
  first.stage.logodds.meta.Grade <-  first.stage.logodds.meta.Grade[p.grade]  
  first.stage.sigma.meta.Grade<- first.stage.sigma.meta.Grade[p.grade,p.grade]  
  
  
  
  
  
  test.result.second.wald.Grade <- DisplaySecondStageTestResult(second.stage.logodds.meta.Grade,second.stage.sigma.meta.Grade)
  test.result.first.wald.Grade <- DisplayFirstStageTestResult(first.stage.logodds.meta.Grade,
                                                              first.stage.sigma.meta.Grade)
  
  
  
  
  
 
  
  
  
  
  
  
  heter.result <- list(data.frame(test.result.second.wald.PR,test.result.second.wald.ER,test.result.second.wald.HER,test.result.second.wald.Grade),
                       data.frame(test.result.first.wald.PR,test.result.first.wald.ER,test.result.first.wald.HER,test.result.first.wald.Grade))
  
  
  
  save(heter.result,file=paste0("./known_SNPs/known_SNPs_analysis_G_revised/ERPosvesNeg/result/heter_result_",i1,".Rdata")) 
  
  
  
  
}
