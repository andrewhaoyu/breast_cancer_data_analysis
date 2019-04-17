#install_github("andrewhaoyu/bc2",args = c('--library="/home/zhangh24/R/x86_64-pc-linux-gnu-library/3.4"'))
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
# load("./known_SNPs/known_SNPs_analysis_G_revised/additive_model/result/z.standard.Rdata")
load("./known_SNPs/known_SNPs_analysis_G_revised/additive_model_third_order/result/z.design.Rdata")
#save(z.design, file = "./known_SNPs/known_SNPs_analysis_G_revised/additive_model_third_order/result/z.design.Rdata")
# 
# idx <- which(z.standard[,1]==0&
#                z.standard[,2]==0&
#                z.standard[,3]==0)
# triple <- rep(1,23)
# triple[idx] <- 0
# z.standard <- cbind(z.standard,triple)
# colnames(z.standard) <- c("ER","PR","HER2","Grade","triple")
# z.design <- cbind(1,z.standard)
library(tidyverse)
#save(z.design,file="./known_SNPs/known_SNPs_analysis_G_revised/additive_model_third_order/result/z.design.Rdata")
#colnames(z.standard) <- c("ER","PR","HER2","Grade")

#third <- rep(0,nrow(z.standard))
# third[c(1,8,16)] <- 1
# z.new.design <- z.standard %>% mutate(
#   ER = 1-ER,
#   PR = 1-PR,
#   HER2 = 1-HER2
# ) %>% mutate(
#   third = third
# )

#z.new.design <- cbind(1,z.new.design)


#z.design <- z.new.design

# colnames(z.design) <- c("Baseline","ER","PR",
#                         "HER2",
#                         "Grade",
#                         "Third")


if(i1<=177){
  data1 <- as.data.frame(fread("./data/iCOGS_euro_v10_10232017.csv",header=T))
  data1 <- as.data.frame(data1)
  y.pheno.mis1 <- cbind(data1$Behaviour1,data1$ER_status1,data1$PR_status1,data1$HER2_status1,data1$Grade1)
  colnames(y.pheno.mis1) = c("Behavior","ER","PR","HER2","Grade")
  # Grade1.fake <- data1$Grade1
  # Grade1.fake[data1$Grade1==2|data1$Grade1==3] <- 1
  # Grade1.fake[data1$Grade1==1] <- 0
  #y.pheno.mis1 <- cbind(data1$Behaviour1,data1$PR_status1,data1$ER_status1,data1$HER2_status1,Grade1.fake)
  # y.pheno.mis1 <- cbind(data1$Behaviour1,data1$PR_status1,data1$ER_status1,data1$HER2_status1)
  
  x.test.all.mis1 <- data1[,c(27:203)]
  ###pc1-10 and age
  x.covar.mis1 <- data1[,c(5:14)]
  
  idx.control <- which(y.pheno.mis1[,1]==0)
  maf <- sum(x.test.all.mis1[idx.control,i1])/(2*length(idx.control))
  if(maf>=0.5){
    x.test.all.mis1[,i1] < 2 - x.test.all.mis1[,i1]
  }
  
  x.all.mis1 <- as.matrix(cbind(x.test.all.mis1[,i1],x.covar.mis1))
  colnames(x.all.mis1)[1] <- "gene"
  
  
  
  Heter.result.Icog = EMmvpolySelfDesign(y.pheno.mis1,x.self.design = x.all.mis1[,1,drop=F],z.design=z.design,baselineonly = NULL,additive = x.all.mis1[,2:ncol(x.all.mis1)],pairwise.interaction = NULL,saturated = NULL,missingTumorIndicator = 888)
  z.standard <- Heter.result.Icog[[12]]
  z.additive.design <- as.matrix(cbind(1,z.standard))
  M <- nrow(z.standard)
  number.of.tumor <- ncol(z.standard)+1
  log.odds.icog <- Heter.result.Icog[[1]][(M+1):(M+1+number.of.tumor)]
  nparm <- length(Heter.result.Icog[[1]])  
  sigma.log.odds.icog <- Heter.result.Icog[[2]][(M+1):(M+1+number.of.tumor),(M+1):(M+1+number.of.tumor)]
  # beta.icog <- z.trans%*%log.odds.icog
  # beta.sigma.icog <- z.trans%*%sigma.log.odds.icog%*%t(z.trans)
  loglikelihood.icog <- Heter.result.Icog[[8]]
  
  ###################Global test for association
  score.test.support.fixed.icog.ga <- 
    ScoreTestSupport(y.pheno.mis1,
                    baselineonly=NULL,
                    additive=x.all.mis1[,2:ncol(x.all.mis1)],
                    missingTumorIndicator = 888)
  score.test.fixed.icog.ga <- 
    ScoreTestSelfDesign(y=y.pheno.mis1,                                                           x=x.all.mis1[,1,drop=F],
                  z.design=z.design[,c(1,2,6)],
                  score.test.support= score.test.support.fixed.icog.ga,                           missingTumorIndicator=888)
  score.fixed.icog.ga <- score.test.fixed.icog.ga[[1]]
  infor.fixed.icog.ga <- score.test.fixed.icog.ga[[2]]
  
  
  score.test.support.random.icog.ga <- 
    ScoreTestSupportSelfDesign(y.pheno.mis1,                                                x.self.design=x.all.mis1[,1,drop=F],
            z.design=z.design[,c(1,2,6),drop=F],                                            additive=x.all.mis1[,2:ncol(x.all.mis1)],                                       missingTumorIndicator = 888)
  
  score.test.random.icog.ga <- 
    ScoreTestSelfDesign(y=y.pheno.mis1,                    
                        x=x.all.mis1[,1,drop=F],
                        z.design=z.design[,3:5],
                        score.test.support= score.test.support.random.icog.ga,                          missingTumorIndicator=888)
  score.random.icog.ga <- score.test.random.icog.ga[[1]]
  infor.random.icog.ga <- score.test.random.icog.ga[[2]]
  
  
   ###################Global test for heterogeneity
  score.test.support.fixed.icog <- ScoreTestSupportSelfDesign(y.pheno.mis1,
                                         x.self.design=x.all.mis1[,1,drop=F],
                                         z.design=z.design[,1,drop=F],
                                         baselineonly=NULL,
                                         additive=x.all.mis1[,2:ncol(x.all.mis1)],
                                         pairwise.interaction=NULL,
                                         saturated=NULL,
                                         missingTumorIndicator = 888)
  score.test.fixed.icog <- 
    ScoreTestSelfDesign(y=y.pheno.mis1,
                        x=x.all.mis1[,1,drop=F],
                        z.design=z.design[,2:ncol(z.design)],
                        score.test.support= score.test.support.fixed.icog,                              missingTumorIndicator=888)
  score.fixed.icog <-   score.test.fixed.icog[[1]]
  infor.fixed.icog <-   score.test.fixed.icog[[2]]
  score.test.support.random.icog <- 
    ScoreTestSupportSelfDesign(y.pheno.mis1,                                                x.self.design=x.all.mis1[,1,drop=F],
            z.design=z.design[,c(1,2,6),drop=F],                                            additive=x.all.mis1[,2:ncol(x.all.mis1)],                                       missingTumorIndicator = 888)
  
  score.test.random.icog <- 
    ScoreTestSelfDesign(y=y.pheno.mis1,
                        x=x.all.mis1[,1,drop=F],
                        z.design=z.design[,3:5],
                       score.test.support= score.test.support.random.icog,                           missingTumorIndicator=888)
  score.random.icog <- score.test.random.icog[[1]]
  infor.random.icog <- score.test.random.icog[[2]]
  
  
  

  data2 <- fread("./data/Onco_euro_v10_10232017.csv",header=T)
  data2 <- as.data.frame(data2)
  y.pheno.mis2 <- cbind(data2$Behaviour1,data2$ER_status1,data2$PR_status1,data2$HER2_status1,data2$Grade1)
  #y.pheno.mis2 <- cbind(data2$Behaviour1,data2$PR_status1,data2$ER_status1,data2$HER2_status1)
  colnames(y.pheno.mis2) = c("Behaviour","PR",
                             "ER","HER2","Grade")
  
  x.test.all.mis2 <- data2[,c(27:203)]
  x.covar.mis2 <- data2[,c(5:14)]
  
  if(maf>=0.5){
    x.test.all.mis2[,i1] <- 2-x.test.all.mis2[,i1]
  }
  
  x.all.mis2 <- as.matrix(cbind(x.test.all.mis2[,i1],x.covar.mis2))
  colnames(x.all.mis2)[1] = "gene"
  
  
  
  
  Heter.result.Onco = EMmvpolySelfDesign(y.pheno.mis2,x.self.design = x.all.mis2[,1,drop=F],z.design = z.design,baselineonly = NULL,additive = x.all.mis2[,2:ncol(x.all.mis2)],pairwise.interaction = NULL,saturated = NULL,missingTumorIndicator = 888)
  z.standard <- Heter.result.Onco[[12]]
  M <- nrow(z.standard)
  number.of.tumor <- ncol(z.standard)+1
  log.odds.onco <- Heter.result.Onco[[1]][(M+1):(M+1+number.of.tumor)]
  nparm <- length(Heter.result.Onco[[1]])
  sigma.log.odds.onco <- Heter.result.Onco[[2]][(M+1):(M+1+number.of.tumor),(M+1):(M+1+number.of.tumor)]
  # beta.onco <- z.trans%*%log.odds.onco
  # beta.sigma.onco <- z.trans%*%sigma.log.odds.onco%*%t(z.trans)
  loglikelihood.onco <- Heter.result.Onco[[8]]
  ###########global test for association
  score.test.support.fixed.onco.ga <- 
    ScoreTestSupport(y.pheno.mis2,
                    additive=x.all.mis2[,2:ncol(x.all.mis2)],
                    missingTumorIndicator = 888)
  score.test.fixed.onco.ga <- 
    ScoreTestSelfDesign(y=y.pheno.mis2,
                        x=x.all.mis2[,1,drop=F],
                        z.design=z.design[,c(1,2,6)],
                        score.test.support= score.test.support.fixed.onco.ga,                           missingTumorIndicator=888)
  score.fixed.onco.ga <-   score.test.fixed.onco.ga[[1]]
  infor.fixed.onco.ga <-   score.test.fixed.onco.ga[[2]]
  score.test.support.random.onco.ga <- 
    ScoreTestSupportSelfDesign(y.pheno.mis2,                                                                  x.self.design=x.all.mis2[,1,drop=F],
                              z.design=z.design[,c(1,2,6),drop=F],                                            additive=x.all.mis2[,2:ncol(x.all.mis2)],                                       missingTumorIndicator = 888)
  
  score.test.random.onco.ga <- 
    ScoreTestSelfDesign(y=y.pheno.mis2,                    
                        x=x.all.mis2[,1,drop=F],
                        z.design=z.design[,3:5],
                        score.test.support= score.test.support.random.onco.ga,                          missingTumorIndicator=888)
  score.random.onco.ga <-   score.test.random.onco.ga[[1]]
  infor.random.onco.ga <-   score.test.random.onco.ga[[2]]
  
  meta.result.score.fix.ga <- 
    ScoreMetaAnalysis(score.fixed.icog.ga,
                      infor.fixed.icog.ga,
                      score.fixed.onco.ga,
                      infor.fixed.onco.ga)
  score.fixed.meta.ga <- meta.result.score.fix.ga[[1]]
  infor.fixed.meta.ga <- meta.result.score.fix.ga[[2]]
  meta.result.score.random.ga <- 
    ScoreMetaAnalysis(score.random.icog.ga,
                      infor.random.icog.ga,
                      score.random.onco.ga,
                      infor.random.onco.ga)
  score.randomed.meta.ga <- meta.result.score.random.ga[[1]]
  infor.randomed.meta.ga <- meta.result.score.random.ga[[2]]
  
  test.result.second.mixed.ga <- 
    DisplayMixedScoreTestResult(score.fixed.meta.ga,                                                            infor.fixed.meta.ga,                                                            score.randomed.meta.ga,                                                         infor.randomed.meta.ga)  
  
  
  ###########global test for heterogneity  
  score.test.support.fixed.onco <- 
    ScoreTestSupportSelfDesign(y.pheno.mis2,
                              x.self.design=x.all.mis2[,1,drop=F],
                              z.design=z.design[,1,drop=F],
                              additive=x.all.mis2[,2:ncol(x.all.mis2)],
                              missingTumorIndicator = 888)
  score.test.fixed.onco <- 
    ScoreTestSelfDesign(y=y.pheno.mis2,
                        x=x.all.mis2[,1,drop=F],
                        z.design=z.design[,2:ncol(z.design)],
                        score.test.support= score.test.support.fixed.onco,                              missingTumorIndicator=888)
  score.fixed.onco <-   score.test.fixed.onco[[1]]
  infor.fixed.onco <-   score.test.fixed.onco[[2]]
  score.test.support.random.onco <- 
    ScoreTestSupportSelfDesign(y.pheno.mis2,                                                                  x.self.design=x.all.mis2[,1,drop=F],
                              z.design=z.design[,c(1,2,6),drop=F],                                            additive=x.all.mis2[,2:ncol(x.all.mis2)],                                       missingTumorIndicator = 888)
  
  score.test.random.onco <- 
    ScoreTestSelfDesign(y=y.pheno.mis2,                    
                        x=x.all.mis2[,1,drop=F],
                        z.design=z.design[,3:5],
                        score.test.support= score.test.support.random.onco,                             missingTumorIndicator=888)
  score.random.onco <- score.test.random.onco[[1]]
  infor.random.onco <- score.test.random.onco[[2]]
  meta.result.score.fix <- ScoreMetaAnalysis(score.fixed.icog,infor.fixed.icog,score.fixed.onco,infor.fixed.onco)
  score.fixed.meta <- meta.result.score.fix[[1]]
  infor.fixed.meta <- meta.result.score.fix[[2]]
  meta.result.score.random <- ScoreMetaAnalysis(score.random.icog,infor.random.icog,score.random.onco,infor.random.onco)
  score.randomed.meta <- meta.result.score.random[[1]]
  infor.randomed.meta <- meta.result.score.random[[2]]
  
  test.result.second.mixed <- DisplayMixedScoreTestResult(score.fixed.meta,                                                        infor.fixed.meta,                                        score.randomed.meta,                                        infor.randomed.meta)  
    
  
  
  
  meta.result <- LogoddsMetaAnalysis(log.odds.icog,
                                     sigma.log.odds.icog,
                                     log.odds.onco,
                                     sigma.log.odds.onco)
  
  second.stage.logodds.meta <- meta.result[[1]]
  second.stage.sigma.meta <- meta.result[[2]]
  
  
  
  test.result.second.wald <- DisplaySecondStageTestResult(second.stage.logodds.meta,second.stage.sigma.meta)
  
  
  heter.result <- list(test.result.second.wald,
                       test.result.second.mixed.ga,
                       test.result.second.mixed)
  
  save(heter.result,file=paste0("./known_SNPs/known_SNPs_analysis_G_revised/additive_model_third_order/result/heter_result_",i1,".Rdata"))
  
  
}else{
  data2 <- fread("./data/Onco_euro_v10_10232017.csv",header=T)
  data2 <- as.data.frame(data2)
  names2 = colnames(data2)
  y.pheno.mis2 <- cbind(data2$Behaviour1,data2$ER_status1,data2$PR_status1,data2$HER2_status1,data2$Grade1)
  #y.pheno.mis2 <- cbind(data2$Behaviour1,data2$PR_status1,data2$ER_status1,data2$HER2_status1)
  colnames(y.pheno.mis2) = c("Behaviour","ER",
                             "PR","HER2","Grade")
  idxi1 = which(names2=="rs554219")
  
  x.test.all.mis2 <- data2
  x.covar.mis2 <- data2[,c(5:14)]
  
  
  x.all.mis2 <- as.matrix(cbind(x.test.all.mis2[,idxi1],x.covar.mis2))
  colnames(x.all.mis2)[1] = "gene"
  
  
  
  
  
  
  Heter.result.Onco = EMmvpolySelfDesign(
    y.pheno.mis2,
    x.self.design = x.all.mis2[, 1, drop = F],
    z.design = z.design,
    baselineonly = NULL,
    additive = x.all.mis2[, 2:ncol(x.all.mis2)],
    pairwise.interaction = NULL,
    saturated = NULL,
    missingTumorIndicator = 888
  )
  z.standard <- Heter.result.Onco[[12]]
  M <- nrow(z.standard)
  number.of.tumor <- ncol(z.standard) + 1
  log.odds.onco <-
    Heter.result.Onco[[1]][(M + 1):(M + 1 + number.of.tumor)]
  nparm = length(Heter.result.Onco[[1]])
  sigma.log.odds.onco <-
    Heter.result.Onco[[2]][(M + 1):(M + 1 + number.of.tumor), (M + 1):(M +
                                                                         1 + number.of.tumor)]
  # beta.onco <- z.trans%*%log.odds.onco
  # beta.sigma.onco <- z.trans%*%sigma.log.odds.onco%*%t(z.trans)
  loglikelihood.onco <- Heter.result.Onco[[8]]
  
  ###########global test for association
  score.test.support.fixed.onco.ga <- 
    ScoreTestSupport(y.pheno.mis2,
                     additive=x.all.mis2[,2:ncol(x.all.mis2)],
                     missingTumorIndicator = 888)
  score.test.fixed.onco.ga <- 
    ScoreTestSelfDesign(y=y.pheno.mis2,
                        x=x.all.mis2[,1,drop=F],
                        z.design=z.design[,c(1,2,6)],
                        score.test.support= score.test.support.fixed.onco.ga,                           missingTumorIndicator=888)
  score.fixed.onco.ga <-   score.test.fixed.onco.ga[[1]]
  infor.fixed.onco.ga <-   score.test.fixed.onco.ga[[2]]
  score.test.support.random.onco.ga <- 
    ScoreTestSupportSelfDesign(y.pheno.mis2,                                                                  x.self.design=x.all.mis2[,1,drop=F],
                               z.design=z.design[,c(1,2,6),drop=F],                                            additive=x.all.mis2[,2:ncol(x.all.mis2)],                                       missingTumorIndicator = 888)
  
  score.test.random.onco.ga <- 
    ScoreTestSelfDesign(y=y.pheno.mis2,                    
                        x=x.all.mis2[,1,drop=F],
                        z.design=z.design[,3:5],
                        score.test.support= score.test.support.random.onco.ga,                          missingTumorIndicator=888)
  score.random.onco.ga <-   score.test.random.onco.ga[[1]]
  infor.random.onco.ga <-   score.test.random.onco.ga[[2]]
  
  score.fixed.meta.ga <-   score.fixed.onco.ga
  infor.fixed.meta.ga <- infor.fixed.onco.ga
  # meta.result.score.random <- ScoreMetaAnalysis(score.random.icog,infor.random.icog,score.random.onco,infor.random.onco)
  score.randomed.meta.ga <- score.random.onco.ga
  infor.randomed.meta.ga <- infor.random.onco.ga
  
  test.result.second.mixed.ga <- 
    DisplayMixedScoreTestResult(score.fixed.meta.ga,                                                            infor.fixed.meta.ga,                                                            score.randomed.meta.ga,                                                         infor.randomed.meta.ga)  
  
    ##############global test for heterogeneity
  score.test.support.fixed.onco <- 
    ScoreTestSupportSelfDesign(y.pheno.mis2,
                               x.self.design=x.all.mis2[,1,drop=F],
                               z.design=z.design[,1,drop=F],
                               additive=x.all.mis2[,2:ncol(x.all.mis2)],
                               missingTumorIndicator = 888)
  score.test.fixed.onco <- 
    ScoreTestSelfDesign(y=y.pheno.mis2,
                        x=x.all.mis2[,1,drop=F],
                        z.design=z.design[,2:ncol(z.design)],
                        score.test.support= score.test.support.fixed.onco,                              missingTumorIndicator=888)
  score.fixed.onco <-   score.test.fixed.onco[[1]]
  infor.fixed.onco <-   score.test.fixed.onco[[2]]
  score.test.support.random.onco <- 
    ScoreTestSupportSelfDesign(y.pheno.mis2,                                                                   x.self.design=x.all.mis2[,1,drop=F],
                               z.design=z.design[,c(1,2,6),drop=F],                                            additive=x.all.mis2[,2:ncol(x.all.mis2)],                                       missingTumorIndicator = 888)
  
  score.test.random.onco <- 
    ScoreTestSelfDesign(y=y.pheno.mis2,
                        x=x.all.mis2[,1,drop=F],
                        z.design=z.design[,3:5],
                        score.test.support= score.test.support.random.onco,                             missingTumorIndicator=888)
  score.random.onco <- score.test.random.onco[[1]]
  infor.random.onco <- score.test.random.onco[[2]]
 # meta.result.score.fix <- ScoreMetaAnalysis(score.fixed.icog,infor.fixed.icog,score.fixed.onco,infor.fixed.onco)
  score.fixed.meta <-   score.fixed.onco 
  infor.fixed.meta <- infor.fixed.onco
 # meta.result.score.random <- ScoreMetaAnalysis(score.random.icog,infor.random.icog,score.random.onco,infor.random.onco)
  score.randomed.meta <- score.random.onco
  infor.randomed.meta <- infor.random.onco
  
  test.result.second.mixed <- DisplayMixedScoreTestResult(score.fixed.meta,                                                        infor.fixed.meta,                                        score.randomed.meta,                                        infor.randomed.meta)  
  
  
  second.stage.logodds.meta <-log.odds.onco
  second.stage.sigma.meta <- sigma.log.odds.onco
  
  
  
  
  
  test.result.second.wald <- DisplaySecondStageTestResult(second.stage.logodds.meta,second.stage.sigma.meta)
  
  heter.result <- list(test.result.second.wald,
                       test.result.second.mixed.ga,
                       test.result.second.mixed)
  
  save(heter.result,file=paste0("./known_SNPs/known_SNPs_analysis_G_revised/additive_model_third_order/result/heter_result_",i1,".Rdata"))
  
 
  
}



