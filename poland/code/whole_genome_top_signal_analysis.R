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
library(nnet)
Generatesubtypes<- function(ER,PR,HER2,Grade){
  n <- length(ER)
  subtypes <- rep("control",n)
  temp = 1
  for(i in 0:1){
    for(j in 0:1){
      for(k in 0:1){
        for(l in 1:3){
          
          idx <- which(ER==i&PR==j&HER2==k&Grade==l)
          if(length(idx)!=0){
            subtypes[idx] <- temp
            temp = temp+1  
          }
          
        }
      }
    }
  }
  
  subtypes <- factor(subtypes,levels=c("control",
                                       c(1:temp)))
  sum <- table(subtypes)
  idx.cat <- which(sum<=10)
  idx.remove <- which((subtypes%in%(unique(idx.cat)-1))==T)
  subtypes <- subtypes[-idx.remove]
  return(list(subtypes,idx.remove))
}



if(i1<=2){
  top_signal <- read.csv("/spin1/users/zhangh24/breast_cancer_data_analysis/poland/result/whole_genome/top_signal_in_poland.csv",header=T)
  data2 <- fread("./data/Onco_euro_v10_10232017.csv",header=T)
  data2 <- as.data.frame(data2)
  idx <- which(data2$StudyCountry=="Poland")
  data2 <- data2[idx,]
  # table(data2$Behaviour1)
  # table(data2$ER_status1)
  # table(data2$PR_status1)
  # table(data2$HER2_status1)
  # table(data2$Grade1)
  y.pheno.mis2 <- cbind(data2$Behaviour1,data2$ER_status1,data2$PR_status1,data2$HER2_status1,data2$Grade1)
  
  #y.pheno.mis2 <- cbind(data2$Behaviour1,data2$PR_status1,data2$ER_status1,data2$HER2_status1)
  colnames(y.pheno.mis2) = c("Behaviour","ER",
                             "PR","HER2","Grade")
  
  x.test.all.mis2 <- top_signal[idx,]
  x.covar.mis2 <- data2[,c(5:8,204)]
  
  
  x.all.mis2 <- as.matrix(cbind(x.test.all.mis2[,i1],x.covar.mis2))
  colnames(x.all.mis2)[1] = "gene"
  
  
  
  
  Heter.result.Onco = EMmvpoly(y.pheno.mis2,baselineonly = NULL,additive = x.all.mis2,pairwise.interaction = NULL,saturated = NULL,missingTumorIndicator = 888)
  z.standard <- Heter.result.Onco[[12]]
  M <- nrow(z.standard)
  number.of.tumor <- ncol(z.standard)
  log.odds.onco <- Heter.result.Onco[[1]][(M+1):(M+1+number.of.tumor)]
  nparm <- length(Heter.result.Onco)
  sigma.log.odds.onco <- Heter.result.Onco[[2]][(M+1):(M+1+number.of.tumor),(M+1):(M+1+number.of.tumor)]
  # beta.onco <- z.additive.design%*%log.odds.onco
  # beta.sigma.onco <- z.additive.design%*%sigma.log.odds.onco%*%t(z.additive.design)
  # loglikelihood.onco <- Heter.result.Onco[[8]]
  #rm(Heter.result.Onco)  
  
  score.test.support.onco <- ScoreTestSupport(
    y.pheno.mis2,
    baselineonly = NULL,
    additive = x.all.mis2[,2:ncol(x.all.mis2)],
    pairwise.interaction = NULL,
    saturated = NULL,
    missingTumorIndicator = 888
  )
  score.test.onco<- ScoreTest(y=y.pheno.mis2,
                              x=x.all.mis2[,1,drop=F],
                              second.stage.structure="additive",
                              score.test.support=score.test.support.onco,
                              missingTumorIndicator=888)
  z.design.score.baseline <- matrix(rep(1,M),ncol=1)
  z.design.score.casecase <-z.standard
  z.design.score.baseline.ER <- cbind(z.design.score.baseline,z.standard[,1])
  z.design.score.casecase.ER <- z.standard[,2:ncol(z.standard)]
  
  score.onco <- score.test.onco[[1]]
  infor.onco <- score.test.onco[[2]]
  score.test.support.onco.baseline <- score.test.support.onco
  #rm(score.test.support.onco)
  rm(score.test.onco)
  
  score.test.onco.baseline<- ScoreTestSelfDesign(y=y.pheno.mis2,
                                                 x=x.all.mis2[,1,drop=F],
                                                 z.design=z.design.score.baseline,
                                                 score.test.support=score.test.support.onco.baseline,
                                                 missingTumorIndicator=888)
  
  score.onco.baseline <- score.test.onco.baseline[[1]]
  infor.onco.baseline <- score.test.onco.baseline[[2]]
  rm(score.test.support.onco.baseline)
  rm(score.test.onco.baseline)
  gc()
  score.test.support.onco.casecase <- ScoreTestSupport(
    y.pheno.mis2,
    baselineonly = x.all.mis2[,1,drop=F],
    additive = x.all.mis2[,2:ncol(x.all.mis2)],
    pairwise.interaction = NULL,
    saturated = NULL,
    missingTumorIndicator = 888
  )
  score.test.onco.casecase<- ScoreTestSelfDesign(y=y.pheno.mis2,
                                                 x=x.all.mis2[,1,drop=F],
                                                 z.design=z.design.score.casecase,
                                                 score.test.support=score.test.support.onco.casecase,
                                                 missingTumorIndicator=888)
  
  score.onco.casecase <- score.test.onco.casecase[[1]]
  infor.onco.casecase <- score.test.onco.casecase[[2]]
  z.design.score.baseline.heterER <- z.standard[,1,drop=F]
  
  
  score.onco.baseline.heterER <-   score.onco.casecase[1]
  infor.onco.baseline.heterER <- infor.onco.casecase[1,1]
  
  
  
  
  rm(score.test.support.onco.casecase) 
  
  
  
  
  
  
  
  score.test.onco.baseline.ER<- ScoreTestSelfDesign(y=y.pheno.mis2,
                                                    x=x.all.mis2[,1,drop=F],
                                                    z.design=z.design.score.baseline.ER,
                                                    score.test.support=score.test.support.onco,
                                                    missingTumorIndicator=888)
  
  score.onco.baseline.ER <- score.test.onco.baseline.ER[[1]]
  infor.onco.baseline.ER <- score.test.onco.baseline.ER[[2]]
  
  score.test.support.onco.casecase.ER <- ScoreTestSupportSelfDesign(
    y.pheno.mis2,
    x.self.design  = x.all.mis2[,1,drop=F],
    z.design = z.design.score.baseline.ER,
    additive = x.all.mis2[,2:ncol(x.all.mis2)],
    pairwise.interaction = NULL,
    saturated = NULL,
    missingTumorIndicator = 888
  )
  
  score.test.onco.casecase.ER<- ScoreTestSelfDesign(y=y.pheno.mis2,
                                                    x=x.all.mis2[,1,drop=F],
                                                    z.design=z.design.score.casecase.ER,
                                                    score.test.support=score.test.support.onco.casecase.ER,
                                                    missingTumorIndicator=888)
  
  score.onco.casecase.ER <- score.test.onco.casecase.ER[[1]]
  infor.onco.casecase.ER <- score.test.onco.casecase.ER[[2]]
  
  rm(score.test.support.onco.casecase.ER)
  rm(score.test.onco.casecase.ER)
  gc()
  
  
  
  
  second.stage.logodds.meta <- log.odds.onco
  second.stage.sigma.meta <- sigma.log.odds.onco
  
  
  
  test.result.second.wald <- DisplaySecondStageTestResult(second.stage.logodds.meta,second.stage.sigma.meta)
  
  
  # beta.meta <- z.additive.design%*%second.stage.logodds.meta
  # beta.sigma.meta <- z.additive.design%*%second.stage.sigma.meta%*%t(z.additive.design)
  # 
  # test.result.first.wald <- DisplayFirstStageTestResult(beta.meta,beta.sigma.meta)
  
  score.meta <- score.onco
  infor.meta <- infor.onco
  
  test.result.second.score <- DisplayFixedScoreTestResult(score.meta,infor.meta)
  
  score.meta.baseline <- score.onco.baseline
  infor.meta.baseline <- infor.onco.baseline
  
  score.meta.casecase <- score.onco.casecase
  infor.meta.casecase <- infor.onco.casecase
  
  score.meta.baseline.heterER <- score.onco.baseline.heterER  
  infor.meta.baseline.heterER <- infor.onco.baseline.heterER
  
  
  score.meta.baseline.ER <-   score.onco.baseline.ER
  infor.meta.baseline.ER <- infor.onco.baseline.ER
  
  score.meta.casecase.ER <-   score.onco.casecase.ER
  infor.meta.casecase.ER <- infor.onco.casecase.ER
  
  
  
  
  test.result.second.mixed <- DisplayMixedScoreTestResult(score.meta.baseline,
                                                          infor.meta.baseline,
                                                          score.meta.casecase,
                                                          infor.meta.casecase)  
  test.result.second.mixed <- data.frame(t(test.result.second.mixed))
  
  
  
  ###########DisplayMixedScoreTestResult is set up for one fixed effect,
  ###########we need to adjust for the heterogeneity test
  test.result.second.mixed.ER <- DisplayMixedScoreTestResult(score.meta.baseline.ER,
                                                             infor.meta.baseline.ER,
                                                             score.meta.casecase.ER,
                                                             infor.meta.casecase.ER)  
  test.result.second.mixed.ER <- data.frame(t(test.result.second.mixed.ER))
  
  test.result.second.mixed.ER.2 <- DisplayMixedScoreTestResult(score.meta.baseline.heterER,
                                                               infor.meta.baseline.heterER,
                                                               score.meta.casecase.ER,
                                                               infor.meta.casecase.ER)  
  test.result.second.mixed.ER.2 <- data.frame(t(test.result.second.mixed.ER.2))
  
  ###Both ER and RANDOM EFFECT should be counted into heterogeneity
  test.result.second.mixed.ER[1,2] <-  test.result.second.mixed.ER.2[1,1]
  
  
  colnames(test.result.second.mixed) <- c("mixed model global test for association(baseline fixed)","mixed model global test for heterogeneity(baseline fixed)")
  colnames(test.result.second.mixed.ER) <- c("mixed model global test for association(baseline and ER fixed)","mixed model global test for heterogeneity(baseline and ER fixed)")
  
  # loglikelihood <- loglikelihood.icog+loglikelihood.onco
  # AIC <- 2*length(Heter.result.Onco[[1]])-2*loglikelihood
  # 
  
  x <- as.matrix(cbind(x.test.all.mis2[,i1],x.covar.mis2))
  model.standard <- glm(y.pheno.mis2[,1]~x,family = binomial(link='logit'))
  
  summary(model.standard)$coefficients[2,4]
  
  heter.result <- data.frame(test.result.second.wald,test.result.second.score, test.result.second.mixed,test.result.second.mixed.ER,summary(model.standard)$coefficients[2,4])
  colnames(heter.result)[18] <- "standard_p"
  data2 <- fread("./data/Onco_euro_v10_10232017.csv",header=T)
  data2 <- as.data.frame(data2)
  idx <- which(data2$StudyCountry=="Poland")
  data2 <- data2[idx,]
  
  idx <- which(data2$ER_status1==888|data2$PR_status1==888|
                 data2$HER2_status1==888|data2$Grade1==888)
  data2 <- data2[-idx,]
  y.pheno.mis2 <- cbind(data2$Behaviour1,data2$ER_status1,data2$PR_status1,data2$HER2_status1,data2$Grade1)
  
  #y.pheno.mis2 <- cbind(data2$Behaviour1,data2$PR_status1,data2$ER_status1,data2$HER2_status1)
  colnames(y.pheno.mis2) = c("Behaviour","ER",
                             "PR","HER2","Grade")
  
  x.test.all.mis2 <- data2[,c(27:203,205)]
  x.test.all.mis2 <- x.test.all.mis2
  x.covar.mis2 <- data2[,c(5:8,204)]
  snpvalue <- x.test.all.mis2[,i1]
  
  
  temp <-  Generatesubtypes(data2$ER_status1,data2$PR_status1,data2$HER2_status1,data2$Grade1)
  subtypes <- temp[[1]]
  idx.remove <- temp[[2]]
  x <- cbind(snpvalue,as.matrix(x.covar.mis2))[-idx.remove,]
  
  
  poly.model <- multinom(subtypes~x)
  poly.model.coef <- coef(poly.model)
  M <- nrow(poly.model.coef)
  p.covariate <- ncol(poly.model.coef)
  snp.cov <- vcov(poly.model)[2+p.covariate*(0:(M-1)),2+p.covariate*(0:(M-1))]
  snp.coef <- poly.model.coef[,2]
  p.poly <- DisplaySecondStageTestResult(snp.coef,snp.cov)[35]
  heter.result <- cbind(heter.result,p.poly)
  colnames(heter.result)[19] <- "polytomous_p"
  
  save(heter.result,file=paste0("./result/whole_genome/heter_result_",i1,".Rdata"))
  
  
  
}