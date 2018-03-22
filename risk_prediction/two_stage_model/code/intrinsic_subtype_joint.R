#install_github("andrewhaoyu/bcutility",args = c('--library="/home/zhangh24/R/x86_64-pc-linux-gnu-library/3.4"'))
rm(list=ls())
arg <- commandArgs(trailingOnly=T)
i1 <- as.numeric(arg[[1]])
library(bcutility)
library(bc2)
library(rstan)
z.design <- matrix(c(
  c(0,1,1,1,0,0,0,0,1,1,1,0,0,0,0,0,1,1,0,0,0,0,0),
  c(0,0,0,0,0,1,1,0,0,0,0,0,1,1,1,0,0,0,0,0,1,1,1),
  c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0),
  c(0,0,0,0,1,0,0,0,0,0,0,1,0,0,0,0,0,0,0,1,0,0,0),
  c(1,0,0,0,0,0,0,1,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0)
),ncol=5)




colnames(z.design) <- c("Luminial A","Luminal B",
                        "Luminal B HER2Neg",
                        "HER2 Enriched",
                        "Triple Negative")







library(data.table)
icog.data <- as.data.frame(fread("/spin1/users/zhangh24/breast_cancer_data_analysis/data/sig_snps_icog.csv",header=T))
onco.data <- as.data.frame(fread("/spin1/users/zhangh24/breast_cancer_data_analysis/data/sig_snps_onco.csv",header=T))
library(tidyverse)
y.pheno.mis1 <- select(icog.data,Behaviour1,ER_status1,PR_status1,HER2_status1,Grade1)

subtypes.icog <- GenerateIntrinsicmis(y.pheno.mis1[,2],y.pheno.mis1[,3],
                                      y.pheno.mis1[,4],y.pheno.mis1[,5])

table(subtypes.icog)+table(subtypes.onco)

x.covar1 <- select(icog.data,5:14)
x.snp.all1 <- select(icog.data,26:230)
colnames(y.pheno.mis1)

icog.test.id <- Generatetestid(subtypes.icog)

y.pheno.mis2 <- select(onco.data,Behaviour1,ER_status1,PR_status1,HER2_status1,Grade1)
subtypes.onco <- GenerateIntrinsicmis(y.pheno.mis2[,2],
                                      y.pheno.mis2[,3],
                                      y.pheno.mis2[,4],
                                      y.pheno.mis2[,5])
x.covar2 <- select(onco.data,5:14)
x.snp.all2 <- select(onco.data,26:230)
colnames(y.pheno.mis2)
subtypes.onco <- GenerateIntrinsicmis(y.pheno.mis2[,2],y.pheno.mis2[,3],
                                      y.pheno.mis2[,4],y.pheno.mis2[,5])
onco.test.id <- Generatetestid(subtypes.onco)



M <- 5
log.odds.standard.result <- rep(0,6)
log.odds.intrinsic.result <- matrix(0,6,M)
sigma.intrinsic.result <- matrix(0,6,M^2)
log.odds.dic.result <-  matrix(0,6,M)
log.odds.intrinsic.eb <-  matrix(0,6,M)
log.odds.intrinsic.la <-  matrix(0,6,M)

heter.variance.estiamte <- rep(0,6)


for(i in 1:6){
  
  ##########split the data
  idx.test.case <-  icog.test.id[[1]][[i]]
  idx.test.control <- icog.test.id[[2]][[i]]
  idx.test1 <- c(idx.test.case,idx.test.control)
  y.pheno.mis1.train <- y.pheno.mis1[-idx.test1,]
  x.covar.train1 <- x.covar1[-idx.test1,]
  x.snp.all.train1 <- x.snp.all1[-idx.test1,]
  idx.test.case <-  onco.test.id[[1]][[i]]
  idx.test.control <- onco.test.id[[2]][[i]]
  idx.test2 <- c(idx.test.case,idx.test.control)
  y.pheno.mis2.train <- y.pheno.mis2[-idx.test2,]
  x.covar.train2 <- x.covar2[-idx.test2,]
  x.snp.all.train2 <- x.snp.all2[-idx.test2,]
  
  ############### standard analysis
  model1 <- glm(y.pheno.mis1.train[,1]~as.matrix(x.snp.all.train1)+as.matrix(x.covar.train1),family = binomial(link = 'logit'))
  log.odds.icog <-as.vector(model1$coefficients)[c(2:206)]
  sigma.icog <- as.matrix(vcov(model1)[c(2:206),c(2:206)])
  
  model2 <- glm(y.pheno.mis2.train[,1]~as.matrix(x.snp.all.train2)+as.matrix(x.covar.train2),family = binomial(link = 'logit'))
  log.odds.onco <-as.vector(model2$coefficients)[c(2:206)]
  sigma.onco <- vcov(model2)[c(2:206),c(2:206)]
  
  meta.result <- LogoddsMetaAnalysis(log.odds.icog,
                                     sigma.icog,
                                     log.odds.onco,
                                     sigma.onco)
  log.odds.meta <- meta.result[[1]][i1]
  log.odds.standard.result[i] <- log.odds.meta
  

  
  
   
  
  #############intrinsic subtypes analysis
  Heter.result.Icog = EMmvpolySelfDesign(y.pheno.mis1.train,x.self.design = x.snp.all.train1[,i1:2],z.design=z.design,baselineonly = NULL,additive = x.covar.train1,pairwise.interaction = NULL,saturated = NULL,missingTumorIndicator = 888)
  z.standard <- Heter.result.Icog[[12]]
  M <- nrow(z.standard)
  number.of.tumor <- ncol(z.standard)
  log.odds.icog <- Heter.result.Icog[[1]][(M+1):(M+1+number.of.tumor)]
  nparm <- length(Heter.result.Icog[[1]])  
  sigma.log.odds.icog <- Heter.result.Icog[[2]][(M+1):(M+1+number.of.tumor),(M+1):(M+1+number.of.tumor)]
  
  Heter.result.Onco = EMmvpolySelfDesign(y.pheno.mis2.train,x.self.design = x.snp.all.train2[,i1],z.design = z.design,baselineonly = NULL,additive = x.covar.train2,pairwise.interaction = NULL,saturated = NULL,missingTumorIndicator = 888)
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
  p.heter <- GlobalTestForHeter(meta.result[[1]],
                                meta.result[[2]],
                                self.design = T)
  
  
  
  log.odds.intrinsic.result[i,] <- meta.result[[1]] 
  sigma.intrinsic.result[i,] <- as.vector(meta.result[[2]])
  log.odds.dic.result[i,] <- dicestiamte(meta.result[[1]],
                                     meta.result[[2]],
                                     log.odds.standard.result[i],
                                     p.heter,
                                     0.05)
  
  heter.variance.estiamte [i] <- 
    HeterVarianceEstimate(meta.result[[1]],
                          meta.result[[2]])
  
  log.odds.intrinsic.eb[i,] <- ebestimate(
    meta.result[[1]],
    meta.result[[2]],
    log.odds.standard.result[i],
    heter.variance.estiamte [i] 
  )
  
  log.odds.intrinsic.la[i,] <- eblaplace(
   as.vector(meta.result[[1]]),
    meta.result[[2]],
    log.odds.standard.result[i],
    heter.variance.estiamte [i]
  )
  
}

all.model.result <- list(log.odds.standard.result
                         =log.odds.standard.result,
                         log.odds.intrinsic.result
                         =log.odds.intrinsic.result,
                         sigma.intrinsic.result
                         =sigma.intrinsic.result,
                         log.odds.dic.result
                         =log.odds.dic.result,
                         log.odds.intrinsic.eb
                         =log.odds.intrinsic.eb,
                         log.odds.intrinsic.la
                         =log.odds.intrinsic.la,
                         heter.variance.estiamte
                         =heter.variance.estiamte)

save(all.model.result,file=paste0("/spin1/users/zhangh24/breast_cancer_data_analysis/risk_prediction/two_stage_model/result/all.model.result",i1,".Rdata"))












# #############random sample 200 cases & 200 controls from icog
# #############random sample 500 cases & 500 controls from onco
# 
# n.test.control.icog <- 200
# n.test.cases.icog <- 200
# n.test.control.onco <- 500
# n.test.cases.onco <- 500
# set.seed(1)
# idx.test.control.icog <- idx.control1[sample(length(idx.control1),n.test.control.icog)]
# idx.test.triple.icog <- idx.triple1[sample(length(idx.triple1),n.test.cases.icog)]
# idx.test1 <- c(idx.test.control.icog,idx.test.triple.icog)
# y.pheno.mis1.test <- y.pheno.mis1[idx.test1,]
# x.covar.test1 <- x.covar1[idx.test1,]
# x.snp.all.test1 <- x.snp.all1[idx.test1,]
# 
# 
# y.pheno.mis1.train <- y.pheno.mis1[-idx.test1,]
# x.covar.train1 <- x.covar1[-idx.test1,]
# x.snp.all.train1 <- x.snp.all1[-idx.test1,]
# 
# 
# 
# 
# idx.test.control.onco <- idx.control2[sample(length(idx.control2),n.test.control.onco)]
# idx.test.triple.onco <- idx.triple2[sample(length(idx.triple2),n.test.cases.onco)]
# idx.test2 <- c(idx.test.control.onco,idx.test.triple.onco)
# y.pheno.mis2.test <- y.pheno.mis2[idx.test2,]
# x.covar.test2 <- x.covar2[idx.test2,]
# x.snp.all.test2 <- x.snp.all2[idx.test2,]
# 
# 
# y.test <- rbind(y.pheno.mis1.test,y.pheno.mis2.test)
# x.snp.all.test <- rbind(as.matrix(x.snp.all.test1),as.matrix(x.snp.all.test2))
# 
# y.pheno.mis2.train <- y.pheno.mis2[-idx.test2,]
# x.covar.train2 <- x.covar2[-idx.test2,]
# x.snp.all.train2 <- x.snp.all2[-idx.test2,]
# 
# 
# 
# 
# 
# 
# save(meta.result,file=paste0("/spin1/users/zhangh24/breast_cancer_data_analysis/risk_prediction/two_stage_model/result/meta.result",i1,".Rdata"))
# HeterVarianceEstimate <- function(log.odds,sigma){
#   M <- length(log.odds)
#   result <- (sum((log.odds-mean(log.odds))^2)-sum(diag(sigma))+sum(sigma)/M)/(M-1)
#   if(result <= 0){
#     result <- 0
#   }
#   return(result)
# }
# 
# 
# load("/spin1/users/zhangh24/breast_cancer_data_analysis/risk_prediction/standard_analysis/result/log.odds.meta.Rdata")
# M <- length(meta.result[[1]])
# beta0 <- log.odds.meta[i1]
# Sigma <- meta.result[[2]]
# betahat <- as.vector(meta.result[[1]])
# heter.sigma <- heter.variance.estimate(meta.result[[1]],meta.result[[2]])
# 
# 
# 
# save(log.odds.meta.la,file=paste0("/spin1/users/zhangh24/breast_cancer_data_analysis/risk_prediction/two_stage_model/result/log.odds.meta.la",i1,".Rdata"))
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# # 
# # Heter.result.Icog = EMmvpolySelfDesign(y.pheno.mis1.train,x.self.design = x.snp.all.train1,z.design=z.design,baselineonly = NULL,additive = x.covar.train1,pairwise.interaction = NULL,saturated = NULL,missingTumorIndicator = 888)
