#install_github("andrewhaoyu/bcutility",args = c('--library="/home/zhangh24/R/x86_64-pc-linux-gnu-library/3.5"'))
#install_github("andrewhaoyu/bcutility")
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

genetic_covariance <- read.csv("/spin1/users/zhangh24/breast_cancer_data_analysis/risk_prediction/two_stage_model/result/genetic_covariance.csv",header=T)
genetic_correlation <- as.matrix(read.csv("/spin1/users/zhangh24/breast_cancer_data_analysis/risk_prediction/two_stage_model/result/genetic_correlation.csv",header=T))


# genetic_correlation <- genetic_correlation[c(2,4,5,3,1),c(2,4,5,3,1)]
#  write.csv(genetic_correlation,file = "/spin1/users/zhangh24/breast_cancer_data_analysis/risk_prediction/two_stage_model/result/genetic_correlation.csv",row.names = F)

#genetic_covariance <- genetic_covariance[,-c(1,2)]
# genetic_covariance <- genetic_covariance[,-6]
# genetic_covariance <- genetic_covariance[c(2,4,5,3,1),c(2,4,5,3,1)]
# write.csv(genetic_covariance,file = "/spin1/users/zhangh24/breast_cancer_data_analysis/risk_prediction/two_stage_model/result/genetic_covariance.csv",row.names = F)

insub.names <- colnames(z.design) <- c("Luminal_A","Luminal_B",
                        "Luminal_B_HER2Neg",
                        "HER2_Enriched",
                        "TripleNeg")

###############transform the genetic correlation



setwd('/spin1/users/zhangh24/breast_cancer_data_analysis/')

library(data.table)
icog.data <- as.data.frame(fread("/spin1/users/zhangh24/breast_cancer_data_analysis/data/sig_snp_icog_prs.csv",header=T))
icog.data <- icog.data[,-1]
onco.data <- as.data.frame(fread("/spin1/users/zhangh24/breast_cancer_data_analysis/data/sig_snp_onco_prs.csv",header=T))
onco.data <- onco.data[,-1]
library(tidyverse)
y.pheno.mis1 <- select(icog.data,Behavior,ER,PR,HER2,Grade)
############put the people with unknown case status as 1
idx.insi.unknown.icog <- which(y.pheno.mis1[,1]==888|
                            y.pheno.mis1[,1]==2) 
y.pheno.mis1[idx.insi.unknown.icog,1] = 1
subtypes.icog <- GenerateIntrinsicmis(y.pheno.mis1[,2],y.pheno.mis1[,3],
                                      y.pheno.mis1[,4],y.pheno.mis1[,5])

#table(subtypes.icog)+table(subtypes.onco)

x.covar1 <- select(icog.data,7:16)
x.snp.all1 <- select(icog.data,17:228)
colnames(y.pheno.mis1)
# split.id <- list(id1.train,
#                  id2.train,
#                  id2.test,
#                  id.cohort.clean1,
#                  id.cohort.clean2)
#load the training and testing data row number
load(paste0("./risk_prediction/result/split.id.rdata"))
#icog.test.id <- Generatetestid(subtypes.icog)
icog.train.id <- split.id[[1]]
onco.train.id <- split.id[[2]]
onco.test.id <- split.id[[3]]
icog.cohort.id <- split.id[[4]]
onco.cohort.id <- split.id[[5]]

y.pheno.mis2 <- select(onco.data,Behavior,ER,PR,HER2,Grade)
############put onco array unknown cases as 1
idx.insi.unknown.onco <- which(y.pheno.mis2[,1]==888|
                                 y.pheno.mis2[,1]==2) 
y.pheno.mis2[idx.insi.unknown.onco,1] <- 1
subtypes.onco <- GenerateIntrinsicmis(y.pheno.mis2[,2],
                                      y.pheno.mis2[,3],
                                      y.pheno.mis2[,4],
                                      y.pheno.mis2[,5])
x.covar2 <- select(onco.data,7:16)
x.snp.all2 <- select(onco.data,17:228)
colnames(y.pheno.mis2)
subtypes.onco <- GenerateIntrinsicmis(y.pheno.mis2[,2],y.pheno.mis2[,3],
                                      y.pheno.mis2[,4],y.pheno.mis2[,5])
#onco.test.id <- Generatetestid(subtypes.onco)



M <- 5
log.odds.standard.result <- 0
log.odds.poly.result <- rep(0,M)
log.odds.intrinsic.result <- rep(0,M)
sigma.intrinsic.result <- rep(0,M^2)
#log.odds.dic.result <-  matrix(0,6,M)
log.odds.intrinsic.eb <-  rep(0,M)
log.odds.intrinsic.ge <-  rep(0,M)
log.odds.intrinsic.ep <-  rep(0,M)
#log.odds.tree <- matrix(0,6,M)

#heter.variance.estiamte <- rep(0,6)


#for(i in 1:6){
  
  ##########split the data
  #idx.test.case <-  icog.test.id[[1]][[i]]
  #idx.test.control <- icog.test.id[[2]][[i]]
  #idx.test1 <- c(idx.test.case,idx.test.control)
  icog.train <- which(icog.data[,1]%in%icog.train.id)
  y.pheno.mis1.train <- y.pheno.mis1[icog.train,]
  x.covar.train1 <- x.covar1[icog.train,]
  x.snp.all.train1 <- x.snp.all1[icog.train,]
  #idx.test.case <-  onco.test.id[[1]][[i]]
  #idx.test.control <- onco.test.id[[2]][[i]]
  #idx.test2 <- c(idx.test.case,idx.test.control)
  onco.train <- which(onco.data[,1]%in%onco.train.id)
  y.pheno.mis2.train <- y.pheno.mis2[onco.train,]
  x.covar.train2 <- x.covar2[onco.train,]
  x.snp.all.train2 <- x.snp.all2[onco.train,]
  
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
  log.odds.standard.result <- rep(log.odds.meta,M)
  

  #############polytomous model
  insub.icog.train <- GenerateIntrinsicmis(y.pheno.mis1.train[,2],
                       y.pheno.mis1.train[,3],
                      y.pheno.mis1.train[,4],
                      y.pheno.mis1.train[,5])
  insub.onco.train <- GenerateIntrinsicmis(y.pheno.mis2.train[,2],
                                           y.pheno.mis2.train[,3],
                                           y.pheno.mis2.train[,4],
                                           y.pheno.mis2.train[,5])
  
  for(j in 1:M){
    ############generate the row id within icog intrinsic subtypes
    icog.train.sub.j <- which(insub.icog.train==insub.names[j]|
                              insub.icog.train=="control")
    n.snp <- ncol(x.snp.all.train1)
       model1 <- glm(y.pheno.mis1.train[icog.train.sub.j,1]~as.matrix(x.snp.all.train1[icog.train.sub.j,])+as.matrix(x.covar.train1[icog.train.sub.j,]),family = binomial(link = 'logit'))
      log.odds.icog <-as.vector(model1$coefficients)[c(1:n.snp)+1]
      sigma.icog <- as.matrix(vcov(model1)[c(1:n.snp)+1,c(1:n.snp)+1])
      onco.train.sub.j <- which(insub.onco.train==insub.names[j]|
                                  insub.onco.train=="control")
      model2 <- glm(y.pheno.mis2.train[onco.train.sub.j,1]~as.matrix(x.snp.all.train2[onco.train.sub.j,])+as.matrix(x.covar.train2[onco.train.sub.j,]),family = binomial(link = 'logit'))
      log.odds.onco <-as.vector(model2$coefficients)[c(1:n.snp)+1]
      sigma.onco <- vcov(model2)[c(1:n.snp)+1,c(1:n.snp)+1]
    
      meta.result <- LogoddsMetaAnalysis(log.odds.icog,
                                       sigma.icog,
                                       log.odds.onco,
                                       sigma.onco)
    log.odds.meta <- meta.result[[1]][i1]
    log.odds.poly.result[j] <- log.odds.meta
    
  }

  
  
  
  # 
  # idx.train.control.onco <- which(y.pheno.mis2.train[,1]==0)
  # freq <- sum(x.snp.all.train2[idx.train.control.onco,i1])/(2*length(idx.train.control.onco))
  #  
  
  #############intrinsic subtypes analysis
  Heter.result.Icog = EMmvpolySelfDesign(y.pheno.mis1.train,x.self.design = x.snp.all.train1[,i1],z.design=z.design,baselineonly = NULL,additive = x.covar.train1,pairwise.interaction = NULL,saturated = NULL,missingTumorIndicator = 888)
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
  # p.heter <- GlobalTestForHeter(meta.result[[1]],
  #                               meta.result[[2]],
  #                               self.design = T)
   logodds <- meta.result[[1]]
  sigma <- meta.result[[2]]
  library(mvtnorm)
  # classtree <- function(logodds,sigma){
  #   M <- length(logodds)
  #   pairs <- combn(M,2)
  #   pos <- ncol(pairs)
  #   logodds.pair <- matrix(0,pos,M)
  #   dis <- rep(0,pos)
  #   ##############check the distance between subtypes
  #   for(i in 1:pos){
  #     dis[i] <- PairDis(logodds[pairs[,i]],sigma[pairs[,i],pairs[,i]])
  #   }
  #   #############find the best distance
  #   idx.col <- which.min(dis)
  #   
  #   
  #       
  # }
  # 
  # PairDis <- function(logodds,sigma){
  #   trans <- c(1,-1)
  #   z = abs(trans%*%logodds)/sqrt(t(trans)%*%sigma%*%trans)
  #   return(z)
  # }
  # 
  
  
  
  
  
  
  log.odds.intrinsic.result <- meta.result[[1]] 
  sigma.intrinsic.result <- as.vector(meta.result[[2]])
  # log.odds.dic.result[i,] <- dicestiamte(meta.result[[1]],
  #                                    meta.result[[2]],
  #                                    log.odds.standard.result[i],
  #                                    p.heter,
  #                                    0.05)
  
  heter.variance.estiamte <- 
    HeterVarianceEstimate(meta.result[[1]],
                          meta.result[[2]])
   log.odds.intrinsic.eb <- ebestimate(
    meta.result[[1]],
    meta.result[[2]],
    log.odds.standard.result,
    heter.variance.estiamte 
  )
   heter.variance.estiamte.R <- HeterVarianceEstimateNew(
     meta.result[[1]],
     meta.result[[2]],
     genetic_correlation
   )
  prior.sigma <-    heter.variance.estiamte.R*as.matrix(genetic_correlation)
  log.odds.intrinsic.ge <- EbestimateNew(
   as.vector(meta.result[[1]]),
    meta.result[[2]],
    log.odds.standard.result,
   prior.sigma
  )
  # log.odds.tree[i,] <- tree_function(as.vector(meta.result[[1]]),
  #                                meta.result[[2]])
  # 
#}

all.model.result <- list(log.odds.standard.result
                         =log.odds.standard.result,
                         log.odds.intrinsic.result
                         =log.odds.intrinsic.result,
                         sigma.intrinsic.result
                         =sigma.intrinsic.result,
                         # log.odds.dic.result
                         # =log.odds.dic.result,
                         log.odds.intrinsic.eb
                         =log.odds.intrinsic.eb,
                         log.odds.intrinsic.ge
                         =log.odds.intrinsic.ge)
                         # heter.variance.estiamte
                         # =heter.variance.estiamte,
                         # log.odds.tree=log.odds.tree)

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
