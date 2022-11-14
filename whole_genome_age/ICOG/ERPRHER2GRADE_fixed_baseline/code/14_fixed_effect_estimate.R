#estimate the additive effect of ER, PR, HER2 and grade
#estimate the pairwise-interaction effect of ER, PR, HER2 and grade

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

#load icogs genotype data
#in your setting, you need to load the genotype data for 8 novel loci
discovery.snp.icog <- as.data.frame(fread("/data/zhangh24/breast_cancer_data_analysis/discovery_SNP/result/discovery_icog_data.csv",header=T))
colnames(discovery.snp.icog)
#load oncoarray genotype data
discovery.snp.onco <- as.data.frame(fread("/data/zhangh24/breast_cancer_data_analysis/discovery_SNP/result/discovery_onco_data.csv",header=T))
x.test.all.mis1 <- discovery.snp.icog
#load the snp information
discovery_snp <- read.csv("/data/zhangh24/breast_cancer_data_analysis/data/discovery_snp_summary_new.csv",header=T)

##analysis for Icog
#load phenotypes data for icogs
data1 <- fread("./data/iCOGS_euro_v10_10232017.csv",header=T)
data1 <- as.data.frame(data1)
y.pheno.mis1 <- cbind(data1$Behaviour1,data1$ER_status1,data1$PR_status1,data1$HER2_status1,data1$Grade1)
colnames(y.pheno.mis1) = c("Behavior","ER","PR","HER2","Grade")
###pc1-10 and age
x.covar.mis1 <- data1[,c(5:14,204)]
#always code the SNPs in terms of minor allele
if(discovery_snp$exp_freq_a1[i1]>0.5){
  x.test.all.mis1[,i1] = 2-x.test.all.mis1[,i1]
}
#prepare the covariate data for fitting two-stage model
#x.covar.mis1 <- as.matrix(cbind(x.test.all.mis1[,i1],x.covar.mis1))

age <- data1[,204]
idx.complete <- which(age!=888)
y.pheno.mis1 <- y.pheno.mis1[idx.complete,]
x.covar.mis1 <- x.covar.mis1[idx.complete,]
SNP = as.matrix(x.test.all.mis1[idx.complete,i1])
colnames(SNP) = "gene"


model.icog = TwoStageModel(y = y.pheno.mis1,
                           additive = x.covar.mis1,
                           pairwise.interaction = SNP,
                           missingTumorIndicator = 888)
model.icog[[4]]
#model.icog[[1]] contains the second.stage.parameter estimate. 
#The total length is 89. It contains intercepts for 23 subtypes, 55 (5*11) parameters for additive covariates
#11 parameters for pairwise.interactions
#model.icog[[2]] contains the corresponding covariance matrix for these 89 parameters
z.design.list = GenerateZDesignCombination(y.pheno.mis1)
z.additive = z.design.list[[1]]
z.interaction = z.design.list[[2]]
#M is the number of subtypes
M <- nrow(z.additive)
#dimension of second stage for additive
p.additive = ncol(z.additive)
#11 covariates with additive structure
n.additive = ncol(x.covar.mis1)
#dimension of second stage for pair-wise interaction
p.interaction = ncol(z.interaction)
n.interaction = 1
#find the log odds ratio for gene
start = M + p.additive*n.additive + 1
end = M + p.additive*n.additive + p.interaction*n.interaction
log.odds.icog <- model.icog[[1]][start:end]
#find covariance matrix for icogs
sigma.log.odds.icog <- model.icog[[2]][start:end,start:end]

#analysis for Onco Array
#data2 <- read.csv("./V10/Onco_euro_v10_05242017.csv",header=T)
data2 <- fread("./data/Onco_euro_v10_10232017.csv",header=T)
data2 <- as.data.frame(data2)
y.pheno.mis2 <- cbind(data2$Behaviour1,data2$ER_status1,data2$PR_status1,data2$HER2_status1,data2$Grade1)
#y.pheno.mis2 <- cbind(data2$Behaviour1,data2$PR_status1,data2$ER_status1,data2$HER2_status1)
colnames(y.pheno.mis2) = c("Behaviour","ER",
                           "PR","HER2","Grade")
#x.test.all.mis2 <- data2[,c(27:203)]
x.test.all.mis2 <- discovery.snp.onco
x.covar.mis2 <- data2[,c(5:14,204)]
if(discovery_snp$exp_freq_a1[i1]>0.5){
  x.test.all.mis2[,i1] = 2-x.test.all.mis2[,i1]
}
#sum(x.test.all.mis2[,i1])/(2*nrow(x.test.all.mis2))

ages <- data2[,204]
idx.complete <- which(ages!=888)
y.pheno.mis2 <- y.pheno.mis2[idx.complete,]
x.covar.mis2 <- x.covar.mis2[idx.complete,]
SNP = as.matrix(x.test.all.mis2[idx.complete,i1])
colnames(SNP)[1] = "gene"


model.onco = TwoStageModel(y = y.pheno.mis2,
                           additive = x.covar.mis2,
                           pairwise.interaction = SNP,
                           missingTumorIndicator = 888)
log.odds.onco <- model.onco[[1]][start:end]
#find covariance matrix for icogs
sigma.log.odds.onco <- model.onco[[2]][start:end,start:end]

meta.result <- LogoddsMetaAnalysis(log.odds.icog,
                                   sigma.log.odds.icog,
                                   log.odds.onco,
                                   sigma.log.odds.onco)

second.stage.logodds.meta <- meta.result[[1]]
second.stage.sigma.meta <- meta.result[[2]]



test.result.second.wald <- DisplaySecondStageTestResult(second.stage.logodds.meta,second.stage.sigma.meta)

model.icog[[4]]

heter.result <- list(test.result.second.wald)
save(heter.result,file=paste0("./discovery_SNP/additive_model/result/heter_result_",i1,".Rdata"))



# }else{
#   data2 <- fread("./data/Onco_euro_v10_10232017.csv",header=T)
#   data2 <- as.data.frame(data2)
#   names2 = colnames(data2)
#   y.pheno.mis2 <- cbind(data2$Behaviour1,data2$ER_status1,data2$PR_status1,data2$HER2_status1,data2$Grade1)
#   #y.pheno.mis2 <- cbind(data2$Behaviour1,data2$PR_status1,data2$ER_status1,data2$HER2_status1)
#   colnames(y.pheno.mis2) = c("Behaviour","ER",
#                              "PR","HER2","Grade")
#   idxi1 = 29
#   
#   #x.test.all.mis2 <- data2
#   x.covar.mis2 <- data2[,c(5:14)]
#   
#   
#   x.all.mis2 <- as.matrix(cbind(x.test.all.mis2[,idxi1],x.covar.mis2))
#   colnames(x.all.mis2)[1] = "gene"
#   
#   ages <- data2[,204]
#   idx.complete <- which(ages!=888)
#   
#   y.pheno.mis2 <- y.pheno.mis2[idx.complete,]
#   x.all.mis2 <- x.all.mis2[idx.complete,]
#   
#   
#   
#   Heter.result.Onco =EMmvpoly(y.pheno.mis2,baselineonly = NULL,additive = x.all.mis2,pairwise.interaction = NULL,saturated = NULL,missingTumorIndicator = 888)
#   z.standard <- Heter.result.Onco[[12]]
#   z.additive.design <- as.matrix(cbind(1,z.standard))
#   M <- nrow(z.standard)
#   number.of.tumor <- ncol(z.standard)
#   log.odds.onco <- Heter.result.Onco[[1]][(M+1):(M+1+number.of.tumor)]
#   nparm = length(Heter.result.Onco[[1]])
#   sigma.log.odds.onco <- Heter.result.Onco[[2]][(M+1):(M+1+number.of.tumor),(M+1):(M+1+number.of.tumor)]
#   beta.onco <- z.additive.design%*%log.odds.onco
#   beta.sigma.onco <- z.additive.design%*%sigma.log.odds.onco%*%t(z.additive.design)
#   loglikelihood.onco <- Heter.result.Onco[[8]]
#   
# 
#   
#   score.test.support.onco <- ScoreTestSupport(
#     y.pheno.mis2,
#     baselineonly = NULL,
#     additive = x.all.mis2[,2:ncol(x.all.mis2)],
#     pairwise.interaction = NULL,
#     saturated = NULL,
#     missingTumorIndicator = 888
#   )
#   score.test.onco<- ScoreTest(y=y.pheno.mis2,
#                               x=x.all.mis2[,1,drop=F],
#                               second.stage.structure="additive",
#                               score.test.support=score.test.support.onco,
#                               missingTumorIndicator=888)
#   z.design.score.baseline <- matrix(rep(1,M),ncol=1)
#   z.design.score.casecase <-z.standard
#   z.design.score.baseline.ER <- cbind(z.design.score.baseline,z.standard[,1])
#   z.design.score.casecase.ER <- z.standard[,2:ncol(z.standard)]
#   
#   score.onco <- score.test.onco[[1]]
#   infor.onco <- score.test.onco[[2]]
#   score.test.support.onco.baseline <- score.test.support.onco
#   #rm(score.test.support.onco)
#   rm(score.test.onco)
#   score.test.onco.baseline<- ScoreTestSelfDesign(y=y.pheno.mis2,
#                                                  x=x.all.mis2[,1,drop=F],
#                                                  z.design=z.design.score.baseline,
#                                                  score.test.support=score.test.support.onco.baseline,
#                                                  missingTumorIndicator=888)
#   
#   score.onco.baseline <- score.test.onco.baseline[[1]]
#   infor.onco.baseline <- score.test.onco.baseline[[2]]
#   rm(score.test.support.onco.baseline)
#   rm(score.test.onco.baseline)
#   
#   score.test.support.onco.casecase <- ScoreTestSupport(
#     y.pheno.mis2,
#     baselineonly = x.all.mis2[,1,drop=F],
#     additive = x.all.mis2[,2:ncol(x.all.mis2)],
#     pairwise.interaction = NULL,
#     saturated = NULL,
#     missingTumorIndicator = 888
#   )
#   score.test.onco.casecase<- ScoreTestSelfDesign(y=y.pheno.mis2,
#                                                  x=x.all.mis2[,1,drop=F],
#                                                  z.design=z.design.score.casecase,
#                                                  score.test.support=score.test.support.onco.casecase,
#                                                  missingTumorIndicator=888)
#   
#   
#   
#   score.onco.casecase <- score.test.onco.casecase[[1]]
#   infor.onco.casecase <- score.test.onco.casecase[[2]]
#   
#   z.design.score.baseline.heterER <- z.standard[,1,drop=F]
#   
#   
#   score.onco.baseline.heterER <- score.onco.casecase[1]
#   infor.onco.baseline.heterER <- infor.onco.casecase[1,1]
#   
#   
#   rm(score.test.support.onco.casecase)
#   rm(score.test.onco.casecase)
#   
#   
#   score.test.onco.baseline.ER<- ScoreTestSelfDesign(y=y.pheno.mis2,
#                                                     x=x.all.mis2[,1,drop=F],
#                                                     z.design=z.design.score.baseline.ER,
#                                                     score.test.support=score.test.support.onco,
#                                                     missingTumorIndicator=888)
#   
#   score.onco.baseline.ER <- score.test.onco.baseline.ER[[1]]
#   infor.onco.baseline.ER <- score.test.onco.baseline.ER[[2]]
#   
#   score.test.support.onco.casecase.ER <- ScoreTestSupportSelfDesign(
#     y.pheno.mis2,
#     x.self.design  = x.all.mis2[,1,drop=F],
#     z.design = z.design.score.baseline.ER,
#     additive = x.all.mis2[,2:ncol(x.all.mis2)],
#     pairwise.interaction = NULL,
#     saturated = NULL,
#     missingTumorIndicator = 888
#   )
#   
#   score.test.onco.casecase.ER<- ScoreTestSelfDesign(y=y.pheno.mis2,
#                                                     x=x.all.mis2[,1,drop=F],
#                                                     z.design=z.design.score.casecase.ER,
#                                                     score.test.support=score.test.support.onco.casecase.ER,
#                                                     missingTumorIndicator=888)
#   
#   score.onco.casecase.ER <- score.test.onco.casecase.ER[[1]]
#   infor.onco.casecase.ER <- score.test.onco.casecase.ER[[2]]
#   
#   rm(score.test.support.onco.casecase.ER)
#   rm(score.test.onco.casecase.ER)
#   gc()
#   
#   
#   
#   
#   
#   
#   second.stage.logodds.meta <-log.odds.onco
#   second.stage.sigma.meta <- sigma.log.odds.onco
#   
#   
#   
#   
#   
#   test.result.second.wald <- DisplaySecondStageTestResult(second.stage.logodds.meta,second.stage.sigma.meta)
#   
#   
#   beta.meta <- z.additive.design%*%second.stage.logodds.meta
#   beta.sigma.meta <- z.additive.design%*%second.stage.sigma.meta%*%t(z.additive.design)
#   
#   test.result.first.wald <- DisplayFirstStageTestResult(beta.meta,beta.sigma.meta)
#   
#   
#   score.meta <- score.onco
#   infor.meta <- infor.onco
#   
#   test.result.second.score <- DisplayFixedScoreTestResult(score.meta,infor.meta)
#   
#   
#   score.meta.baseline <- score.onco.baseline
#   infor.meta.baseline <- infor.onco.baseline
#   
#   score.meta.casecase <- score.onco.casecase
#   infor.meta.casecase <- infor.onco.casecase
#   
#   score.meta.baseline.heterER <- score.onco.baseline.heterER
#   infor.meta.baseline.heterER <- infor.onco.baseline.heterER
#   
#   test.result.second.mixed <- DisplayMixedScoreTestResult(score.meta.baseline,
#                                                           infor.meta.baseline,
#                                                           score.meta.casecase,
#                                                           infor.meta.casecase)  
#   test.result.second.mixed <- data.frame(t(test.result.second.mixed))
#   
#   
#   
#   
#   score.meta.baseline.ER <-     score.onco.baseline.ER
#   infor.meta.baseline.ER <- infor.onco.baseline.ER
#   
#   
#   score.meta.casecase.ER <-   score.onco.casecase.ER
#   infor.meta.casecase.ER <- infor.onco.casecase.ER
#   
#   
#   
#   
#   ###########DisplayMixedScoreTestResult is set up for one fixed effect,
#   ###########we need to adjust for the heterogeneity test
#   test.result.second.mixed.ER <- DisplayMixedScoreTestResult(score.meta.baseline.ER,
#                                                              infor.meta.baseline.ER,
#                                                              score.meta.casecase.ER,
#                                                              infor.meta.casecase.ER)  
#   test.result.second.mixed.ER <- data.frame(t(test.result.second.mixed.ER))
#   
#   test.result.second.mixed.ER.2 <- DisplayMixedScoreTestResult(score.meta.baseline.heterER,
#                                                                infor.meta.baseline.heterER,
#                                                                score.meta.casecase.ER,
#                                                                infor.meta.casecase.ER)  
#   test.result.second.mixed.ER.2 <- data.frame(t(test.result.second.mixed.ER.2))
#   
#   ###Both ER and RANDOM EFFECT should be counted into heterogeneity
#   test.result.second.mixed.ER[1,2] <-  test.result.second.mixed.ER.2[1,1]
#   
#   
#   colnames(test.result.second.mixed) <- c("mixed model global test for association(baseline fixed)","mixed model global test for heterogeneity(baseline fixed)")
#   colnames(test.result.second.mixed.ER) <- c("mixed model global test for association(baseline and ER fixed)","mixed model global test for heterogeneity(baseline and ER fixed)")
#   
#   
#   
#   
#   
#   
#   
#   loglikelihood <- loglikelihood.onco
#   AIC <- 2*length(Heter.result.Onco[[1]])-2*loglikelihood
#   
#   heter.result <- list(data.frame(test.result.second.wald,test.result.second.score, test.result.second.mixed,test.result.second.mixed.ER,loglikelihood = loglikelihood,AIC=AIC),
#                        data.frame(test.result.first.wald))
#   
#   save(heter.result,file=paste0("./discovery_SNP/additive_model/result/heter_result_",i1,".Rdata"))
#   
#   
#   
#   
# }