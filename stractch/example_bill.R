library(devtools)
install_github("andrewhaoyu/bc2", args = c('--library="/home/zhangh24/R/x86_64-pc-linux-gnu-library/3.4"'))
library(bc2)
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
i1 = 1
#if(i1<=177){
  ##analysis for Icog
  ##data1 contains the diseases, tumor characteristics, principle components and some known SNPs information
  data1 <- fread("./data/iCOGS_euro_v10_10232017.csv",header=T)
  data1 <- as.data.frame(data1)
  ##y pheno mis1 is the matrix of disease and tumor characteristics
  ##we need y pheno mis1 as the y side of regression
  y.pheno.mis1 <- cbind(data1$Behaviour1,data1$ER_status1,data1$PR_status1,data1$HER2_status1,data1$Grade1)
  colnames(y.pheno.mis1) = c("Behavior","ER","PR","HER2","Grade")
  # Grade1.fake <- data1$Grade1
  # Grade1.fake[data1$Grade1==2|data1$Grade1==3] <- 1
  # Grade1.fake[data1$Grade1==1] <- 0
  #y.pheno.mis1 <- cbind(data1$Behaviour1,data1$PR_status1,data1$ER_status1,data1$HER2_status1,Grade1.fake)
  # y.pheno.mis1 <- cbind(data1$Behaviour1,data1$PR_status1,data1$ER_status1,data1$HER2_status1)
  #x test all mis1 include the 177 known SNPs snpvalue
  x.test.all.mis1 <- data1[,c(27:203)]
  #x.test.all.mis1 <- 2-x.test.all.mis1
  ###pc1-10
  x.covar.mis1 <- data1[,c(5:14)]
  
  
  ####this x all mis1 include the all the covariate for a single snp
  x.all.mis1 <- as.matrix(cbind(x.test.all.mis1[,i1],x.covar.mis1))
  colnames(x.all.mis1)[1] <- "gene"
  
  ####EM mvpoly runs the two-stage model regression
  ####under current setting, we keep every covaraite second stage structure as additive model
  Heter.result.Icog = EMmvpoly(y.pheno.mis1,baselineonly = NULL,additive = x.all.mis1,pairwise.interaction = NULL,saturated = NULL,missingTumorIndicator = 888)
########EMmvpoly has multiple output
########EMmvpoly[[1]] are the estimate of second stage parameters for all of the covariates. When we use additve model, the number of second stage parameter is 5 for a single covariate, including baseline effecct, four main effect for tumor characteristics. For EMmvpoly[[1]], 1:M are the estimates for intercept, (M+1):(M+5) are the estimates for SNP, (M+6):(M+11) are estimates for priciple component1 and so on. 
########EMmvpoly[[2]] are the covariate matrix estimate for second stage parameters.
########EMmvpoly[[3]] are the data.frame format for second stage parameters
#######EMmvpoly[[4]] are the odds ratio, 95% CI and p value for the second stage parameters
#######EMmvpoly[[5]] are the global association test and global heterogneity test results
#######EMmvpoly[[6]] are the first stage paraters estimate under data frame format
#######EMmvpoly[[7]] are the odds ratio, 95% CI and p-value for the first stage parameters
#######EMmvpoly[[8]] are the log likelihood of the MLE
#######EMmvpoly[[9]] are the AIC of the model
#######EMmvpoly[[12]] are the standard z design matrix
  z.standard <- Heter.result.Icog[[12]]
  
  M <- nrow(z.standard)
  number.of.tumor <- ncol(z.standard)
  log.odds.icog <- Heter.result.Icog[[1]][(M+1):(M+1+number.of.tumor)]
  
  sigma.log.odds.icog <- Heter.result.Icog[[2]][(M+1):(M+1+number.of.tumor),(M+1):(M+1+number.of.tumor)]
 

  
  
########There are two versions of score tests
########first version is fixed effect score test
########second version is mixed effect score test
########To get the mixed effect score test, we need to get the fixed effect score results and random effect score result
########when we calculate the score test, we first need the score support
########this function is fitting the two-stage model under the null hypothesis 
########x.all.mis1[,2:ncol(x.all.mis1)] are the PC1-PC10
########It's slow and taking a lot of memory
########But under fixed effect score test, we only need to do it for one time, since the null hypothesis is the same for all SNPs
score.test.support.fixed.icog <- ScoreTestSupport(
    y.pheno.mis1,
    baselineonly = NULL,
    additive = x.all.mis1[,2:ncol(x.all.mis1)],
    pairwise.interaction = NULL,
    saturated = NULL,
    missingTumorIndicator = 888
  )
######Method_Nov21
######Under page 5, formula 4.1 is the formula for fixed effect score test statistics
######ScoreTestSupport[[1]] is called "inv_info_vec", it corresponds to (t(Z_C)*I_{etaeta}*Z_C)^(-1) in the formula, I saved it as a vector
######ScoreTestSupport[[2]] is called "YminusP", it corresponds to (Y-P)
######ScoreTestSupport[[3]] is called "Wobs", it corresponds to W
######ScoreTestSupport[[4]] is called "WXZ_vec", it corresponds to WCZ_C, and I saved it as a vector
######ScoreTestSupport[[5]] is called "zc", it corresponds to Z_C
######ScoreTestSuuport[[6]] to [[10]] are some supporting matrix I saved in R just to avoid recomputing




###now we have the score test support, we can run the score test
###we need three key elements for score test
###1.)the SNP, here the code is x=x.all.mis1[,1,drop=F]
###2.)the second stage design matrix for SNP, here we use additive
###second.stage.structure="additive"
###3.) score support we calculated last step
score.test.fixed.icog<- ScoreTest(y=y.pheno.mis1,
                              x=x.all.mis1[,1,drop=F],
                              second.stage.structure="additive",
                              score.test.support=score.test.support.fixed.icog,
                              missingTumorIndicator=888)
  ##the scoretest[[1]] are the score estimates for SNP
  ##the scoretest[[2]] are the information matrix for SNP
  score.fixed.icog <- score.test.fixed.icog[[1]]
  infor.fixed.icog <- score.test.fixed.icog[[2]]
  #DisplayFixedScoreTestResult will give you the corresponding score test p-value based on the score and information matrix
  DisplayFixedScoreTestResult(score.fixed.icog,infor.fixed.icog) 

  
  
############now let's go to the mixed effect score test, since we already got the fixed effect score test result, we only need to compute the random effect model
############The first thing is to generate score support
############The key difference here: the null hypothesis is different
############we allow some of the SNPs effect to be non zero here
############Under this example, we will choose baseline effect and ER effect as non zero effect
############PR, HER2, GRADE as random effect
############We need to create the corresponding matrix for ER
############Since we define some of the second stage parameters to be nonzero, it's not the standard way of two-stage model, we need to write the second stage design matrix by ourself
############z.design.ER is the second stage design matrix when we keep baseline effect and ER effect as non zero

  z.design.ER <- cbind(1,z.standard[,1])
########### x.self.design means the covairate you want to apply your self design second stage matrix, here x.all.mis1[,1] is the SNP
########### For all the other covariates (PC1-PC10), we used additive second stage design matrix
########### This ScoreTestSupportSelfDesign is esentially the same function ScoreTestSupport in previous section, it just add some freedom to second stage design matrix. But when we implement it in C, it's the same
    score.test.support.random.icog <- ScoreTestSupportSelfDesign(
    y.pheno.mis1,
    x.self.design  = x.all.mis1[,1,drop=F],
    z.design = z.design.ER,
    additive = x.all.mis1[,2:ncol(x.all.mis1)],
    pairwise.interaction = NULL,
    saturated = NULL,
    missingTumorIndicator = 888
  )
########For the Score test, we need the second stage design matrix for other tumor characteristics PR,HER2 and Grade
  z.design.other <-   z.standard[,2:4]
  
  ########### This ScoreTestSelfDesign is esentially the same function ScoreTest in previous section, it just add some freedom to second stage design matrix. But when we implement it in C, it's the same  
  
  score.test.random.icog<- ScoreTestSelfDesign(y=y.pheno.mis1,
                                                 x=x.all.mis1[,1,drop=F],
                                                 z.design=z.design.other,
                                                 score.test.support=score.test.support.random.icog,
                                                 missingTumorIndicator=888)
  score.random.icog <- score.test.random.icog[[1]]
  infor.random.icog <- score.test.random.icog[[2]]
######after we get the the fixed effect score, infor and random effect score, infor, we can combine them through the following function. Although mathmatically, they are very complex. But coding is very fast, this part doesn't need optimize
  test.result.second.mixed <- DisplayMixedScoreTestResult(
    score.fixed.icog,
    infor.fixed.icog,
    score.random.icog,
    infor.random.icog

                                                              )  
  
  #############This ScoreTestSupportSelfDesign was originally developed for fixed effect score test. It's slow and takes a lot of memory. But originally I thought it only need to be used once. But for random effect model, we need to use it for whole genome, so I did a little bit update on ScoreTestSupportSelfDesign and change it to ScoreTestSupportMixedModelSelfDesign. The corresponding ScoreTestSelfDesign function becomes ScoreTestMixedModel. ScoreTestSupportMixedModelSelfDesign will be faster and uses less memory compared to ScoreTestSupportSelfDesign. ScoreTestMixedModel will be slower compared to ScoreTestSelfDesign. In total, the new function will be faster than the old function. So I used the following function when I am doing whole genome analysis
#########I took out all the print function in   ScoreTestSupportMixedModelSelfDesign. Since running GWAS, printing would print too much in .o file

  
  score.test.support.random.icog <- ScoreTestSupportMixedModelSelfDesign(
    y.pheno.mis1,
    x.self.design  = x.all.mis1[,1,drop=F],
    z.design = z.design.ER,
    additive = x.all.mis1[,2:ncol(x.all.mis1)],
    pairwise.interaction = NULL,
    saturated = NULL,
    missingTumorIndicator = 888
  )
  ########For the Score test, we need the second stage design matrix for other tumor characteristics PR,HER2 and Grade
  z.design.other <-   z.standard[,2:4]
  
  ########### This ScoreTestSelfDesign is esentially the same function ScoreTest in previous section, it just add some freedom to second stage design matrix. But when we implement it in C, it's the same  
  
  score.test.random.icog<- ScoreTestMixedModel(y=y.pheno.mis1,
                                               x=x.all.mis1[,1,drop=F],
                                               z.design=z.design.other,
                                               score.test.support=score.test.support.random.icog,
                                               missingTumorIndicator=888)
  score.random.icog <- score.test.random.icog[[1]]
  infor.random.icog <- score.test.random.icog[[2]]
  ######after we get the the fixed effect score, infor and random effect score, infor, we can combine them through the following function. Although mathmatically, they are very complex. But coding is very fast, this part doesn't need optimize
  test.result.second.mixed <- DisplayMixedScoreTestResult(
    score.fixed.icog,
    infor.fixed.icog,
    score.random.icog,
    infor.random.icog
  )  
  