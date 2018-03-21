rm(list=ls())
arg <- commandArgs(trailingOnly=T)
i1 <- as.numeric(arg[[1]])

library(bc2)
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
x.covar1 <- select(icog.data,5:14)
x.snp.all1 <- select(icog.data,26:230)
colnames(y.pheno.mis1)
y.pheno.mis2 <- select(onco.data,Behaviour1,ER_status1,PR_status1,HER2_status1,Grade1)
x.covar2 <- select(onco.data,5:14)
x.snp.all2 <- select(onco.data,26:230)
colnames(y.pheno.mis2)

idx.control1 <- which(icog.data$Behaviour1==0)
length(idx.control1)
idx.triple1 <- which(y.pheno.mis1[,2]==0&y.pheno.mis1[,3]==0&y.pheno.mis1[,4]==0)
length(idx.triple1)
idx.control2 <- which(onco.data$Behaviour1==0)
length(idx.control2)
idx.triple2 <- which(y.pheno.mis2[,2]==0&y.pheno.mis2[,3]==0&y.pheno.mis2[,4]==0)
length(idx.triple2)
#############random sample 200 cases & 200 controls from icog
#############random sample 500 cases & 500 controls from onco

n.test.control.icog <- 200
n.test.cases.icog <- 200
n.test.control.onco <- 500
n.test.cases.onco <- 500
set.seed(1)
idx.test.control.icog <- idx.control1[sample(length(idx.control1),n.test.control.icog)]
idx.test.triple.icog <- idx.triple1[sample(length(idx.triple1),n.test.cases.icog)]
idx.test1 <- c(idx.test.control.icog,idx.test.triple.icog)
y.pheno.mis1.test <- y.pheno.mis1[idx.test1,]
x.covar.test1 <- x.covar1[idx.test1,]
x.snp.all.test1 <- x.snp.all1[idx.test1,]


y.pheno.mis1.train <- y.pheno.mis1[-idx.test1,]
x.covar.train1 <- x.covar1[-idx.test1,]
x.snp.all.train1 <- x.snp.all1[-idx.test1,]




idx.test.control.onco <- idx.control2[sample(length(idx.control2),n.test.control.onco)]
idx.test.triple.onco <- idx.triple2[sample(length(idx.triple2),n.test.cases.onco)]
idx.test2 <- c(idx.test.control.onco,idx.test.triple.onco)
y.pheno.mis2.test <- y.pheno.mis2[idx.test2,]
x.covar.test2 <- x.covar2[idx.test2,]
x.snp.all.test2 <- x.snp.all2[idx.test2,]


y.test <- rbind(y.pheno.mis1.test,y.pheno.mis2.test)
x.snp.all.test <- rbind(as.matrix(x.snp.all.test1),as.matrix(x.snp.all.test2))

y.pheno.mis2.train <- y.pheno.mis2[-idx.test2,]
x.covar.train2 <- x.covar2[-idx.test2,]
x.snp.all.train2 <- x.snp.all2[-idx.test2,]


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
save(meta.result,file=paste0("/spin1/users/zhangh24/breast_cancer_data_analysis/risk_prediction/two_stage_model/result/meta.result",i1,".Rdata"))
heter.variance.estimate <- function(log.odds,sigma){
  M <- length(log.odds)
  result <- (sum((log.odds-mean(log.odds))^2)-sum(diag(sigma))+sum(sigma)/M)/(M-1)
  if(result <= 0){
    result <- 0
  }
  return(result)
}

load("/spin1/users/zhangh24/breast_cancer_data_analysis/risk_prediction/standard_analysis/result/log.odds.meta.Rdata")
M <- length(meta.result[[1]])
beta0 <- log.odds.meta[i1]
Sigma <- meta.result[[2]]
betahat <- as.vector(meta.result[[1]])
heter.sigma <- heter.variance.estimate(meta.result[[1]],meta.result[[2]])

if(heter.sigma==0){
  log.odds.meta.la <- rep(beta0,M)
}else{
  b <- sqrt(heter.sigma/2)
  data <- list(M=M,beta0=beta0,betahat=betahat,Sigma=Sigma,
               b=b)
  library(rstan)
  stan.model2 <- '
  data{
  //define data
  int<lower=1> M;
  real beta0;
  vector[M] betahat;
  matrix[M,M] Sigma;
  real<lower=0> b;
  }
  parameters{
  vector[M] beta;
  }
  model{
  //prior
  for(i in 1:M){
  beta[i] ~ double_exponential(beta0,b);
  }
  //data
  betahat ~ multi_normal(beta,Sigma);
  }
  '
  smodel <- stan_model(model_code = stan.model2)
  fit1 <- sampling(smodel,
                   data=data,
                   warmup=5000,
                   iter=10000,
                   control = list(adapt_delta=0.95),
                   chains=4)
  #traceplot(fit1,pars=c('beta'))
  
  log.odds.meta.la <- colMeans(extract(fit1)[[1]])
}


save(log.odds.meta.la,file=paste0("/spin1/users/zhangh24/breast_cancer_data_analysis/risk_prediction/two_stage_model/result/log.odds.meta.la",i1,".Rdata"))
















# 
# Heter.result.Icog = EMmvpolySelfDesign(y.pheno.mis1.train,x.self.design = x.snp.all.train1,z.design=z.design,baselineonly = NULL,additive = x.covar.train1,pairwise.interaction = NULL,saturated = NULL,missingTumorIndicator = 888)
