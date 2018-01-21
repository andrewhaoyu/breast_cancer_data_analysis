rm(list=ls())
args = commandArgs(trailingOnly = T)
i1 = as.numeric(args[[1]])



setwd("/spin1/users/zhangh24/breast_cancer_data_analysis/")
library(readr)
library(devtools)
library(CompQuadForm)
library(bc2)
library(data.table)
library(bigmemory)
load("/spin1/users/zhangh24/breast_cancer_data_analysis/discovery_SNP/conditional_analysis/result/support.matrix.Rdata")
z.standard <- support.matrix[[1]]
z.additive.design <- support.matrix[[2]]
M <- support.matrix[[3]]
number.of.tumor <- support.matrix[[4]]
z.design.score.baseline <- support.matrix[[5]]
z.design.score.casecase <- support.matrix[[6]]
z.design.score.baseline.ER <- support.matrix[[7]]
z.design.score.casecase.ER <- support.matrix[[8]]

data1 <- fread("./data/iCOGS_euro_v10_10232017.csv",header=T)
data1 <- as.data.frame(data1)
age1 <- as.vector(data1[,204])
idx.complete1 <- which(age1!=888)
age1 <- age1[idx.complete1]


y.pheno.mis1 <- cbind(data1$Behaviour1,data1$ER_status1,data1$PR_status1,data1$HER2_status1,data1$Grade1)
colnames(y.pheno.mis1) = c("Behavior","ER","PR","HER2","Grade")

x.covar.mis1 <- data1[,c(5:14)]

y.pheno.mis1 <- y.pheno.mis1[idx.complete1,]
x.covar.mis1 <- x.covar.mis1[idx.complete1,]
x.covar.mis1 <- cbind(x.covar.mis1,age1)
data1.com <- data1[idx.complete1,]

data2 <- fread("./data/Onco_euro_v10_10232017.csv",header=T)
data2 <- as.data.frame(data2)
age2 <- data2[,204]
idx.complete2 <- which(age2!=888)
age2 <- age2[idx.complete2]
y.pheno.mis2 <- cbind(data2$Behaviour1,data2$ER_status1,data2$PR_status1,data2$HER2_status1,data2$Grade1)

#y.pheno.mis2 <- cbind(data2$Behaviour1,data2$PR_status1,data2$ER_status1,data2$HER2_status1)
colnames(y.pheno.mis2) = c("Behaviour","ER",
                           "PR","HER2","Grade")

y.pheno.mis2 <- y.pheno.mis2[idx.complete2,]

x.covar.mis2 <- data2[,c(5:14)]
x.covar.mis2 <- x.covar.mis2[idx.complete2,]
x.covar.mis2 <- cbind(x.covar.mis2,age2)
data2.com <- data2[idx.complete2,]
discovery.snp.icog <- as.data.frame(fread("/spin1/users/zhangh24/breast_cancer_data_analysis/data/discovery_icog_data.csv",header=T))

discovery.snp.icog.complete <- discovery.snp.icog[idx.complete1,]

discovery.snp.onco <- as.data.frame(fread("/spin1/users/zhangh24/breast_cancer_data_analysis/data/discovery_onco_data.csv"))
discovery.snp.onco.complete <- discovery.snp.onco[idx.complete2,]



n.coun <- length(all.counries)
all.countries <- unique(c(data1$StudyCountry,data2$StudyCountry))
for(i2 in 1:n.coun){
i2 = i2+1
  idx.icog <- which(data1.com$StudyCountry==all.countries[i2])
  print(length(idx.icog))
  y.pheno.mis1.sub <- y.pheno.mis1[idx.icog,]
  x.covar.mis1.sub <- x.covar.mis1[idx.icog,]
  discovery.snp.icog.sub <- discovery.snp.icog.complete[idx.icog,i1]
  x.all.mis1.sub <- cbind(discovery.snp.icog.sub,x.covar.mis1.sub)
  Heter.result.Icog = TwoStageModel(y.pheno.complete,baselineonly = NULL,additive = x.all.mis1.sub.complete,pairwise.interaction = NULL,saturated = NULL,missingTumorIndicator = NULL)
  
  delta0 <- Heter.result.Icog[[1]]
  Heter.result.Icog = EMmvpoly(y.pheno.mis1.sub,baselineonly = NULL,additive = x.all.mis1.sub,pairwise.interaction = NULL,saturated = NULL,missingTumorIndicator = 888,delta0=delta0)
  
  idx.onco <- which(data2.com$StudyCountry==all.countries[i2])
  print(length(idx.onco))
  y.pheno.mis2.sub <- y.pheno.mis2[idx.onco,]
  x.covar.mis2.sub <- x.covar.mis2[idx.onco,]
  discovery.snp.onco.sub <- discovery.snp.onco.complete[idx.onco,i1]
  x.all.mis2.sub <- cbind(discovery.snp.onco.sub,x.covar.mis2.sub)
  
  Heter.result.Onco = EMmvpoly(y.pheno.mis2.sub,baselineonly = NULL,additive = x.all.mis2.sub,pairwise.interaction = NULL,saturated = NULL,missingTumorIndicator = 888)
  
  
}




delta0 <-   Heter.result.Onco[[1]]




idx.onco <- which(data2$StudyCountry=="Greece")

Heter.result.Onco = EMmvpoly(y.pheno.mis2.sub,baselineonly = NULL,additive = x.all.mis2.sub,pairwise.interaction = NULL,saturated = NULL,missingTumorIndicator = 888)

z.standard <- Heter.result.Onco[[12]]
M <- nrow(z.standard)
number.of.tumor <- ncol(z.standard)
log.odds.onco <- Heter.result.Onco[[1]][(M+1):(M+1+number.of.tumor)]
nparm <- length(Heter.result.Onco)
sigma.log.odds.onco <- Heter.result.Onco[[2]][(M+1):(M+1+number.of.tumor),(M+1):(M+1+number.of.tumor)]
beta.onco <- z.additive.design%*%log.odds.onco
beta.sigma.onco <- z.additive.design%*%sigma.log.odds.onco%*%t(z.additive.design)
loglikelihood.onco <- Heter.result.Onco[[8]]
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







for(i2 in 1:(end-start+1)){
  print(i2)
  snp.name.icog <- snp.name.all.icog[i2]
  
  snp.icog <- x.test.all.mis1[,i2]
  snp.onco <- x.test.all.mis2[,i2]
  snp.name.onco <- snp.name.all.onco[i2]
  snp.icog <- snp.icog[idx.complete1]
  snp.onco <- snp.onco[idx.complete2]
  known.flag <- known.flag.all[start+i2-1]
  known.flag.new<- known.flag
  if(known.flag.new!=known.flag.last){
    load(paste0("/spin1/users/zhangh24/breast_cancer_data_analysis/discovery_SNP/conditional_analysis/result/score.test.support.icog",known.flag,".Rdata"))
    load(paste0("/spin1/users/zhangh24/breast_cancer_data_analysis/discovery_SNP/conditional_analysis/result/score.test.support.onco",known.flag,".Rdata"))
    known.flag.last <- known.flag.new
  }
  
  p.value.all[i2] <- condition_additive_model(y.pheno.mis1,
                                              x.covar.mis1,
                                              snp.name.icog,
                                              snp.icog,
                                              y.pheno.mis2,
                                              x.covar.mis2,
                                              snp.name.onco,
                                              snp.onco,
                                              known.flag,
                                              known.all.mis1,
                                              known.all.mis2,
                                              z.standard,
                                              z.additive.design,
                                              M,
                                              number.of.tumor,
                                              z.design.score.baseline,
                                              z.design.score.casecase,
                                              z.design.score.baseline.ER,
                                              z.design.score.casecase.ER,
                                              score.test.support.icog = score.test.support.icog,
                                              score.test.support.onco = score.test.support.onco,
                                              region.all = region.all
  )
  
  
}


save(p.value.all,file=paste0("/spin1/users/zhangh24/breast_cancer_data_analysis/discovery_SNP/conditional_analysis/result/psub",i1,".Rdata"))

