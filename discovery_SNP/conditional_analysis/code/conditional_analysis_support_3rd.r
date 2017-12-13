

rm(list=ls())
args = commandArgs(trailingOnly = T)
i1 = as.numeric(args[[1]])

setwd("/spin1/users/zhangh24/breast_cancer_data_analysis/")
library(readr)
library(devtools)
library(CompQuadForm)
library(bc2)
library(data.table)
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
age1 <- data1[,204]
idx.complete1 <- which(age1!=888)
age1 <- age1[idx.complete1]

y.pheno.mis1 <- cbind(data1$Behaviour1,data1$ER_status1,data1$PR_status1,data1$HER2_status1,data1$Grade1)
colnames(y.pheno.mis1) = c("Behavior","ER","PR","HER2","Grade")

known.all.mis1 <- data1[,c(27:203)]
###fake the onco array only snp



x.covar.mis1 <- data1[,c(5:14)]


data2 <- fread("./data/Onco_euro_v10_10232017.csv",header=T)
data2 <- as.data.frame(data2)
age2 <- data2[,204]
idx.complete2 <- which(age2!=888)
age2 <- age2[idx.complete2]
y.pheno.mis2 <- cbind(data2$Behaviour1,data2$ER_status1,data2$PR_status1,data2$HER2_status1,data2$Grade1)
known.all.mis2 <- data2[,c(27:203,205)]
#y.pheno.mis2 <- cbind(data2$Behaviour1,data2$PR_status1,data2$ER_status1,data2$HER2_status1)
colnames(y.pheno.mis2) = c("Behaviour","ER",
                           "PR","HER2","Grade")

x.covar.mis2 <- data2[,c(5:14)]



icog.julie <- fread("/spin1/users/zhangh24/breast_cancer_data_analysis/data/Julie_snp_icog.csv")
icog.julie <- icog.julie[,-1]
discovery.snp.icog <- fread("/spin1/users/zhangh24/breast_cancer_data_analysis/data/discovery_icog_data.csv",header=T)
onco.julie <- fread("/spin1/users/zhangh24/breast_cancer_data_analysis/data/Julie_snp_onco.csv")
onco.julie <- onco.julie[,-1]
discovery.snp.onco <- fread("/spin1/users/zhangh24/breast_cancer_data_analysis/data/discovery_onco_data.csv")
x.discovery.mis1 <- as.data.frame(cbind(icog.julie,discovery.snp.icog))
x.discovery.mis2 <- as.data.frame(cbind(onco.julie,discovery.snp.onco))
###two snps onco array only
sudo.icog.na <- rep(NA,nrow(data1))

known.all.mis1 <- cbind(known.all.mis1,sudo.icog.na,x.discovery.mis1)
known.all.mis2 <- cbind(known.all.mis2,x.discovery.mis2)




y.pheno.mis1 <- y.pheno.mis1[idx.complete1,]
x.covar.mis1 <- x.covar.mis1[idx.complete1,]
x.covar.mis1 <- cbind(x.covar.mis1,age1)

known.all.mis1 <- known.all.mis1[idx.complete1,]
y.pheno.mis2 <- y.pheno.mis2[idx.complete2,]
x.covar.mis2 <- x.covar.mis2[idx.complete2,]
x.covar.mis2 <- cbind(x.covar.mis2,age2)

known.all.mis2 <- known.all.mis2[idx.complete2,]


load("/spin1/users/zhangh24/breast_cancer_data_analysis/discovery_SNP/conditional_analysis/result/icog.first.Rdata")
load("/spin1/users/zhangh24/breast_cancer_data_analysis/discovery_SNP/conditional_analysis/result/onco.first.Rdata")
load("/spin1/users/zhangh24/breast_cancer_data_analysis/discovery_SNP/conditional_analysis/result/conditional.results.first")

load("/spin1/users/zhangh24/breast_cancer_data_analysis/discovery_SNP/conditional_analysis/result/icog.2nd.Rdata")
load("/spin1/users/zhangh24/breast_cancer_data_analysis/discovery_SNP/conditional_analysis/result/onco.2nd.Rdata")
load("/spin1/users/zhangh24/breast_cancer_data_analysis/discovery_SNP/conditional_analysis/result/conditional.results.2nd")

icog.first.snpvalue <- icog.first[[2]][idx.complete1,]
onco.first.snpvalue <- onco.first[[2]][idx.complete2,]

icog.2nd.snpvalue <- icog.2nd[[2]][idx.complete1,]
onco.2nd.snpvalue <- onco.2nd[[2]][idx.complete2,]





n.first <- nrow(conditional.results.first)
first.known.flag <- conditional.results.first$known.flag
n.2nd <- nrow(conditional.results.2nd)
known.flag.2nd <- conditional.results.2nd$known.flag

if(i1 %in%known.flag.2nd){
  if(i1 ==178|i1==207|i1==172){
    
    idx.first <- which(i1==first.known.flag)
    idx.2nd <- which(i1==known.flag.2nd)
    snp.value.onco.2nd <- onco.2nd.snpvalue[,idx.2nd]
    
    first.snp.value.onco <- onco.first.snpvalue[,idx.first]
    snp.onco <- known.snp.value.onco <- known.all.mis2[,i1]
    
    x.all.mis2 <- cbind(snp.onco,known.snp.value.onco,
                        first.snp.value.onco,snp.value.onco.2nd,x.covar.mis2)
    
    
    
    
    score.test.support.icog <- NULL
    
    
    score.test.support.onco <- ScoreTestSupport(
      y.pheno.mis2,
      baselineonly = NULL,
      additive = x.all.mis2[,2:ncol(x.all.mis2)],
      pairwise.interaction = NULL,
      saturated = NULL,
      missingTumorIndicator = 888
    )
    save(score.test.support.icog,file=paste0("/spin1/users/zhangh24/breast_cancer_data_analysis/discovery_SNP/conditional_analysis/result/score.test.support.icog.2nd",i1,".Rdata"))
    save(score.test.support.onco,file=paste0("/spin1/users/zhangh24/breast_cancer_data_analysis/discovery_SNP/conditional_analysis/result/score.test.support.onco.2nd",i1,".Rdata"))
    
    
  }else{
    
    idx.first <- which(i1==first.known.flag)
    
    first.snp.value.icog <- icog.first.snpvalue[,idx.first]
    
    first.snp.value.onco <- onco.first.snpvalue[,idx.first]
    
    
    snp.icog <- known.snp.value.icog <- known.all.mis1[,i1]
    snp.onco <- known.snp.value.onco <- known.all.mis2[,i1]
    idx.2nd <- which(i1==known.flag.2nd)
    snp.value.icog.2nd <- icog.2nd.snpvalue[,idx.2nd]
    snp.value.onco.2nd <- onco.2nd.snpvalue[,idx.2nd]
    
    
    x.all.mis1 <- cbind(snp.icog,known.snp.value.icog,first.snp.value.icog,
                        snp.value.icog.2nd,
                        x.covar.mis1)
    x.all.mis2 <- cbind(snp.onco,known.snp.value.onco,first.snp.value.onco,
                        snp.value.onco.2nd,
                        x.covar.mis2)
    
    
    score.test.support.icog <- ScoreTestSupport(
      y.pheno.mis1,
      baselineonly = NULL,
      additive = x.all.mis1[,2:ncol(x.all.mis1)],
      pairwise.interaction = NULL,
      saturated = NULL,
      missingTumorIndicator = 888
    )
    score.test.support.onco <- ScoreTestSupport(
      y.pheno.mis2,
      baselineonly = NULL,
      additive = x.all.mis2[,2:ncol(x.all.mis2)],
      pairwise.interaction = NULL,
      saturated = NULL,
      missingTumorIndicator = 888
    )
    save(score.test.support.icog,file=paste0("/spin1/users/zhangh24/breast_cancer_data_analysis/discovery_SNP/conditional_analysis/result/score.test.support.icog.2nd",i1,".Rdata"))
    save(score.test.support.onco,file=paste0("/spin1/users/zhangh24/breast_cancer_data_analysis/discovery_SNP/conditional_analysis/result/score.test.support.onco.2nd",i1,".Rdata"))
    
    
  }
  
  
  
  
  
}else{
  
}



