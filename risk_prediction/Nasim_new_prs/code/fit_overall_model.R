rm(list=ls())
arg <- commandArgs(trailingOnly=T)
i1 <- as.numeric(arg[[1]])
print(i1)

library(R.utils)

setwd("/data/zhangh24/breast_cancer_data_analysis/")

#subject.file <- "/gpfs/gsfs4/users/NC_BW/icogs_onco/genotype/imputed2/icogs_order.txt.gz"
subject.file <- "/gpfs/gsfs4/users/NC_BW/icogs_onco/genotype/imputed2/icogs_order.txt"
load(paste0("./risk_prediction/result/split.id.rdata"))
icog.train.id <- split.id[[1]]
onco.train.id <- split.id[[2]]
#onco.test.id <- split.id[[3]]
#icog.cohort.id <- split.id[[4]]
#onco.cohort.id <- split.id[[5]]
load("/data/zhangh24/breast_cancer_data_analysis/risk_prediction/Nasim_new_prs/result/icog.dis.snp")
library(data.table)
data1 <- as.data.frame(fread("/data/zhangh24/breast_cancer_data_analysis/data/sig_snp_icog_prs.csv",header=T))
data1 <- as.data.frame(data1[,-1])
icog.train <- which(data1[,1]%in%icog.train.id)
data1 <- data1[icog.train,]
y.pheno.mis1 <- cbind(data1$Behavior,data1$ER,data1$PR,data1$HER2,data1$Grade)
table(y.pheno.mis1[,1])
##########clean phenotype file
idx <- which(y.pheno.mis1[,1]==888|y.pheno.mis1[,1]==2)
y.pheno.mis1[idx,1] <- 1
colnames(y.pheno.mis1) = c("Behavior","ER","PR","HER2","Grade")

x.covar.mis1 <- data1[,c(7:16)]



library(bc2, lib.loc ="/home/zhangh24/R/x86_64-pc-linux-gnu-library/3.6/")

#load("./risk_prediction/intrinsic_subtypes_whole_genome/ICOG/result/delta0.icog.Rdata")

library(dplyr)
ID <- data1[,1,drop=F]
snp_test_all1 <- left_join(ID,icog.dis.snp,by="ID")
snpvalue_all = snp_test_all1[,2:ncol(snp_test_all1)]
snpvalue = snpvalue_all[,i1]

StandardLogisticmodel <- function(y.pheno.mis,snpvalue,x.covar.mis,idx){
  model <- glm(y.pheno.mis[idx,1]~snpvalue[idx]+as.matrix(x.covar.mis)[idx,],family=binomial())
  model.summary <- summary(model)
  log.odds <- coef(model.summary)[2,1]
  var.log.odds <- coef(model.summary)[2,2]^2
  return(list(log.odds,var.log.odds))
}
result.overall.icog <- StandardLogisticmodel(y.pheno.mis1,snpvalue,x.covar.mis1,idx=c(1:nrow(y.pheno.mis1)))




load("/data/zhangh24/breast_cancer_data_analysis/risk_prediction/Nasim_new_prs/result/onco.dis.snp.rdata")
data2 <- as.data.frame(fread("/data/zhangh24/breast_cancer_data_analysis/data/sig_snp_onco_prs.csv",header=T))
data2 <- data2[,-1]
#load("./risk_prediction/intrinsic_subtypes_whole_genome/ONCO/result/delta0.Rdata")

#onco.train <- which(data2[,1]%in%onco.train.id)
onco.train <- which(data2[,1]%in%onco.train.id)
data2 <- data2[onco.train,]
y.pheno.mis2 <- cbind(data2$Behavior,data2$ER,data2$PR,data2$HER2,data2$Grade)
#idx <- which(data2[,1]%in%split.id[[3]])
#######clean y.phneo.mis2 
idx <- which(y.pheno.mis2[,1]==888|y.pheno.mis2[,1]==2)
y.pheno.mis2[idx,1] <- 1
x.covar.mis2 <- data2[,c(7:16)]


library(dplyr)
ID <- data2[,1,drop=F]
snp_test_all2 <- left_join(ID,onco.dis.snp,by="ID")
snpvalue_all = snp_test_all2[,2:ncol(snp_test_all2)]
snpvalue = snpvalue_all[,i1]

result.overall.onco <- StandardLogisticmodel(y.pheno.mis2,snpvalue,x.covar.mis2,idx=c(1:nrow(y.pheno.mis2)))
meta.result.overall <- LogoddsMetaAnalysis(result.overall.icog[[1]],
                                           result.overall.icog[[2]],
                                           result.overall.onco[[1]],
                                           result.overall.onco[[2]])

logodds.overall <- meta.result.overall[[1]]
sigma.overall <- meta.result.overall[[2]]




result <- list(logodds.overall=logodds.overall,sigma.overall=sigma.overall)
save(result, file = paste0("/data/zhangh24/breast_cancer_data_analysis/risk_prediction/Nasim_new_prs/result/overall_log_odds_",i1,".Rdata"))
