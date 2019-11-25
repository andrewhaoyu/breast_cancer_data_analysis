rm(list=ls())
arg <- commandArgs(trailingOnly=T)
i1 <- as.numeric(arg[[1]])
print(i1)

library(R.utils)

setwd("/data/zhangh24/breast_cancer_data_analysis/")

#subject.file <- "/gpfs/gsfs4/users/NC_BW/icogs_onco/genotype/imputed2/icogs_order.txt.gz"
subject.file <- "/gpfs/gsfs4/users/NC_BW/icogs_onco/genotype/imputed2/icogs_order.txt"
z.design <- matrix(c(
  c(0,1,1,1,0,0,0,0,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0),
  c(0,0,0,0,0,1,1,0,0,0,0,0,1,1,1,0,0,0,0,0,1,1,1),
  c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,0,0,0,0),
  c(0,0,0,0,1,0,0,0,0,0,0,1,0,0,0,0,0,0,0,1,0,0,0),
  c(1,0,0,0,0,0,0,1,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0)
),ncol=5)
rowSums(z.design)

colnames(z.design) <- c("Luminial A","Luminal B",
                        "Luminal B HER2Neg",
                        "HER2 Enriched",
                        "Triple Negative")
load(paste0("./risk_prediction/result/split.id.rdata"))
icog.train.id <- split.id[[1]]
onco.train.id <- split.id[[2]]
#onco.test.id <- split.id[[3]]
#icog.cohort.id <- split.id[[4]]
#onco.cohort.id <- split.id[[5]]
load("/data/zhangh24/breast_cancer_data_analysis/risk_prediction/Nasim_prs/result/icog.nasim.snp")
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

gc()

library(bc2, lib.loc ="/home/zhangh24/R/x86_64-pc-linux-gnu-library/3.6/")

load("./risk_prediction/intrinsic_subtypes_whole_genome/ICOG/result/delta0.icog.Rdata")

library(dplyr)
ID <- data1[,1,drop=F]
snp_test_all1 <- left_join(ID,icog.nasim.snp,by="ID")
snpvalue_all = snp_test_all1[,2:ncol(snp_test_all1)]
snpvalue = snpvalue_all[,i1]

Heter.result.Icog = EMmvpolySelfDesign(y.pheno.mis1,x.self.design = snpvalue,z.design=z.design,baselineonly = NULL,additive = x.covar.mis1,pairwise.interaction = NULL,saturated = NULL,missingTumorIndicator = 888,delta0 = delta0)
z.standard <- Heter.result.Icog[[12]]
M <- nrow(z.standard)
number.of.tumor <- ncol(z.standard)
log.odds.icog <- Heter.result.Icog[[1]][(M+1):(M+1+number.of.tumor)]
nparm <- length(Heter.result.Icog[[1]])  
sigma.log.odds.icog <- Heter.result.Icog[[2]][(M+1):(M+1+number.of.tumor),(M+1):(M+1+number.of.tumor)]



load("/data/zhangh24/breast_cancer_data_analysis/risk_prediction/Nasim_prs/result/onco.nasim.snp.rdata")
data2 <- as.data.frame(fread("/data/zhangh24/breast_cancer_data_analysis/data/sig_snp_onco_prs.csv",header=T))
data2 <- data2[,-1]
load("./risk_prediction/intrinsic_subtypes_whole_genome/ONCO/result/delta0.Rdata")

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
snp_test_all2 <- left_join(ID,onco.nasim.snp,by="ID")
snpvalue_all = snp_test_all2[,2:ncol(snp_test_all2)]
snpvalue = snpvalue_all[,i1]

Heter.result.Onco = EMmvpolySelfDesign(y.pheno.mis2,x.self.design = snpvalue,z.design = z.design,baselineonly = NULL,additive = x.covar.mis2,pairwise.interaction = NULL,saturated = NULL,missingTumorIndicator = 888, delta0=delta0)
nz.standard <- Heter.result.Onco[[12]]
M <- nrow(z.standard)
number.of.tumor <- ncol(z.standard)
log.odds.onco <- Heter.result.Onco[[1]][(M+1):(M+1+number.of.tumor)]
nparm <- length(Heter.result.Onco[[1]])
sigma.log.odds.onco <- Heter.result.Onco[[2]][(M+1):(M+1+number.of.tumor),(M+1):(M+1+number.of.tumor)]

y.pheno.mis2.temp <- y.pheno.mis2
x.covar.mis2.temp <- x.covar.mis2
z.design.temp = z.design
snpvalue.temp = snpvalue
log.odds.onco.temp = log.odds.onco

meta.result <- LogoddsMetaAnalysis(log.odds.icog,
                                   sigma.log.odds.icog,
                                   log.odds.onco,
                                   sigma.log.odds.onco)

second.stage.logodds.meta <- meta.result[[1]]
second.stage.sigma.meta <- meta.result[[2]]



test.result.second.wald <- DisplaySecondStageTestResult(second.stage.logodds.meta,second.stage.sigma.meta)

save( test.result.second.wald,file=paste0("/data/zhangh24/breast_cancer_data_analysis/risk_prediction/Nasim_prs/result/intrinsic_subtype_",i1,".Rdata"))

result <- list(second.stage.logodds.meta,second.stage.sigma.meta)
save(result, file = paste0("/data/zhangh24/breast_cancer_data_analysis/risk_prediction/Nasim_prs/result/intrinsic_subtype_logodds",i1,".Rdata"))



