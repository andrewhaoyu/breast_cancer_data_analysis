# Goal: estimate log odds ratio and var of discovery SNPs only use complete data

rm(list=ls())
#commandarg <- commandArgs(trailingOnly=F)
#myarg <- commandarg[length(commandarg)]
#myarg <- sub("-","",myarg)
#i1 <- as.numeric(myarg)
args = commandArgs(trailingOnly = T)
i1 = as.numeric(args[[1]])

setwd("/data/zhangh24/breast_cancer_data_analysis/")
#setwd("/dcl01/chatterj/data/hzhang1/breast_cancer_data_analysis/")
library(bcutility,lib.loc = "/spin1/home/linux/zhangh24/R/x86_64-pc-linux-gnu-library/3.6/")
library(readr)
library(devtools)
library(bc2, lib.loc ="/home/zhangh24/R/x86_64-pc-linux-gnu-library/3.6/")
library(data.table)
discovery.snp.icog <- as.data.frame(fread("/data/zhangh24/breast_cancer_data_analysis/discovery_SNP/result/discovery_icog_data.csv",header=T))

#discovery.snp.icog <- as.data.frame(fread("./data/discovery_icog_data.csv",header=T))
colnames(discovery.snp.icog)
#onco.julie <- fread("/data/zhangh24/breast_cancer_data_analysis/data/Julie_snp_onco.csv")
#onco.julie <- onco.julie[,-1]
x.test.all.mis1 <- discovery.snp.icog

#sum(x.test.all.mis2[,11])/(2*nrow(x.test.all.mis2))

discovery_snp <- read.csv("./data/discovery_snp_summary_new.csv",header=T)

#x.test.all.mis1 <- as.data.frame(cbind(icog.julie,discovery.snp.icog))
#x.test.all.mis2 <- as.data.frame(cbind(onco.julie,discovery.snp.onco))


data1 <- fread("./data/iCOGS_euro_v10_10232017.csv",header=T)
data1 <- as.data.frame(data1)
y.pheno.mis1 <- cbind(data1$Behaviour1,data1$ER_status1,data1$PR_status1,data1$HER2_status1,data1$Grade1)
colnames(y.pheno.mis1) = c("Behavior","ER","PR","HER2","Grade")
x.covar.mis1 <- data1[,c(5:14,204)]
if(discovery_snp$exp_freq_a1[i1]>0.5){
  x.test.all.mis1[,i1] = 2-x.test.all.mis1[,i1]
}

age <- data1[,204]
snpvalue <- x.test.all.mis1[,i1]
#x.all.mis1 <- as.matrix(cbind(x.test.all.mis1[,i1],x.covar.mis1))
age <- data1[,204]
idx.complete <- which(age!=888)
snpvalue <- snpvalue[idx.complete]
y.pheno.mis1 <- y.pheno.mis1[idx.complete,]
#x.all.mis1 <- x.all.mis1[idx.complete,]
x.covar.mis1 <- x.covar.mis1[idx.complete,]
#colnames(x.all.mis1)[1] <- "gene"
subtypes.icog <- GenerateIntrinsicmis(y.pheno.mis1[,2],
                                      y.pheno.mis1[,3],
                                      y.pheno.mis1[,4],
                                      y.pheno.mis1[,5])
idx <- which(subtypes.icog!="mis")
y.pheno.mis1 = y.pheno.mis1[idx,]
x.covar.mis1 = x.covar.mis1[idx,]
snpvalue = snpvalue[idx]
subtypes.icog = as.character(subtypes.icog[idx])
subtypes.icog = as.factor(subtypes.icog)
levels(subtypes.icog) = c("control","Luminal_A",
                          "Luminal_B","Luminal_B_HER2Neg","HER2Enriched","TripleNeg")
library(nnet)
model1 <- multinom(subtypes.icog~ snpvalue+
                     as.matrix(x.covar.mis1), maxit= 500)
coef.1 <- coef(model1)
covar.1 <- vcov(model1)
jdx <- grep("snpvalue",colnames(covar.1))


log.odds.icog <- coef(model1)[,2]

sigma.log.odds.icog <- covar.1[jdx,jdx]


data2 <- fread("./data/Onco_euro_v10_10232017.csv",header=T)
data2 <- as.data.frame(data2)
y.pheno.mis2 <- cbind(data2$Behaviour1,data2$ER_status1,data2$PR_status1,data2$HER2_status1,data2$Grade1)
#y.pheno.mis2 <- cbind(data2$Behaviour1,data2$PR_status1,data2$ER_status1,data2$HER2_status1)
colnames(y.pheno.mis2) = c("Behaviour","ER",
                           "PR","HER2","Grade")
#x.test.all.mis2 <- data2[,c(27:203)]
# discovery.snp.onco <- as.data.frame(fread("./data/discovery_onco_data.csv",header=T))
discovery.snp.onco <- as.data.frame(fread("/data/zhangh24/breast_cancer_data_analysis/discovery_SNP/result/discovery_onco_data.csv",header=T))

x.test.all.mis2 <- discovery.snp.onco
x.covar.mis2 <- data2[,c(5:14,204)]
if(discovery_snp$exp_freq_a1[i1]>0.5){
  x.test.all.mis2[,i1] = 2-x.test.all.mis2[,i1]
}

# sum(x.test.all.mis1[,i1])/(2*nrow(x.test.all.mis1))
x.all.mis2 <- as.matrix(cbind(x.test.all.mis2[,i1],x.covar.mis2))
ages <- data2[,204]
idx.complete <- which(ages!=888)

y.pheno.mis2 <- y.pheno.mis2[idx.complete,]
x.all.mis2 <- x.all.mis2[idx.complete,]
x.covar.mis2 <- x.covar.mis2[idx.complete,]
colnames(x.all.mis2)[1] = "gene"

snpvalue = x.all.mis2[,1,drop=F]


subtypes.onco <- GenerateIntrinsicmis(y.pheno.mis2[,2],
                                      y.pheno.mis2[,3],
                                      y.pheno.mis2[,4],
                                      y.pheno.mis2[,5])
idx <- which(subtypes.onco!="mis")
y.pheno.mis2 = y.pheno.mis2[idx,]
x.covar.mis2 = x.covar.mis2[idx,]
snpvalue = snpvalue[idx]
subtypes.onco = as.character(subtypes.onco[idx])
subtypes.onco = as.factor(subtypes.onco)
levels(subtypes.onco) = c("control","Luminal_A",
                          "Luminal_B","Luminal_B_HER2Neg","HER2Enriched","TripleNeg")
library(nnet)
model2 <- multinom(subtypes.onco~ snpvalue+
                     as.matrix(x.covar.mis2), maxit= 500)
coef.2 <- coef(model2)
covar.2 <- vcov(model2)
jdx <- grep("snpvalue",colnames(covar.2))


log.odds.onco <- coef.2[,2]

sigma.log.odds.onco <- covar.2[jdx,jdx]


meta.result <- LogoddsMetaAnalysis(log.odds.icog,
                                   sigma.log.odds.icog,
                                   log.odds.onco,
                                   sigma.log.odds.onco)

second.stage.logodds.meta <- meta.result[[1]]
second.stage.sigma.meta <- meta.result[[2]]



test.result.second.wald <- DisplaySecondStageTestResult(second.stage.logodds.meta,second.stage.sigma.meta)

save( test.result.second.wald,file=paste0("./discovery_SNP/additive_model/result/intrinsic_subtype_logodds_complete_data_analysis_",i1,".Rdata"))
result <- list(second.stage.logodds.meta,second.stage.sigma.meta)
save(result, file = paste0("./discovery_SNP/additive_model/result/intrinsic_subtype_logodds_complete_data_analysis",i1,".Rdata"))



#}
