#install_github("andrewhaoyu/bc2",ref='development', args = c('--library="/home/zhangh24/R/x86_64-pc-linux-gnu-library/3.4"'))
###1 represent Icog
###2 represent Onco

rm(list=ls())
#commandarg <- commandArgs(trailingOnly=F)
#myarg <- commandarg[length(commandarg)]
#myarg <- sub("-","",myarg)
#i1 <- as.numeric(myarg)
args = commandArgs(trailingOnly = T)
i1 = as.numeric(args[[1]])


load(paste0("/spin1/users/zhangh24/breast_cancer_data_analysis/discovery_SNP/conditional_analysis/result/conditional.snp.list.icog.clean.sub",i1,".Rdata"))

load(paste0("/spin1/users/zhangh24/breast_cancer_data_analysis/discovery_SNP/conditional_analysis/result/conditional.snp.list.onco.clean.sub",i1,".Rdata"))
load("/spin1/users/zhangh24/breast_cancer_data_analysis/discovery_SNP/conditional_analysis/result/all.conditional.snps.Rdata")

setwd("/spin1/users/zhangh24/breast_cancer_data_analysis/")
library(readr)
library(devtools)
library(CompQuadForm)
library(bc2)
library(data.table)



n.condition <- nrow(all.conditional.snps)

start.end <- startend(n.condition,1000,i1)
start <- start.end[1]
end <- start.end[2]

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
snp.name.all.icog <- conditional.snp.list.icog.clean.sub[[1]]
x.test.all.mis1 <- conditional.snp.list.icog.clean.sub[[2]]
x.all.mis1 <- as.matrix(cbind(snp.icog,x.covar.mis1))

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

x.test.all.mis2 <- conditional.snp.list.onco.clean.sub[[2]]
snp.name.all.onco <- conditional.snp.list.onco.clean.sub[[1]]
x.covar.mis2 <- data2[,c(5:14)]

x.all.mis2 <- as.matrix(cbind(snp.onco,x.covar.mis2))
colnames(x.all.mis2)[1] = "gene"
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

known.flag.all <- all.conditional.snps$known.flag



for(i2 in 1:(end-start+1)){
  snp.name.icog <- snp.name.all.icog[i2]
  snp.icog <- x.test.all.mis1[,i2]
  snp.onco <- x.test.all.mis2[,i2]
  snp.name.onco <- snp.name.all.onco[i2]
  snp.icog <- snp.icog[idx.complete1]
  snp.onco <- snp.onco[idx.complete2]
  known.flag <- known.flag.all[start+i2-1]
    
  
  
}





