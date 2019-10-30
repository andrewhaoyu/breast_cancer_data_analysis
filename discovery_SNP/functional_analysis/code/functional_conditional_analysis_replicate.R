#install_github("andrewhaoyu/bc2", args = c('--library="/home/zhangh24/R/x86_64-pc-linux-gnu-library/3.5"'))
###1 represent Icog
###2 represent Onco

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
library(bcutility)
library(data.table)
#############get the conditional analysis results for the SNPs with nearby known SNPs
load(paste0("./discovery_SNP/functional_analysis/result/functional_snp_conditional_replicate.Rdata"))
load("./discovery_SNP/functional_analysis/result/ICOG/functional_conditional_icog_snpvalue_replicate.Rdata")
load("./discovery_SNP/functional_analysis/result/ONCO/functional_conditional_onco_snpvalue_replicate.Rdata")
new_filter <- fread("/data/zhangh24/breast_cancer_data_analysis/data/two_new_SNPs_for_replicate.csv",header=T,stringsAsFactors = F)
idx.temp <- which(is.na(functional_snp_conditional$SNP.ICOG))
idx.temp2 <- which(is.na(functional_snp_conditional$SNP.ONCO))
#########only analysis the shared SNPs
if(length(idx.temp)!=0){
  functional_snp_conditional <- functional_snp_conditional[-idx.temp,]
  snpvalue.result.icog <- snpvalue.result.icog[,-idx.temp]
  snpvalue.result.onco <- snpvalue.result.onco[,-idx.temp]
}

num = nrow(functional_snp_conditional)
temp = diff(functional_snp_conditional$position)
size = 100
start.end <- startend(num,size,i1)
start = start.end[1]
end  = start.end[2]




p.value.result <- rep(0,(end-start+1))
temp = 1


data1 <- as.data.frame(fread("./data/icogs_overall.csv",header=T))
idx.case1 <- which(data1$Behaviour1==2|data1$Behaviour1==888)
data1$Behaviour1[idx.case1] <- 1
data1 <- as.data.frame(data1)
y.pheno.mis1 <- cbind(data1$Behaviour1,data1$ER_status1,data1$PR_status1,data1$HER2_status1,data1$Grade1)
colnames(y.pheno.mis1) = c("Behavior","ER","PR","HER2","Grade")

country1 <- as.factor(data1[,3])
data2 <- fread("./data/oncoarray_overall.csv",header=T)
data2 <- as.data.frame(data2)
idx.case2 <- which(data2$Behaviour1==2|data2$Behaviour1==888)
data2$Behaviour1[idx.case2] <- 1
y.pheno.mis2 <- cbind(data2$Behaviour1,data2$ER_status1,data2$PR_status1,data2$HER2_status1,data2$Grade1)
#y.pheno.mis2 <- cbind(data2$Behaviour1,data2$PR_status1,data2$ER_status1,data2$HER2_status1)
colnames(y.pheno.mis2) = c("Behaviour","ER",
                           "PR","HER2","Grade")

for(i in start:end){
  print(i)
  #####only two snps in the anlaysis
  #####chr 5 nearby known SNP is rs6596100
  #####chr 5 nearby known SNP is rs6507583
  chr = functional_snp_conditional$CHR[i]
  if(chr==5){
    snp.known.name <- "rs6596100"
  }else{
    snp.known.name <- "rs6507583"
  }
  idx.known1 <- which(colnames(data1)==snp.known.name)
  #known.idx <- DisKnown(dis.idx)
  SNP.ICOG <- functional_snp_conditional$SNP.ICOG[i]
  
  
  x.covar1 <- cbind(data1[,c(6:15)],data1[,idx.known1],country1)
  gene1 <- snpvalue.result.icog[,i]
  
  
  idx.known2 <- which(colnames(data2)%in%snp.known.name)
  country2 <- as.factor(data2[,4])
  #x.covar2 <- cbind(data2[,c(5:14)],country2)
  x.covar2 <- cbind(data2[,c(5:14)],data2[,idx.known2],country2)
  gene2 <- snpvalue.result.onco[,i]
  #age <- data2[,204]
  
  
  p.value.result[temp] = 
    two_data_standard_anlysis(y.pheno.mis1,
                              gene1,
                              x.covar1,
                              y.pheno.mis2,
                              gene2,
                              x.covar2)
  temp = temp+1
  
}

save(p.value.result,
     file = paste0("./discovery_SNP/functional_analysis/result/ICOG/p.value.result_replicate",i1,".Rdata"))


