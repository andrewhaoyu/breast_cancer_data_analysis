#install_github("andrewhaoyu/bc2",ref='development', args = c('--library="/home/zhangh24/R/x86_64-pc-linux-gnu-library/3.4"'))
#install_github("andrewhaoyu/bc2", args = c('--library="/home/zhangh24/R/x86_64-pc-linux-gnu-library/3.4"'))
#install_github("andrewhaoyu/bc2", args = c('--library="/home/zhangh24/R/x86_64-pc-linux-gnu-library/3.4"'))
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
load("/spin1/users/zhangh24/breast_cancer_data_analysis/whole_genome/ICOG/ERPRHER2_fixed/result/extract_list.Rdata")
load("/spin1/users/zhangh24/breast_cancer_data_analysis/whole_genome/ONCO/ERPRHER2_fixed/result/extract_result.Rdata")

#idx <- which(extract.result[[1]]=="rs372562666:1:120561314:G:A")


dim(extract.result[[2]])


discovery.snp.onco <- as.data.frame(fread("/spin1/users/zhangh24/breast_cancer_data_analysis/discovery_SNP/result/discovery_onco_data.csv",header=T))
#x.test.all.mis1 <- discovery.snp.icog


discovery_snp <- read.csv("/spin1/users/zhangh24/breast_cancer_data_analysis/data/discovery_snp_summary_new.csv",header=T)

data2 <- fread("./data/Onco_euro_v10_10232017.csv",header=T)
data2 <- as.data.frame(data2)
LDandDp <- function(snp1,snp2){
  LD <- cor(snp1,snp2)^2
  p1 <- mean(snp1)
  p2 <- mean(snp2)
  #p12 <- cov(snp1,snp2)+p1*p2
  d <- cov(snp1,snp2)
  if(d<0){
    dmax <- min(p1*p2,(1-p1)*(1-p2))
  }else{
    dmax <- min(p1*(1-p2),(1-p1)*p2)
  }
  Dp <- d/dmax
  return(list(LD,Dp))
  
}





extract.result.onco <- extract.result[[2]]
idx.match <- match(extract.list$SNP.ONCO,extract.result[[1]])
extract.result.onco <- extract.result.onco[,idx.match]
idx.control <- which(data2$Behaviour1==0)
extract.result.onco.control <- extract.result.onco[idx.control,]
data_onco <- as.data.frame(fread("./data/ONCO_pruning.csv",header=T))
#x.test.all.mis2 <- data2[,c(27:205)]


# snp 5:45333860 and snp rs10941679
i.dis <- which(extract.list$CHR==12&
                 extract.list$position==115108136)
i.known <- which(colnames(data_onco)=="rs1292011")
idx.control <- which(data_onco$Behaviour1==0)
snp1 <- extract.result.onco.control[,i.dis]
mean(snp1)/2
snp2 <- as.vector(data_onco[idx.control,i.known])
LDandDp(snp1,snp2)
#colnames(data2)[i.known]
























#######snp 1:145126177 and snp rs12405132
i.dis <- which(colnames(discovery.snp.onco)=="1:145126177:G:A")
i.known <- which(colnames(data2)=="rs12405132")
idx.control <- which(data2$Behaviour1==0)
snp1 <- discovery.snp.onco[idx.control,i.dis]
mean(snp1)/2
snp2 <- as.vector(data2[idx.control,i.known])
LDandDp(snp1,snp2)


#######snp rs6677545 and snp rs35383942
i.dis <- which(colnames(discovery.snp.onco)=="rs6677545:200342046:A:C")
i.known <- which(colnames(data2)=="rs35383942")
idx.control <- which(data2$Behaviour1==0)
snp1 <- discovery.snp.onco[idx.control,i.dis]
snp2 <- as.vector(data2[idx.control,i.known])
LDandDp(snp1,snp2)

#######snp rs6677545 and snp rs6678914
i.dis <- which(colnames(discovery.snp.onco)=="rs6677545:200342046:A:C")
i.known <- which(colnames(data2)=="rs6678914")
idx.control <- which(data2$Behaviour1==0)
snp1 <- discovery.snp.onco[idx.control,i.dis]
snp2 <- as.vector(data2[idx.control,i.known])
LDandDp(snp1,snp2)



#######snp rs17743054:42900892:T:C and snp rs6507583
i.dis <- which(colnames(discovery.snp.onco)=="rs6677545:200342046:A:C")
i.known <- which(colnames(data2)=="rs6678914")
idx.control <- which(data2$Behaviour1==0)
snp1 <- discovery.snp.onco[idx.control,i.dis]
snp2 <- as.vector(data2[idx.control,i.known])
LDandDp(snp1,snp2)


#######snp rs16988381 and snp rs17879961
i.dis <- which(colnames(discovery.snp.onco)=="rs16988381")
i.known <- which(colnames(data2)=="rs17879961")
idx.control <- which(data2$Behaviour1==0)
snp1 <- discovery.snp.onco[idx.control,i.dis]
snp2 <- as.vector(data2[idx.control,i.known])
LDandDp(snp1,snp2)


#######snp rs16988381 and snp rs132390
i.dis <- which(colnames(discovery.snp.onco)=="rs16988381")
i.known <- which(colnames(data2)=="rs132390")
idx.control <- which(data2$Behaviour1==0)
snp1 <- discovery.snp.onco[idx.control,i.dis]
snp2 <- as.vector(data2[idx.control,i.known])
LDandDp(snp1,snp2)


#######snp rs16988381 and snp rs132390
i.dis <- which(colnames(discovery.snp.onco)=="rs16988381")
i.known <- which(colnames(data2)=="rs132390")
idx.control <- which(data2$Behaviour1==0)
snp1 <- discovery.snp.onco[idx.control,i.dis]
snp2 <- as.vector(data2[idx.control,i.known])
LDandDp(snp1,snp2)


#######snp 12:121435475:G:A and snp rs206966
i.dis <- which(colnames(discovery.snp.onco)=="12:121435475:G:A")
i.known <- which(colnames(data2)=="rs206966")
idx.control <- which(data2$Behaviour1==0)
snp1 <- discovery.snp.onco[idx.control,i.dis]
snp2 <- as.vector(data2[idx.control,i.known])
LDandDp(snp1,snp2)





i.dis <- which(colnames(discovery.snp.onco)=="1:145126177:G:A")
i.known <- which(colnames(data2)=="rs12405132")
idx.control <- which(data2$Behaviour1==0)
snp1 <- discovery.snp.onco[idx.control,i.dis]
snp2 <- as.vector(data2[idx.control,i.known])
LDandDp(snp1,snp2)

#######snp 12:29140260  and snp rs7297051
i.dis <- which(colnames(discovery.snp.onco)=="12:29140260:G:A")
i.known <- which(colnames(data2)=="rs7297051")
idx.control <- which(data2$Behaviour1==0)
snp1 <- discovery.snp.onco[idx.control,i.dis]
snp2 <- as.vector(data2[idx.control,i.known])
LDandDp(snp1,snp2)

