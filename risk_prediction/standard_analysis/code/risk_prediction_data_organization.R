setwd("/spin1/users/zhangh24/breast_cancer_data_analysis/")

library(readr)
library(devtools)
library(CompQuadForm)
library(bc2)
library(data.table)
icog.julie <- fread("/spin1/users/zhangh24/breast_cancer_data_analysis/data/Julie_snp_icog.csv")
icog.julie <- icog.julie[,-1]
discovery.snp.icog <- fread("/spin1/users/zhangh24/breast_cancer_data_analysis/data/discovery_icog_data.csv",header=T)
onco.julie <- fread("/spin1/users/zhangh24/breast_cancer_data_analysis/data/Julie_snp_onco.csv")
onco.julie <- onco.julie[,-1]
discovery.snp.onco <- fread("/spin1/users/zhangh24/breast_cancer_data_analysis/data/discovery_onco_data.csv")
x.test.all.mis1 <- as.data.frame(cbind(icog.julie,discovery.snp.icog))
x.test.all.mis2 <- as.data.frame(cbind(onco.julie,discovery.snp.onco))

x.test.all.mis1 <- x.test.all.mis1[,-29]
x.test.all.mis2 <- x.test.all.mis2[,-29]


data1 <- fread("./data/iCOGS_euro_v10_10232017.csv",header=T)
x.covar.mis1 <- data1[,c(5:14)]




data2 <- fread("./data/Onco_euro_v10_10232017.csv",header=T) 
x.covar.mis2 <- data2[,c(5:14)]

newdata1 <- cbind(data1[,c(1:19,21:25,204,27:203)],x.test.all.mis1)
newdata2 <- cbind(data2[,c(1:19,21:25,204,27:203)],x.test.all.mis2)

write.csv(newdata1,file="/spin1/users/zhangh24/breast_cancer_data_analysis/data/sig_snps_icog.csv",row.names=F,quote = F)
write.csv(newdata2,file="/spin1/users/zhangh24/breast_cancer_data_analysis/data/sig_snps_onco.csv",
          row.names=F,quote = F)
