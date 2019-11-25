#check the tumor characteristics correlations
setwd("/Users/zhangh24/GoogleDrive/breast_cancer_data_analysis")
#get the correlation between tumor characteristics for invasive breast cancer
library(data.table)
data1 <- as.data.frame(fread("./data/iCOGS_euro_v10_10232017.csv",header=T))
data1 <- as.data.frame(data1)
y.pheno.mis1 <- data.frame(Behaviour1=data1$Behaviour1,ER=data1$ER_status1,PR=data1$PR_status1,HER2=data1$HER2_status1,Grade=data1$Grade1,stringsAsFactors = F)
library(dplyr)
y.pheno1.complete <- y.pheno.mis1 %>% 
  filter(Behaviour1==1&
           (ER!=888&
              PR!=888&
              HER2!=888&
              Grade!=888))
data2 <- fread("./data/Onco_euro_v10_10232017.csv",header=T)
data2 <- as.data.frame(data2)
y.pheno.mis2 <- data.frame(Behaviour1=data2$Behaviour1,ER=data2$ER_status1,PR=data2$PR_status1,HER2=data2$HER2_status1,Grade=data2$Grade1,stringsAsFactors = F)
y.pheno2.complete <- y.pheno.mis2 %>% 
  filter(Behaviour1==1&
           (ER!=888&
              PR!=888&
              HER2!=888&
              Grade!=888))
y.pheno.complete <- rbind(y.pheno1.complete,
                          y.pheno2.complete)
cor(y.pheno.complete[,2:5])

#get the tumor correlation for in-situ and unknown invasiveness
data1_all <- as.data.frame(fread("./data/sig_snp_icog_prs.csv",header=T))
