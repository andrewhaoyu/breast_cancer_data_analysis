#check the tumor characteristics correlations
setwd("/data/zhangh24/breast_cancer_data_analysis")
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
y.pheno1.all <- y.pheno.mis1 %>% 
  filter(Behaviour1==1)
data2 <- fread("./data/Onco_euro_v10_10232017.csv",header=T)
data2 <- as.data.frame(data2)
y.pheno.mis2 <- data.frame(Behaviour1=data2$Behaviour1,ER=data2$ER_status1,PR=data2$PR_status1,HER2=data2$HER2_status1,Grade=data2$Grade1,stringsAsFactors = F)
y.pheno2.complete <- y.pheno.mis2 %>% 
  filter(Behaviour1==1&
           (ER!=888&
              PR!=888&
              HER2!=888&
              Grade!=888))
y.pheno2.all <- y.pheno.mis2 %>% 
  filter(Behaviour1==1)
y.pheno.complete <- rbind(y.pheno1.complete,
                          y.pheno2.complete)
y.pheno.all <- rbind(y.pheno1.all,
                          y.pheno2.all)
invasvie.incomplete.rate <- 1-nrow(y.pheno.complete)/nrow(y.pheno.all)
cor_invasive <- cor(y.pheno.complete[,2:5])

#get the tumor correlation for in-situ and unknown invasiveness
data1 <- as.data.frame(fread("./data/sig_snp_icog_prs.csv",header=T))
data1 <- as.data.frame(data1)
y.pheno.mis1 <- data.frame(Behaviour1=data1$Behavior,ER=data1$ER,PR=data1$PR,HER2=data1$HER2,Grade=data1$Grade,stringsAsFactors = F)
library(dplyr)
y.pheno1.complete <- y.pheno.mis1 %>% 
  filter((Behaviour1==888|Behaviour1==2)&
           (ER!=888&
              PR!=888&
              HER2!=888&
              Grade!=888))
y.pheno1.all <- y.pheno.mis1 %>% 
  filter((Behaviour1==888|Behaviour1==2))
data2 <- as.data.frame(fread("./data/sig_snp_onco_prs.csv",header=T))
data2 <- as.data.frame(data2)
y.pheno.mis2 <- data.frame(Behaviour1=data2$Behavior,ER=data2$ER,PR=data2$PR,HER2=data2$HER2,Grade=data2$Grade,stringsAsFactors = F)
y.pheno2.complete <- y.pheno.mis2 %>% 
  filter((Behaviour1==888|Behaviour1==2)&
           (ER!=888&
              PR!=888&
              HER2!=888&
              Grade!=888))

y.pheno2.all <- y.pheno.mis2 %>% 
  filter((Behaviour1==888|Behaviour1==2))
y.pheno.complete <- rbind(y.pheno1.complete,
                          y.pheno2.complete)
cor_in_situ <- cor(y.pheno.complete[,2:5])
y.pheno.all <- rbind(y.pheno1.all,
                     y.pheno2.all)
insitu.incomplete.rate <- 1-nrow(y.pheno.complete)/nrow(y.pheno.all)


correlation_table <-rbind(cor_in_situ,cor_invasive)
write.csv(correlation_table,file = paste0("./discovery_SNP/tumor_characteristics_correlation/result/tumor_correlation_table.csv"))
