library(R.utils)
library(data.table)
setwd("/spin1/users/zhangh24/breast_cancer_data_analysis/")



data1 <- fread("./data/PRS_subtype_icgos_pheno_v10_euro.csv",header=T)
data1 <- as.data.frame(data1)
data2 <- fread("./data/PRS_subtype_Onco_euro_v10_08012018.csv",
               header=T)

library(xlsx)
table(data1$study)
#write.xlsx(table(data1$study),file = "./data/study_sample_size_bcac.xlsx",sheetName="icog")
#write.xlsx(table(data2$study),file = "./data/study_sample_size_bcac.xlsx",sheetName="onco",append=T)
#cohort.study <- c("GS","KARMA","MMHS","Sisters",
 #                 "WGHS","PLCO")
cohort.study <- c("AHS","EPIC","FHRISK","KARMA",
                   "NHS","NHS2",
                  "PLCO","PROCAS",
                  "SISTER",
                  "UKBGS")
idx1 <- which(data1$study%in%cohort.study)
idx2 <- which(data2$study%in%cohort.study)
data2.cohort <-data2[idx2,]

y.pheno.mis2 <- cbind(data2.cohort$Behaviour1,data2.cohort$ER_status1,data2.cohort$PR_status1,data2.cohort$HER2_status1,data2.cohort$Grade1)
idx.mis <- which(y.pheno.mis2[,1]!=0&
                   y.pheno.mis2[,1]!=1)
y.pheno.mis2 <- y.pheno.mis2[-idx.mis,]
idx.mis <- which(y.pheno.mis2[,2]==888|y.pheno.mis2[,3]==888|y.pheno.mis2[,4]==888|y.pheno.mis2[,5]==888)
y.pheno.mis2 <- y.pheno.mis2[-idx.mis,]
idx.LA <- which((y.pheno.mis2[,2]==1|y.pheno.mis2[,3]==1)&y.pheno.mis2[,4]==0)
idx.LB <- which((y.pheno.mis2[,2]==1|y.pheno.mis2[,3]==1)&y.pheno.mis2[,4]==1&y.pheno.mis2[,5]!=3)
idx.LB.HERneg <- which((y.pheno.mis2[,2]==1|y.pheno.mis2[,3]==1)&y.pheno.mis2[,4]==1&y.pheno.mis2[,5]==3)
idx.TN <- which(y.pheno.mis2[,2]==0&y.pheno.mis2[,3]==0&y.pheno.mis2[,4]==0)
idx.HER2 <- which(y.pheno.mis2[,2]==0&y.pheno.mis2[,3]==0&y.pheno.mis2[,4]==1)
length(idx.TN)+length(idx.HER2)+length(idx.LB.HERneg)+length(idx.LA)+length(idx.LB)
table(y.pheno.mis2[,1])
  
  
enriched.study <- c("2SISTER","ABCS-F","BBCS","BCFR-NY", "BCFR-PA", "BCFR-UT", "BOCS", "CNIO-BCS", "FHRISK","GC-HBOC", "HEBCS", "HEBON", "HKBCS", "IPOBCS", "KARBAC", "kConFab/AOCS", "KOHBRA", "MBCSG", "MSKCC","MYBRCA", "NBCS", "NC-BCFR", "OFBCR", "RBCS", "SUCCESSB", "SUCCESSC")

y.pheno.mis1 <- cbind(data1$Behaviour1,data1$ER_status1,data1$PR_status1,data1$HER2_status1,data1$Grade1)
colnames(y.pheno.mis1) = c("Behavior","ER","PR","HER2","Grade")
#x.test.all.mis1 <- data1[,c(27:206)]
SG_ID <- data1$SG_ID
x.covar.mis1 <- data1[,c(5:14,204)]
age <- data1[,204]
idx.complete <- which(age!=888)
y.pheno.mis1 <- y.pheno.mis1[idx.complete,]
x.covar.mis1 <- x.covar.mis1[idx.complete,]
SG_ID <- SG_ID[idx.complete]




idx.fil <- Icog.order[,1]%in%SG_ID
idx.match <- match(SG_ID,Icog.order[idx.fil,1])
#Icog.order.match <- Icog.order[idx.fil,1][idx.match]
library(bc2)
