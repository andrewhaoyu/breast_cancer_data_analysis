########Generate intrinsic subtypes function
SubtypesTrans <- function(casecon,ER,PR,HER2,grade){
  n <- length(casecon)
  result <- rep("unknown",n)
  idx.con <- which(casecon==0)
  result[idx.con] <- "control"
  idx.LA <- which((ER==1|PR==1)&HER2==0)
  idx.LB <- which((ER==1|PR==1)&HER2==1&grade!=3)
  idx.LB.HERneg <- which((ER==1|PR==1)&HER2==1&grade==3)
  idx.TN <- which(ER==0&PR==0&HER2==0)
  idx.HER2 <- which(ER==0&PR==0&HER2==1)
  result[idx.LA] <- "Luminal_A"
  result[idx.LB] <- "Luminal_B"
  result[idx.LB.HERneg] <- "Luminal_B_HER2neg"
  result[idx.HER2] <- "HER2_enriched"
  result[idx.TN] <- "triple_neg"
  return(result)
  
}
library(R.utils)
library(data.table)
setwd("/spin1/users/zhangh24/breast_cancer_data_analysis/")



data1 <- fread("./data/PRS_subtype_icgos_pheno_v10_euro.csv",header=T)
data1 <- as.data.frame(data1)
data2 <- fread("./data/PRS_subtype_Onco_euro_v10_08012018.csv",
               header=T)
data2 <- as.data.frame(data2)
library(xlsx)
#table(data1$study)
#write.xlsx(table(data1$study),file = "./data/study_sample_size_bcac.xlsx",sheetName="icog")
#write.xlsx(table(data2$study),file = "./data/study_sample_size_bcac.xlsx",sheetName="onco",append=T)
#cohort.study <- c("GS","KARMA","MMHS","Sisters",
 #                 "WGHS","PLCO")
cohort.study <- c("KARMA",
                  "MMHS",
                  "PLCO",
                  "SISTER",
                  "UKBGS")
############take out the cohort studies for validation
############leave other data for training and testing
idx1.co <- which(data1$study%in%cohort.study)
idx2.co <- which(data2$study%in%cohort.study)
data1.cohort <- data1[idx1.co,]
data2.cohort <-data2[idx2.co,]
data1.other <- data1[-idx1.co,]
data2.other <- data2[-idx2.co,]
###########take out the in-situ and unknown invasive people
idx1.co.clean <- which(data1.cohort$Behaviour1==0|
                         data1.cohort$Behaviour1==1)
idx2.co.clean <- which(data2.cohort$Behaviour1==0|
                         data2.cohort$Behaviour1==1)
data1.cohort.clean <- data1.cohort[idx1.co.clean,]
data2.cohort.clean <- data2.cohort[idx2.co.clean,]
subtypes.cohort1 <- SubtypesTrans(data1.cohort.clean$Behaviour1,
                                  data1.cohort.clean$ER_status1,
                                  data1.cohort.clean$PR_status1,
                                  data1.cohort.clean$HER2_status1,
                                  data1.cohort.clean$Grade1)

subtype2.cohort2 <- SubtypesTrans(data2.cohort.clean$Behaviour1,
                                  data2.cohort.clean$ER_status1,
                                  data2.cohort.clean$PR_status1,
                                  data2.cohort.clean$HER2_status1,
                                  data2.cohort.clean$Grade1)

t(table(subtype2.cohort2,data2.cohort.clean$study))
#write.csv(t(table(subtype2.cohort2,data2.cohort.clean$study)),file = "./risk_prediction/result/cohort_study_sample_size.csv")
idx.unknown <- which(subtype2.cohort2=="unknown")
head(data2.cohort.clean[idx.unknown,])

##########select testing dataset and training dataset

y.pheno.mis.other1 <- cbind(data1.other$Behaviour1,data1.other$ER_status1,data1.other$PR_status1,data1.other$HER2_status1,data1.other$Grade1)
y.pheno.mis.other2 <- cbind(data2.other$Behaviour1,data2.other$ER_status1,data2.other$PR_status1,data2.other$HER2_status1,data2.other$Grade1)
y.pheno.mis.other <- rbind(y.pheno.mis.other1,
                           y.pheno.mis.other2)
casecon.other <- y.pheno.mis.other[,1]
ER.other <- y.pheno.mis.other[,2]
PR.other <- y.pheno.mis.other[,3]
HER2.other <- y.pheno.mis.other[,4]
grade.other <- y.pheno.mis.other[,5]
subtypes.other <- SubtypesTrans(casecon.other,
                                ER.other,
                                PR.other,
                                HER2.other,
                                grade.other)
##########always rank the subtypes with the following orders:
##########Lu_A, Lu_B，Lu_B_HER2_neg, HER2_en_TN
subtypes.ratio <- table(subtypes.other)[c(3,4,5,2,6)]/sum( table(subtypes.other)[c(3,4,5,2,6)])
n.test.case <- round(table(casecon.other)[2]*0.1,0)
n.test.control <- round(table(casecon.other)[1]*0.1,0)
n.test.subtype <- round(n.test.case*subtypes.ratio,0)

##########test dataset: 1. Onco array 2. invasive 3. no oversampling for family history 4. individuals with unknown intrinsic subtypes 5. no studies of bilateral breast cancer
enriched.study <- c("2SISTER","ABCS-F","BBCS","BCFR-NY", "BCFR-PA", "BCFR-UT", "BOCS", "CNIO-BCS", "FHRISK","GC-HBOC", "HEBCS", "HEBON", "HKBCS", "IPOBCS", "KARBAC", "kConFab/AOCS", "KOHBRA", "MBCSG", "MSKCC","MYBRCA", "NBCS", "NC-BCFR", "OFBCR", "RBCS", "SUCCESSB", "SUCCESSC")
#########remove the data from onco array match the criteria
idx2.remove.test <- which((data2.other$study%in%enriched.study)|
                            data2.other$Behaviour1==2|
                            data2.other$Behaviour1==888)
data2.other.remove <- data2.other[-idx2.remove.test,]





y.pheno.mis2 <- cbind(data2.cohort$Behaviour1,data2.cohort$ER_status1,data2.cohort$PR_status1,data2.cohort$HER2_status1,data2.cohort$Grade1)
idx.mis <- which(y.pheno.mis2[,1]!=0&
                   y.pheno.mis2[,1]!=1)
y.pheno.mis2 <- y.pheno.mis2[-idx.mis,]
idx.mis <- which(y.pheno.mis2[,2]==888|y.pheno.mis2[,3]==888|y.pheno.mis2[,4]==888|y.pheno.mis2[,5]==888)
y.pheno.mis2 <- y.pheno.mis2[-idx.mis,]
length(idx.TN)+length(idx.HER2)+length(idx.LB.HERneg)+length(idx.LA)+length(idx.LB)
table(y.pheno.mis2[,1])
  
  


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
