#-------------------------------------------------------------------
# Update Date: 11/26/2018
# Create Date: 11/26/2018
# Goal: split the BCAC data into traning, testing and cohort study design
# Author: Haoyu Zhang
#-------------------------------------------------------------------
#---------------------------------------#---------------------------------------
########Generate intrinsic subtypes function
SubtypesTrans <- function(casecon,ER,PR,HER2,grade){
  n <- length(casecon)
  result <- rep("unknown",n)
  idx.con <- which(casecon==0)
  result[idx.con] <- "control"
  idx.LA <- which(HER2==0&(ER==1|PR==1)&(grade==1|grade==2))
  idx.LB <- which(HER2==1&(ER==1|PR==1))
  idx.LUBHER2 <- which(HER2==0&(ER==1|PR==1)&grade==3)
  idx.HER2 <- which(HER2==1&ER==0&PR==0)
  idx.TN <- which(HER2==0&ER==0&PR==0)
  #idx.mis <- which(HER2==888|ER==888|PR==888|Grade==888)
  result[idx.LA] <- "Luminal_A"
  result[idx.LB] <- "Luminal_B"
  result[idx.LUBHER2] <- "Luminal_B_HER2Neg"
  result[idx.HER2] <- "HER2_Enriched"
  result[idx.TN] <- "TripleNeg"
  result <- factor(result,levels=c("control",
                                       "Luminal_A",
                                       "Luminal_B",
                                       "Luminal_B_HER2Neg",
                                       "HER2_Enriched",
                                       "TripleNeg",
                                       "unknown"))
  
  return(result)
  
}
library(R.utils)
library(data.table)
setwd("/data/zhangh24/breast_cancer_data_analysis/")



data1 <- fread("./data/PRS_subtype_icgos_pheno_v10_euro.csv",header=T)
data1 <- as.data.frame(data1)
data2 <- fread("./data/PRS_subtype_Onco_euro_v10_08012018.csv",
               header=T)
#temp =data2[,1]%in%onco.data[,1]
#idx <- which(data2[,1]%in%onco.data[,1]==F)
#data2[idx,]
data2 <- as.data.frame(data2)
library(xlsx)
################take out all the people with in-situ status
data1 = data1[data1$Behaviour1!=2,]
data2 = data2[data2$Behaviour1!=2,]

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
############
idx1.co <- which(data1$study%in%cohort.study)
idx2.co <- which(data2$study%in%cohort.study)
data1.cohort <- data1[idx1.co,]
data2.cohort <-data2[idx2.co,]
data1.other <- data1[-idx1.co,]
data2.other <- data2[-idx2.co,]
###########take out the unknown invasive people
idx1.co.clean <- which(data1.cohort$Behaviour1==0|
                         data1.cohort$Behaviour1==1)
idx2.co.clean <- which(data2.cohort$Behaviour1==0|
                         data2.cohort$Behaviour1==1)
data1.cohort.clean <- data1.cohort[idx1.co.clean,]
data2.cohort.clean <- data2.cohort[idx2.co.clean,]
id.cohort.clean1 <- data1.cohort.clean[,1]
id.cohort.clean2 <- data2.cohort.clean[,1]
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
write.csv(t(table(subtype2.cohort2,data2.cohort.clean$study)),file = "./risk_prediction/result/cohort_study_sample_size.csv")
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
#subtypes.ratio <- table(subtypes.other)[c(3,4,5,2,6)]/sum( table(subtypes.other)[c(3,4,5,2,6)])
subtypes.ratio <- table(subtypes.other)[c(2:6)]/sum(table(subtypes.other)[c(2:6)])

n.test.case <- round(table(casecon.other)[2]*0.1,0)
n.test.control <- round(table(casecon.other)[1]*0.1,0)
n.test.subtype <- round(n.test.case*subtypes.ratio,0)
n.test <- c(n.test.control,n.test.subtype)
names(n.test)[1] <- "control"
##########test dataset: 1. Onco array 2. invasive 3. no oversampling for family history 4. individuals with unknown intrinsic subtypes 5. no studies of bilateral breast cancer
enriched.study <- c("2SISTER","ABCS-F","BBCS","BCFR-NY", "BCFR-PA", "BCFR-UT", "BOCS", "CNIO-BCS", "FHRISK","GC-HBOC", "HEBCS", "HEBON", "HKBCS", "IPOBCS", "KARBAC", "kConFab/AOCS", "KOHBRA", "MBCSG", "MSKCC","MYBRCA", "NBCS", "NC-BCFR", "OFBCR", "RBCS", "SUCCESSB", "SUCCESSC")
#########remove the data from onco array match the criteria
idx2.remove.test <- which((data2.other$study%in%enriched.study)|
                            data2.other$Behaviour1==888)
data2.other.clean <- data2.other[-idx2.remove.test,]
id2.other.clean <- data2.other.clean[,1]
pheno.other.clean <-  cbind(data2.other.clean$Behaviour1,data2.other.clean$ER_status1,data2.other.clean$PR_status1,data2.other.clean$HER2_status1,data2.other.clean$Grade1)
casecon.other.clean <- pheno.other.clean[,1]
ER.other.clean <- pheno.other.clean[,2]
PR.other.clean <- pheno.other.clean[,3]
HER2.other.clean <- pheno.other.clean[,4]
grade.other.clean <- pheno.other.clean[,5]
subtypes.other.clean <- SubtypesTrans(casecon.other.clean,
                                      ER.other.clean,
                                      PR.other.clean,
                                      HER2.other.clean,
                                      grade.other.clean)
############set up training and testing
############the training and testing ratio is 9:1

set.seed(6)
subtypes.names <- names(n.test)

idx.test <- NULL
for(i in 1:length(subtypes.names)){
  idx <- which(subtypes.other.clean==subtypes.names[i])
  idx.test <- c(idx.test,sample(idx,n.test[i]))
}
id2.test <- id2.other.clean[idx.test]

#########all of the data in icog not in the cohort study are training
#########all of the data in onco array not in the cohort study, not in the testing data are training
id1.train <- data1.other[,1]
id2.train <- data2.other[,1][data2.other[,1]%in%id2.test!=T]

(length(id1.train)+length(id2.train))/length(id2.test)
data1.train <- data1[data1[,1]%in%id1.train,]
data2.train <- data2[data2[,1]%in%id2.train,]
data2.test <- data2[data2[,1]%in%id2.test,]
subtypes1.train <- SubtypesTrans(data1.train$Behaviour1,
                                 data1.train$ER_status1,
                                 data1.train$PR_status1,
                                 data1.train$HER2_status1,
                                 data1.train$Grade1)
subtype2.train <- SubtypesTrans(data2.train$Behaviour1,
                                data2.train$ER_status1,
                                data2.train$PR_status1,
                                data2.train$HER2_status1,
                                data2.train$Grade1)
subtypes2.test <- SubtypesTrans(data2.test$Behaviour1,
                                data2.test$ER_status1,
                                data2.test$PR_status1,
                                data2.test$HER2_status1,
                                data2.test$Grade1)
########number of studies and countries
########number of cases and controls
studies <- unique(c(data1.train$study,data2.train$study,data2.test$study))
countries <- unique(c(data1.train$StudyCountry,data2.train$StudyCountry,data2.test$StudyCountry))
casecon.icog <- c(data1.train$Behaviour1)
casecon.onco <- c(data2.train$Behaviour1,data2.test$Behaviour1)
table(casecon.icog)
table(casecon.onco)
sum(table(casecon.icog))
sum(table(casecon.onco))
sum(table(casecon.icog))+sum(table(casecon.onco))
table(casecon.icog)+table(casecon.onco)
casecon.onco.test <- data2.test$Behaviour1
casecon.onco.train <- data2.train$Behaviour1
sum(table(casecon.onco.test))
sum(table(casecon.onco.train)+table(casecon.icog))

############write out the training and testing sample size by study
write.xlsx(cbind(table(data1.train$study,subtypes1.train),table(data1.train$study,data1.train$Behaviour1)),
           file = "./risk_prediction/result/training_testing_data_sample_size_by_study.xlsx",
           sheetName="icog_training")
temp = cbind(table(data2.train$study,subtype2.train),table(data2.train$study,data2.train$Behaviour1))
colnames(temp)
sum(temp[,8:10])
sum(table(casecon.onco))
write.xlsx(cbind(table(data2.train$study,subtype2.train),table(data2.train$study,data2.train$Behaviour1)),
           file = "./risk_prediction/result/training_testing_data_sample_size_by_study.xlsx",
           sheetName = "onco_training",
           append=T)
write.xlsx(cbind(table(data2.test$study,subtypes2.test),table(data2.test$study,data2.test$Behaviour1)),
           file = "./risk_prediction/result/training_testing_data_sample_size_by_study.xlsx",
           sheetName="onco_testing",
           append=T
)
###############write out the id for training, testing and cohort study
split.id <- list(id1.train,
                 id2.train,
                 id2.test,
                 id.cohort.clean1,
                 id.cohort.clean2)
save(split.id,file = paste0("./risk_prediction/result/split.id.rdata"))



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
