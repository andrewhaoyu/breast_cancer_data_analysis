#goal: the mean and sd of age in controls and cases
setwd("/spin1/users/zhangh24/breast_cancer_data_analysis/")
data1 <- fread("./data/iCOGS_euro_v10_10232017.csv",header=T)
data1 <- as.data.frame(data1)
y.pheno.mis1 <- cbind(data1$Behaviour1,data1$ER_status1,data1$PR_status1,data1$HER2_status1,data1$Grade1)
colnames(y.pheno.mis1) = c("Behavior","ER","PR","HER2","Grade")
age <- data1[,204]
idx.complete <- which(age!=888)
y.pheno.mis1 <- y.pheno.mis1[idx.complete,]
##########
idx.case <- which(y.pheno.mis1[,1]==1)
age = age[idx.complete]
mean(age[idx.case])

data2 <- fread("./data/Onco_euro_v10_10232017.csv",header=T)
data2 <- as.data.frame(data2)
y.pheno.mis2 <- cbind(data2$Behaviour1,data2$ER_status1,data2$PR_status1,data2$HER2_status1,data2$Grade1)
#y.pheno.mis2 <- cbind(data2$Behaviour1,data2$PR_status1,data2$ER_status1,data2$HER2_status1)
colnames(y.pheno.mis2) = c("Behaviour","ER",
                           "PR","HER2","Grade")
#x.test.all.mis2 <- data2[,c(27:203)]
age2 <- data2[,204]
idx.complete2 <- which(age2!=888)
age2 <- age2[idx.complete2]
y.pheno.mis2 <- y.pheno.mis2[idx.complete2,]
age.all <- c(age,age2)
y.pheno.all <- rbind(y.pheno.mis1,y.pheno.mis2)
idx.case <- which(y.pheno.all[,1]==1)
mean(age.all[idx.case])
sd(age.all[idx.case])
idx.control <- which(y.pheno.all[,1]==0)
mean(age.all[idx.control])
sd(age.all[idx.control])
