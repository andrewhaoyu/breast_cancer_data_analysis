library(data.table)
library(bc2)
#data <- fread("./data/dataset_montse_20180522.txt")
setwd('/Users/haoyuzhang/GoogleDrive/breast_cancer_data_analysis (1)')
data <- fread("./data/Dataset_Montse_2018-10-10.txt")
##############we only focus on the invasive breast cancer cases
idx.invasive <- which(data$status==0|data$status==1)
data <- data[idx.invasive]


#############put the missing tumor characteristics as 888
idx.ER.mis <- which(data$status==1&is.na(data$ER_status1))
data$ER_status1[idx.ER.mis] <- 888
idx.PR.mis <- which(data$status==1&is.na(data$PR_status1))
data$PR_status1[idx.PR.mis] <- 888
idx.HER2.mis <- which(data$status==1&is.na(data$HER2_status1))
data$HER2_status1[idx.HER2.mis] <- 888
idx.grade.mis <- which(data$status==1&is.na(data$Grade1))
data$Grade1[idx.grade.mis] <- 888
#############put the subject BCAC-16687664 HER2 status as 1
idx <- which(data$HER2_status1==2)
data$HER2_status1[idx] <- 1
#############check the result
table(data$status,data$ER_status1)
table(data$status,data$PR_status1)
table(data$status,data$HER2_status1)
table(data$status,data$Grade1)
library(nnet)
############collapes breast mos cat 0,1
############collapes age fftp cat 0,1
data$breastmos_cat[data$breastmos_cat==0] <- 1
data$agefftp_cat[data$agefftp_cat==0] <- 1
###########create the dummy variable for breastmos_cat
###########create the dummy variable for parity_cat
###########create the dummy variable for agefftp_cat
breast_mat <- model.matrix(~as.factor(data$breastmos_cat)-1)[,-1]
colnames(breast_mat) <- paste0("breast_cat",c(2:5,9))
parity_mat <- model.matrix(~as.factor(data$parity_cat)-1)[,-1]
colnames(parity_mat) <- paste0("parity_cat",c(1:4,9))
agefftp_mat <- model.matrix(~as.factor(data$agefftp_cat)-1)[,-1]
colnames(agefftp_mat) <- paste0("agefftp",c(2:4,9))
refage <- data$refage
###########create the phenotype file
y <- cbind(data$status,data$ER_status1,data$PR_status1,
           data$HER2_status1,data$Grade1)
colnames(y) <- c("casecontrol",
                 "ER",
                 "PR",
                 "HER2",
                 "Grade")

############don't adjust for study
############population based study
idx1 <- which(data$design_cat==0)
###########two-stage model based on population based study
model.1 <- TwoStageModel(y=y[idx1,],
                         additive = cbind(agefftp_mat,
                                         breast_mat,
                                         parity_mat,
                                         refage)[idx1,],
                         missingTumorIndicator = 888
                         )
write.xlsx(model.1[[4]]
  ,file = "risk_factor_result_110118.xlsx",
           sheetName = "population_based_second_stage")
write.xlsx(model.1[[5]]
           ,file = "risk_factor_result_110118.xlsx",
           sheetName = "population_based_second_stage")






# model.1 <- multinom(molgroup~as.factor(agefftp_cat)+as.factor(breastmos_cat)+as.factor(parity_cat)+study+refage,data=data1, maxit= 500)
# coef(model.1)





idx2 <- which(data$design_cat==1)
data2 <- data[idx2,]
model.2 <- multinom(molgroup~as.factor(breastmos_cat)+study+refage,data=data2, maxit= 500)
coef(model.2)






model1 <- multinom(molgroup~as.factor(breastmos_cat)+as.factor(parity_cat)+as.factor(agefftp_cat)+study+refage,data=data1, maxit= 500)
coef(model1)
