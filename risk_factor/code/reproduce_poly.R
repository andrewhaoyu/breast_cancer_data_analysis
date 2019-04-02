library(data.table)
library(bc2)
#data <- fread("./data/dataset_montse_20180522.txt")
setwd('/Users/zhangh24/GoogleDrive/breast_cancer_data_analysis')
data <- fread("./data/Dataset_Montse_20190322.txt")
#data <- fread("./data/Dataset_Montse_2018-10-10.txt")




##reproduce Jenny's group polytomous model
############collapes breast mos cat 0,1
############collapes age fftp cat 0,1
##############we only focus on the invasive breast cancer cases
data$staus[data$status==2|data$status==3] <- 1
data$breastmos_cat[data$breastmos_cat==0] <- 1
data$agefftp_cat[data$agefftp_cat==0] <- 1
data$lastchildage_cat[data$lastchildage_cat==0] <- 1


##############only focus on population based study


library(dplyr)
data1 = data %>% filter(design_cat==0)
###############polytomous model 
library(nnet)
model1 <- multinom(molgroup~as.factor(agemenarche_cat)+
                     as.factor(parity_cat)+
                     as.factor(mensagelast_cat)+
                     as.factor(agefftp_cat)+
                     as.factor(breastmos_cat)+
                     as.factor(lastchildage_cat)+
                     +study+refage,data=data1, maxit= 500)
coef.1 <- coef(model1)
covar.1 <- vcov(model1)
colnames(coef.1)[1:26] <- c("Intercept",
                            paste0("agemenarche_cat",c(1,2,3,9)),
                            paste0("parity_cat",c(1,2,3,4,9)),
                            paste0("mensagelast_cat",c(1,2,9)),
                            paste0("agefftp_cat",c(2,3,4,9)),
                            paste0("breastmos_cat",c(2,3,4,5,9)),
                            paste0("lastchildage_cat",c(2,3,4,9)))



#reorganize the covariance matrix to match coef data.frame
#columns are the variables
n.var <- ncol(coef.1)
#rows are the 5 intrinsic subtypes plus 9 
n.sub <- 6
var.1 <- coef.1
for(i in 1:n.var){
  var.1[,i] <- diag(covar.1[c(i+n.var*(0:(n.sub-1))),c(i+n.var*(0:(n.sub-1)))])
}

GenerateP <- function(coef,)













#delete <- na.action(model1)
#print(delete)


######puting all of the missing data as NA for mice package to run
agemenarche_cat <- data1$agemenarche_cat
idx <- which(agemenarche_cat==9)
agemenarche_cat[idx] <- NA
parity_cat <- data1$parity_cat
idx <- which(parity_cat==9)
parity_cat[idx] <- NA
mensagelast_cat <- data1$mensagelast_cat
idx <- which(mensagelast_cat==9)
mensagelast_cat[idx] <- NA
agefftp_cat <- data1$agefftp_cat
idx <- which(agefftp_cat==9)
agefftp_cat[idx] <- NA
breastmos_cat <- data1$breastmos_cat
idx <- which(breastmos_cat==9)
breastmos_cat[idx] <- NA
lastchildage_cat <- data1$lastchildage_cat
idx <- which(lastchildage_cat==9)
lastchildage_cat[idx] <- NA
study <- data1$study
refage <- data1$refage
all.covariates <- data.frame(as.factor(agemenarche_cat),
                               as.factor(parity_cat),
                               as.factor(mensagelast_cat),
                               as.factor(agefftp_cat),
                               as.factor(breastmos_cat),
                               as.factor(lastchildage_cat),
                               study,
                              refage)

# time1 = proc.time()
# imp <- mice(all.covariates,m=1,seed=1,print=FALSE)
# time = proc.time()-time1
# head(imp$imp$as.factor.agemenarche_cat.)
# all.covariates.c <- complete(imp,1)
# 
# x <- seq(0,1,0.01)
# y <- 1+sin(x)
# plot(x,y)
# data <- data.frame(x,y)
# library(ggplot2)
# ggplot(data)+geom_point(aes(x,y))
# M <- 10
# y_new <- 1+y/M
# ggplot(data)+geom_point(aes(x,y_new))
# 
# 













###########create the dummy variable for agemenarchecat
###########create the dummy variable for parity_cat
###########create the dummy variable for mensagelast_cat
###########create the dummy variable for breastmos_cat
###########create the dummy variable for agefftp_cat
###########create the dummy variable for lastchildage_mat
agemenarche_mat <- model.matrix(~as.factor(data1$agemenarche_cat)-1)[,-1]
colnames(agemenarche_mat) <- paste0("agemenarche_cat",c(1:3,9))
parity_mat <- model.matrix(~as.factor(data1$parity_cat)-1)[,-1]
colnames(parity_mat) <- paste0("parity_cat",c(1:4,9))
mensagelast_mat <- model.matrix(~as.factor(data1$mensagelast_cat)-1)[,-1]
colnames(mensagelast_mat) <- paste0("mensagelast_cat",c(1:2,9))
agefftp_mat <- model.matrix(~as.factor(data1$agefftp_cat)-1)[,-1]
colnames(agefftp_mat) <- paste0("agefftp",c(1:4,9))
breastmos_mat <- model.matrix(~as.factor(data1$breastmos_cat)-1)[,-1]
colnames(breastmos_mat) <- paste0("breastmos_cat",c(1:5,9))
lastchildage_mat <- model.matrix(~as.factor(data1$lastchildage_cat)-1)[,-1]
colnames(lastchildage_mat) <- paste0("lastchildage_cat",c(1:4,9))
refage <- data$refage





# idx.try <- which(data$design_cat==0&data$molgroup==2)
# data.new <- data[idx.try,]
# table(data.new$ER_status1,data.new$PR_status1,
#       data.new$HER2_status1,data.new$Grade1)

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


############ population based study




# model.1 <- multinom(molgroup~as.factor(agefftp_cat)+as.factor(breastmos_cat)+as.factor(parity_cat)+study+refage,data=data1, maxit= 500)
# coef(model.1)




############non population based study
idx2 <- which(data$design_cat==1)
data2 <- data[idx2,]
model.2 <- multinom(molgroup~as.factor(breastmos_cat)+study+refage,data=data2, maxit= 500)
coef(model.2)




