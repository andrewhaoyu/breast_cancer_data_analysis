args = commandArgs(trailingOnly = T)
i1 = as.numeric(args[[1]])


library(data.table)
library(bc2)
#data <- fread("./data/dataset_montse_20180522.txt")
setwd('/dcl01/chatterj/data/hzhang1/breast_cancer_data_analysis')
data <- as.data.frame(fread("./data/Dataset_Montse_2018-10-10.txt"))
##############we only focus on the invasive breast cancer cases
idx.invasive <- which(data$status==0|data$status==1)
data <- data[idx.invasive,]
##############take out all the people with missing parous
idx <- which(is.na(data$parous)!=T)
data <- data[idx,]
#############take out all the people with missing parity, breastmos, or agfftp
idx <- which((is.na(data$parity)!=T)&
               is.na(data$breastMos)!=T&
               is.na(data$ageFFTP)!=T)
data <- data[idx,]

# try <- which(is.na(data$parity)==T)
# try2 <- which(is.na(data$breastMos)==T)
# try3 <- which(is.na(data$ageFFTP)==T)

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
data$HER2_status1[idx] <- 888
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
agefftp <- data$ageFFTP
breastmos <- data$breastMos/6
parity <- data$parity
###########create the phenotype file
y <- cbind(data$status,data$ER_status1,data$PR_status1,
           data$HER2_status1,data$Grade1)
colnames(y) <- c("casecontrol",
                 "ER",
                 "PR",
                 "HER2",
                 "Grade")
idx <- which(data$status==0)
y[idx,2:ncol(y)] <- NA

###########create study mat
idx <- which(data$design_cat==0)
study1 <- data[idx,2]
idx <- which(data$design_cat==1)
study2 <- data[idx,2]
study_mat1 <- model.matrix(~as.factor(study1)-1)[,-1]
study_mat2 <- model.matrix(~as.factor(study2)-1)[,-1]
############don't adjust for study
############population based study
idx <- which(data$design_cat==0)
study_names <- names(table(data[idx,2]))
###########create the phenotype file
z.standard <- GenerateZstandard(y)
z.ER <- cbind(1,z.standard[,1])

if(i1==7){
  idx <- which(data$design_cat==0)
  z.additive = cbind(1,z.standard)
  model <- EMmvpolySelfDesignnew(y=y[idx,],
                         z.design = z.additive,
                         x.self.design = cbind(agefftp,
                                               breastmos,
                                               parity,
                                               refage)[idx,],
                         z.design2 = z.ER,
                         x.self.design2 = study_mat1,
                         missingTumorIndicator = 888
  )
  save(model,file=paste0("./risk_factor/result/model",i1,".Rdata"))
  
}else if(i1 ==8){
  z.standard <- GenerateZstandard(y)
  ############find triple negative row
  idx <- which(z.standard[,1]==0&
                 z.standard[,2]==0&
                 z.standard[,3]==0)
  z.triple <- rep(1,nrow(z.standard))
  z.triple[idx] <- 0
  z.design <- cbind(1,z.standard,z.triple)
  colnames(z.design) <- c("baseline", "ER", "PR", "HER2", "grade",
                          "TN_interaction")
  idx <- which(data$design_cat==0)
  model <- EMmvpolySelfDesignnew(y=y[idx,],
                              z.design = z.design,
                              x.self.design = cbind(agefftp,
                                                    breastmos,
                                                    parity,
                                                    refage)[idx,],
                              z.design2 = z.ER,
                              x.self.design2 = study_mat1,
                              missingTumorIndicator = 888)
  save(model,file=paste0("./risk_factor/result/model",i1,".Rdata"))
}else if(i1==9){
  idx <- which(data$design_cat==0)
  z.additive = cbind(1,z.standard)
  model <- EMmvpolySelfDesignnew(y=y[idx,],
                                 z.design = z.additive,
                                 x.self.design = cbind(agefftp,
                                                       breastmos,
                                                       parity,
                                                       refage)[idx,],
                                 z.design2 = z.ER,
                                 x.self.design2 = study_mat2,
                                 missingTumorIndicator = 888
  )
  save(model,file=paste0("./risk_factor/result/model",i1,".Rdata"))
  
}else if(i1==10){
  z.standard <- GenerateZstandard(y)
  ############find triple negative row
  idx <- which(z.standard[,1]==0&
                 z.standard[,2]==0&
                 z.standard[,3]==0)
  z.triple <- rep(1,nrow(z.standard))
  z.triple[idx] <- 0
  z.design <- cbind(1,z.standard,z.triple)
  colnames(z.design) <- c("baseline", "ER", "PR", "HER2", "grade",
                          "TN_interaction")
  idx <- which(data$design_cat==1)
  model <- EMmvpolySelfDesignnew(y=y[idx,],
                                 z.design = z.design,
                                 x.self.design = cbind(agefftp,
                                                       breastmos,
                                                       parity,
                                                       refage)[idx,],
                                 z.design2 = z.ER,
                                 x.self.design2 = study_mat1,
                                 missingTumorIndicator = 888)
  save(model,file=paste0("./risk_factor/result/model",i1,".Rdata"))
}







# write.xlsx(model[[4]]
#            ,file = "risk_factor_result_110118.xlsx",
#            sheetName = "pb_tr_study_continous",
#            append = T)
# write.xlsx(model[[5]]
#            ,file = "risk_factor_result_110118.xlsx",
#            sheetName = "pb_test_study_continous_tr",
#            append=T)
# 
# 
# 
# 
# ############non population based study
# idx <- which(data$design_cat==1)
# ###########two-stage model based on non-population based study
# model <- TwoStageModel(y=y[idx,],
#                        baselineonly = study_mat2,
#                        additive = cbind(agefftp,
#                                         breastmos,
#                                         parity,
#                                         refage)[idx,],
#                        missingTumorIndicator = 888
# )
# write.xlsx(model[[4]]
#            ,file = "risk_factor_result_110118.xlsx",
#            sheetName = "npb_ss_continous",
#            append = T)
# write.xlsx(model[[5]]
#            ,file = "risk_factor_result_110118.xlsx",
#            sheetName = "npb_test_continous",
#            append=T)
# 
# 
# z.standard <- GenerateZstandard(y)
# ############find triple negative row
# idx <- which(z.standard[,1]==0&
#                z.standard[,2]==0&
#                z.standard[,3]==0)
# z.triple <- rep(1,nrow(z.standard))
# z.triple[idx] <- 0
# z.design <- cbind(1,z.standard,z.triple)
# colnames(z.design) <- c("baseline", "ER", "PR", "HER2", "grade",
#                         "TN_interaction")
# 
# 
# idx <- which(data$design_cat==1)
# model <- EMmvpolySelfDesign(y=y[idx,],
#                             z.design = z.design,
#                             x.self.design = cbind(agefftp,
#                                                   breastmos,
#                                                   parity,
#                                                   refage)[idx,],
#                             baselineonly = study_mat2,
#                             missingTumorIndicator = 888)
# 
# write.xlsx(model[[4]]
#            ,file = "risk_factor_result_110118.xlsx",
#            sheetName = "npb_tr_continous",
#            append = T)
# write.xlsx(model[[5]]
#            ,file = "risk_factor_result_110118.xlsx",
#            sheetName = "npb_test_continous_tr",
#            append=T)
# 
# 
