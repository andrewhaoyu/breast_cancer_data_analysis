args = commandArgs(trailingOnly = T)
i1 = as.numeric(args[[1]])
library(devtools)
############install the development R package bc2
############the repository is called bc3, but the package is called bc2
install_github("andrewhaoyu/bc3")

library(data.table)
library(bc2)
#data <- fread("./data/dataset_montse_20180522.txt")
setwd('/dcl01/chatterj/data/hzhang1/breast_cancer_data_analysis')
#setwd('/Users/haoyuzhang/GoogleDrive/breast_cancer_data_analysis (1)')
data <- as.data.frame(fread("./data/Dataset_Montse_2018-10-10.txt"))
##############we only focus on the invasive breast cancer cases
idx.invasive <- which(data$status==0|data$status==1)
data <- data[idx.invasive,]

idx <- which(data$status==0)

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


##########create study mat
idx <- which(data$design_cat==0)
data1 <- data[idx,]
study_mat1 <- model.matrix(~as.factor(data1$study)-1)[,-1]
idx <- which(data$design_cat==1)
data2 <- data[idx,]
study_mat2 <- model.matrix(~as.factor(data2$study)-1)[,-1]


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
z.standard <- GenerateZstandard(y)
z.ER <- cbind(1,z.standard[,1])
if(i1==1){
  ############don't adjust for study
  ############population based study
  idx <- which(data$design_cat==0)
  z.additive <- cbind(1,z.standard)
  ###########two-stage model based on population based study
  model <- EMmvpolySelfDesignnew(y=y[idx,],
                               x.self.design = cbind(agefftp_mat,
                                                     breast_mat,
                                                     parity_mat,
                                                     refage)[idx,],
                               z.design = z.additive,
                               x.self.design2 = study_mat1,
                               z.design2 = z.ER,
                         missingTumorIndicator = 888
  )
save(model,file=paste0("./risk_factor/result/model",i1,".Rdata"))
    # write.xlsx(model[[4]]
  #            ,file = "risk_factor_result_110118.xlsx",
  #            sheetName = "population_based_second_stage")
  # write.xlsx(model[[5]]
  #            ,file = "risk_factor_result_110118.xlsx",
  #            sheetName = "population_based_test",
  #            append=T)
  
}else if(i1 ==2){
  ############two-stage model with triple negative interaction
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
                                 x.self.design = cbind(agefftp_mat,
                                                       breast_mat,
                                                       parity_mat,
                                                       refage)[idx,],
                              z.design = z.design,
                             
                              z.design2 = z.ER,
                              x.self.design2 = study_mat1,
                              missingTumorIndicator = 888)
  save(model,file=paste0("./risk_factor/result/model",i1,".Rdata"))
  # write.xlsx(model[[4]]
  #            ,file = "risk_factor_result_110118.xlsx",
  #            sheetName = "population_based_second_stage_triple",
  #            append = T)
  # write.xlsx(model[[5]]
  #            ,file = "risk_factor_result_110118.xlsx",
  #            sheetName = "population_based_test_triple",
  #            append=T)
  
}else if(i1==3){
  ############two-stage model with intrinsic subtypes
  z.design <- matrix(c(
    c(0,1,1,1,0,0,0,0,1,1,1,0,0,0,0,0,1,1,0,0,0,0,0),
    c(0,0,0,0,0,1,1,0,0,0,0,0,1,1,1,0,0,0,0,0,1,1,1),
    c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0),
    c(0,0,0,0,1,0,0,0,0,0,0,1,0,0,0,0,0,0,0,1,0,0,0),
    c(1,0,0,0,0,0,0,1,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0)
  ),ncol=5)
  rowSums(z.design)
  
  
  
  colnames(z.design) <- c("Luminal A","Luminal B",
                          "Luminal B HER2Neg",
                          "HER2 Enriched",
                          "Triple Negative")
  
  
  idx <- which(data$design_cat==0)
  model <- EMmvpolySelfDesignnew(y=y[idx,],
                              z.design = z.design,
                              x.self.design = cbind(agefftp_mat,
                                                    breast_mat,
                                                    parity_mat,
                                                    refage)[idx,],
                              z.design2 = z.ER,
                              x.self.design2 = study_mat1,
                              missingTumorIndicator = 888)
  save(model,file=paste0("./risk_factor/result/model",i1,".Rdata"))
  # write.xlsx(model[[4]]
  #            ,file = "risk_factor_result_110118.xlsx",
  #            sheetName = "pb_second_stage_in",
  #            append = T)
  # 
}else if(i1 ==4){
  idx <- which(data$design_cat==1)
  z.additive <- cbind(1,z.standard)
  ###########two-stage model based on non population based study
  model <- EMmvpolySelfDesignnew(y=y[idx,],
                                 x.self.design = cbind(agefftp_mat,
                                                       breast_mat,
                                                       parity_mat,
                                                       refage)[idx,],
                                 z.design = z.additive,
                                 x.self.design2 = study_mat2,
                                 z.design2 = z.ER,
                                 missingTumorIndicator = 888
  )
  save(model,file=paste0("./risk_factor/result/model",i1,".Rdata"))
}else if(i1==5){
  ############two-stage model with triple negative interaction
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
                                 x.self.design = cbind(agefftp_mat,
                                                       breast_mat,
                                                       parity_mat,
                                                       refage)[idx,],
                                 z.design = z.design,
                                 
                                 z.design2 = z.ER,
                                 x.self.design2 = study_mat2,
                                 missingTumorIndicator = 888)
  save(model,file=paste0("./risk_factor/result/model",i1,".Rdata"))
  
}else if(i1==6){
  ############two-stage model with intrinsic subtypes
  z.design <- matrix(c(
    c(0,1,1,1,0,0,0,0,1,1,1,0,0,0,0,0,1,1,0,0,0,0,0),
    c(0,0,0,0,0,1,1,0,0,0,0,0,1,1,1,0,0,0,0,0,1,1,1),
    c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0),
    c(0,0,0,0,1,0,0,0,0,0,0,1,0,0,0,0,0,0,0,1,0,0,0),
    c(1,0,0,0,0,0,0,1,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0)
  ),ncol=5)
  rowSums(z.design)
  
  
  
  colnames(z.design) <- c("Luminal A","Luminal B",
                          "Luminal B HER2Neg",
                          "HER2 Enriched",
                          "Triple Negative")
  
  
  idx <- which(data$design_cat==0)
  model <- EMmvpolySelfDesignnew(y=y[idx,],
                                 z.design = z.design,
                                 x.self.design = cbind(agefftp_mat,
                                                       breast_mat,
                                                       parity_mat,
                                                       refage)[idx,],
                                 z.design2 = z.ER,
                                 x.self.design2 = study_mat1,
                                 missingTumorIndicator = 888)
  save(model,file=paste0("./risk_factor/result/model",i1,".Rdata"))
}
