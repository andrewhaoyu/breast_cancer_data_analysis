args = commandArgs(trailingOnly = T)
i1 = as.numeric(args[[1]])
library(devtools)
############install the development R package bc2
############the repository is called bc3, but the package is called bc2
#install_github("andrewhaoyu/bc3")

library(data.table)
library(bc3)
#data <- fread("./data/dataset_montse_20180522.txt")
#setwd('/dcl01/chatterj/data/hzhang1/breast_cancer_data_analysis')
setwd('/spin1/users/zhangh24/breast_cancer_data_analysis/')
data <- as.data.frame(fread("./data/Dataset_Montse_20190322.txt"))
##############we only focus on the invasive breast cancer cases
##############collapse breastmost_cat 0 agefftp 0 and lastchildage 0 as 1, otherwise there will be collinearity issue with parity. Since women with no children will also fall into these categoriess
data$status[data$status==2|data$status==3] <- 1
data$breastmos_cat[data$breastmos_cat==0] <- 1
data$agefftp_cat[data$agefftp_cat==0] <- 1
data$lastchildage_cat[data$lastchildage_cat==0] <- 1


##############only focus on population based study


library(dplyr)
data1 = data %>% filter(design_cat==0)

#############put the missing tumor characteristics as 888
idx.ER.mis <- which(data1$status==1&is.na(data1$ER_status1))
data1$ER_status1[idx.ER.mis] <- 888
idx.PR.mis <- which(data1$status==1&is.na(data1$PR_status1))
data1$PR_status1[idx.PR.mis] <- 888
idx.HER2.mis <- which(data1$status==1&is.na(data1$HER2_status1))
data1$HER2_status1[idx.HER2.mis] <- 888
idx.grade.mis <- which(data1$status==1&is.na(data1$Grade1))
data1$Grade1[idx.grade.mis] <- 888
#############put the subject BCAC-16687664 HER2 status as 1
idx <- which(data1$HER2_status1==2)
data1$HER2_status1[idx] <- 1
#############check the result
table(data1$status,data1$ER_status1)
table(data1$status,data1$PR_status1)
table(data1$status,data1$HER2_status1)
table(data1$status,data1$Grade1)
library(nnet)


########putting missing variable as NA
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
                             #as.factor(lastchildage_cat),
                             study,
                             refage)
library(mice)
imp <- mice(all.covariates,m=1,seed=i1,print=FALSE)

all.covariates.c <- complete(imp,1)
new.data1 <- data.frame(data1$molgroup,all.covariates.c)
colnames(new.data1) <- c("molgroup",
                         "agemenarche_cat",
                         "parity_cat",
                         "mensagelast_cat",
                         "agefftp_cat",
                         "breastmos_cat",
                         #"lastchildage_cat",
                         "study",
                         "refage")

###########create the dummy variable for agemenarche_cat
###########create the dummy variable for parity_cat
###########create the dummy variable for mensagelat_cat
###########create the dummy variable for agefftp_cat
###########create the dummy variable for breastmos_cat
###########create the dummy variable for study

agemenarche_mat <- model.matrix(~as.factor(new.data1$agemenarche_cat)-1)[,-1]
colnames(agemenarche_mat) <- paste0("agemenarche_cat",c(1:3))
parity_mat <- model.matrix(~as.factor(new.data1$parity_cat)-1)[,-1]
colnames(parity_mat) <- paste0("parity_cat",c(1:4))
mensagelat_mat <- model.matrix(~as.factor(new.data1$mensagelast_cat)-1)[,-1]
colnames(mensagelat_mat) <- paste0("mensagelat_cat",c(1:2))
agefftp_mat <- model.matrix(~as.factor(new.data1$agefftp_cat)-1)[,-1]
colnames(agefftp_mat) <- paste0("agefftp",c(2:4))
breastmos_mat <- model.matrix(~as.factor(new.data1$breastmos_cat)-1)[,-1]
colnames(breastmos_mat) <- paste0("breastmos_cat",c(2:5))
study_mat <- model.matrix(~as.factor(new.data1$study)-1)[,-1]
refage <- new.data1$refage


###########create the phenotype file
y <- cbind(data1$status,data1$ER_status1,data1$PR_status1,
           data1$HER2_status1,data1$Grade1)
colnames(y) <- c("casecontrol",
                 "ER",
                 "PR",
                 "HER2",
                 "Grade")
# idx <- which(data$status==0)
# y[idx,2:ncol(y)] <- NA
z.standard <- GenerateZstandard(y)
z.ER <- cbind(1,z.standard[,1])
  ############two-stage model with intrinsic subtypes
  z.design <- matrix(c(
    c(0,1,1,1,0,0,0,0,1,1,1,0,0,0,0,0,1,1,0,0,0,0,0),
    c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0),
    c(0,0,0,0,0,1,1,0,0,0,0,0,1,1,1,0,0,0,0,0,1,1,1),
    c(0,0,0,0,1,0,0,0,0,0,0,1,0,0,0,0,0,0,0,1,0,0,0),
    c(1,0,0,0,0,0,0,1,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0)
  ),ncol=5)
  rowSums(z.design)
  
  
  
  colnames(z.design) <- c("Luminal A-like","Luminal B,HER2-negative-like",
                          "Luminal B-like",
                          "HER2 enriched-like",
                          "TN")
  
  #idx <- which(data$design_cat==0)
  model <- EMmvpolySelfDesignnew(y=y,
                                 z.design = z.design,
                                 x.self.design = cbind(
                                   agemenarche_mat,
                                   parity_mat,
                                   mensagelat_mat,
                                   agefftp_mat,
                                                       breastmos_mat,
                                                       refage),
                                 z.design2 = z.ER,
                                 x.self.design2 = study_mat,
                                 missingTumorIndicator = 888)
  save(model,file=paste0("./risk_factor/result/intrinsic_imp",i1,".Rdata"))
