library(data.table)
library(bc2)
#data <- fread("./data/dataset_montse_20180522.txt")
data <- fread("./data/dataset_montse_20180524.txt")


transform <- function(status,tumor){
  n <- length(status)
  result <- tumor
  idx.mis <- which(status==1&is.na(tumor))
  result[idx.mis] <- 888
  idx.control <- which(status==0)
  result[idx.control] <- NA
  return(result)
}
##########remove status 2 in-situ
##########remove status 3 unknown people
idx.remove <- which(data$status==2|data$status==3)
data <- data[-idx.remove,]
#########remove people with missing parity_cat


idx.remove <- which(data$parity_cat==9)
data <- data[-idx.remove,]


status <- data$status
ER <- data$ER_status1
ER.update <- transform(status,ER)
PR <- data$PR_status1
PR.update <- transform(status,PR)
HER2 <- data$HER2_status1
HER2.update <- transform(status,HER2)
Grade <- data$Grade1
Grade.update <- transform(status,Grade)
table(data$parity_cat)
table(data$study)
table(data$ethnicity)
table(data$refage)
y.pheno.mis <- cbind(status,ER.update,PR.update,HER2.update
                     ,Grade.update)
colnames(y.pheno.mis) <- c("casecon","ER","PR","HER","Grade")
parity.mat <- model.matrix(~as.factor(data$parity_cat)-1)
colnames(parity.mat) <- paste0("parity_cat",c(0:4))
study.mat <- model.matrix(~as.factor(data$study)-1)
colnames(study.mat) <- unique(data$study)
ethnicity.mat <- model.matrix(~as.factor(data$ethnicity)-1)
colnames(ethnicity.mat) <- paste0("ethnicity",c(1:6))
refage <- data$refage
x.covar <- cbind(parity.mat,ethnicity.mat,refage)
#study.mat
model <- TwoStageModel(y = y.pheno.mis,baselineonly = study.mat,additive = x.covar,
                       missingTumorIndicator = 888)

y = y.pheno.mis
additive = x.covar


# data.com <- fread("./data/concept_542-543-change-claude_bcac_pheno_v10_SPT_100217.txt")
# 
# 
# data.com.new <- data.com[,c(2,40,48,56,10)]
# 
# 
# new.
# 
# 
# getZstandard()