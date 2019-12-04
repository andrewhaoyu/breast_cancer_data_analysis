idxi1 = 1
library(devtools)
#install_github("andrewhaoyu/bc2", ref = "development",args = c('--library="/home/zhangh24/R/x86_64-pc-linux-gnu-library/3.4"'))
library(bc2, lib.loc ="/home/zhangh24/R/x86_64-pc-linux-gnu-library/3.6/")
library(data.table)
setwd("/data/zhangh24/breast_cancer_data_analysis/")
load(paste0("./risk_prediction/result/split.id.rdata"))
#icog.test.id <- Generatetestid(subtypes.icog)
#icog.train.id <- split.id[[1]]
onco.train.id <- split.id[[2]]
#onco.test.id <- split.id[[3]]
#icog.cohort.id <- split.id[[4]]
#onco.cohort.id <- split.id[[5]]

data2 <- as.data.frame(fread("/data/zhangh24/breast_cancer_data_analysis/data/sig_snp_onco_prs.csv",header=T))
data2 <- data2[,-1]
#onco.train <- which(data2[,1]%in%onco.train.id)
onco.train <- which(data2[,1]%in%onco.train.id)
data2 <- data2[onco.train,]
y.pheno.mis2 <- cbind(data2$Behavior,data2$ER,data2$PR,data2$HER2,data2$Grade)
#######clean y.phneo.mis2 
idx <- which(y.pheno.mis2[,1]==888|y.pheno.mis2[,1]==2)
y.pheno.mis2[idx,1] <- 1
#y.pheno.mis2 <- cbind(data2$Behaviour1,data2$PR_status1,data2$ER_status1,data2$HER2_status1)
colnames(y.pheno.mis2) = c("Behaviour","ER",
                           "PR","HER2","Grade")

x.covar.mis2 <- data2[,c(7:16)]
#ages <- data2[,230]
#idx.complete <- which(ages!=888)

#y.pheno.mis2 <- y.pheno.mis2[idx.complete,]
#x.covar.mis2 <- x.covar.mis2[idx.complete,]
Onc_ID <- data2$ID
#Onc_ID <- Onc_ID[idx.complete]



score.test.support.onco.ERPRHER2Grade <- ScoreTestSupportMixedModel(
  y.pheno.mis2,
  baselineonly = NULL,
  additive = x.covar.mis2,
  pairwise.interaction = NULL,
  saturated = NULL,
  missingTumorIndicator = 888
)

save(score.test.support.onco.ERPRHER2Grade,file="./risk_prediction/FTOP_whole_genome/ONCO/result/score.test.support.onco.ERPRHER2Grade.Rdata")

