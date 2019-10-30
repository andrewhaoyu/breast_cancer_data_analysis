i1 = 1
library(devtools)
library(data.table)
#install_github("andrewhaoyu/bc2", ref = "development",args = c('--library="/home/zhangh24/R/x86_64-pc-linux-gnu-library/3.4"'))
library(bc2)
setwd("/data/zhangh24/breast_cancer_data_analysis/")
load(paste0("./risk_prediction/result/split.id.rdata"))
#icog.test.id <- Generatetestid(subtypes.icog)
icog.train.id <- split.id[[1]]
#onco.train.id <- split.id[[2]]
#onco.test.id <- split.id[[3]]
#icog.cohort.id <- split.id[[4]]
#onco.cohort.id <- split.id[[5]]
#Icog.order <- read.table(gzfile(subject.file))

data1 <- as.data.frame(fread("/data/zhangh24/breast_cancer_data_analysis/data/sig_snp_icog_prs.csv",header=T))
data1 <- as.data.frame(data1[,-1])
icog.train <- which(data1[,1]%in%icog.train.id)
data1 <- data1[icog.train,]
y.pheno.mis1 <- cbind(data1$Behavior,data1$ER,data1$PR,data1$HER2,data1$Grade)
##########clean phenotype file
idx <- which(y.pheno.mis1[,1]==888|y.pheno.mis1[,1]==2)
y.pheno.mis1[idx,1] <- 1
colnames(y.pheno.mis1) = c("Behavior","ER","PR","HER2","Grade")
#x.test.all.mis1 <- data1[,c(27:206)]

x.covar.mis1 <- data1[,c(7:16)]

#x.all.mis1 <- as.matrix(cbind(x.test.all.mis1[,i1],x.covar.mis1))
#colnames(x.all.mis1)[1] <- "gene"

score.test.support.icog.ERPRHER2Grade <- ScoreTestSupport(
  y.pheno.mis1,
  baselineonly = NULL,
  additive = x.covar.mis1,
  pairwise.interaction = NULL,
  saturated = NULL,
  missingTumorIndicator = 888
)

save(score.test.support.icog.ERPRHER2Grade,file="./risk_prediction/FTOP_whole_genome/ICOG/result/score.test.support.icog.ERPRHER2Grade.Rdata")


