i1 = 1
library(devtools)
library(data.table)
#install_github("andrewhaoyu/bc2", ref = "development",args = c('--library="/home/zhangh24/R/x86_64-pc-linux-gnu-library/3.4"'))
library(bc2)
setwd("/data/zhangh24/breast_cancer_data_analysis/")
data1 <- fread("./data/iCOGS_euro_v10_10232017.csv",header=T)
data1 <- as.data.frame(data1)
y.pheno.mis1 <- cbind(data1$Behaviour1,data1$ER_status1,data1$PR_status1,data1$HER2_status1,data1$Grade1)
colnames(y.pheno.mis1) = c("Behavior","ER","PR","HER2","Grade")
#x.test.all.mis1 <- data1[,c(27:206)]

x.covar.mis1 <- data1[,c(5:14,204)]
age <- data1[,204]
idx.complete <- which(age!=888)
y.pheno.mis1 <- y.pheno.mis1[idx.complete,]
x.covar.mis1 <- x.covar.mis1[idx.complete,]

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

save(score.test.support.icog.ERPRHER2Grade,file="./whole_genome_age/ICOG/ERPRHER2GRADE_fixed_baseline/result/score.test.support.icog.ERPRHER2Grade.Rdata")

