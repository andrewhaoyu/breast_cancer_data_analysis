i1 = 1
library(devtools)
install_github("andrewhaoyu/bc2", ref = "development",args = c('--library="/home/zhangh24/R/x86_64-pc-linux-gnu-library/3.4"'))
library(bc2)
setwd("/data/zhangh24/breast_cancer_data_analysis/")
data1 <- read.csv("./data/iCOGS_euro_v10_05242017.csv",header=T)
y.pheno.mis1 <- cbind(data1$Behaviour1,data1$PR_status1,data1$ER_status1,data1$HER2_status1)
colnames(y.pheno.mis1) = c("Behavior","PR","ER","HER2")
# Grade1.fake <- data1$Grade1
# Grade1.fake[data1$Grade1==2|data1$Grade1==3] <- 1
# Grade1.fake[data1$Grade1==1] <- 0
#y.pheno.mis1 <- cbind(data1$Behaviour1,data1$PR_status1,data1$ER_status1,data1$HER2_status1,Grade1.fake)
# y.pheno.mis1 <- cbind(data1$Behaviour1,data1$PR_status1,data1$ER_status1,data1$HER2_status1)

x.test.all.mis1 <- data1[,c(27:206)]

x.covar.mis1 <- data1[,5:14]
x.all.mis1 <- as.matrix(cbind(x.test.all.mis1[,i1],x.covar.mis1))
colnames(x.all.mis1)[1] <- "gene"

score.test.support.icog.ERPRHER2 <- ScoreTestSupport(
  y.pheno.mis1,
    baselineonly = NULL,
  additive = x.all.mis1[,2:11],
  pairwise.interaction = NULL,
  saturated = NULL,
  missingTumorIndicator = 888
)

save(score.test.support.icog.ERPRHER2,file="./whole_genome/ICOG/ERPRHER2_fixed/result/score.test.support.icog.ERPRHER2.Rdata")

##create pheno.file 

pheno <- data1[,c(1,5:14)]
pheno <- cbind(pheno,data1$Behaviour1,data1$PR_status1,
               data1$ER_status1,data1$HER2_status1,
               data1$Grade1)
colnames(pheno)[12:16] = c("Behaviour1","PR_status1",
                           "ER_status1","HER2_status1",
                           "Grade1")
save(pheno,file="./data/pheno.Icog")
