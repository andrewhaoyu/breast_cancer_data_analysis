i1 = 1
library(devtools)
library(data.table)
#install_github("andrewhaoyu/bc2", ref = "development",args = c('--library="/home/zhangh24/R/x86_64-pc-linux-gnu-library/3.4"'))
library(bc2)
setwd("/data/zhangh24/breast_cancer_data_analysis/")
data2 <- fread("./data/Onco_euro_v10_10232017.csv",header=T)
data2 <- as.data.frame(data2)
y.pheno.mis2 <- cbind(data2$Behaviour1,data2$ER_status1,data2$PR_status1,data2$HER2_status1,data2$Grade1)
#y.pheno.mis2 <- cbind(data2$Behaviour1,data2$PR_status1,data2$ER_status1,data2$HER2_status1)
colnames(y.pheno.mis2) = c("Behaviour","ER",
                           "PR","HER2","Grade")

x.covar.mis2 <- data2[,c(5:14,204)]
ages <- data2[,204]
idx.complete <- which(ages!=888)

y.pheno.mis2 <- y.pheno.mis2[idx.complete,]
x.covar.mis2 <- x.covar.mis2[idx.complete,]
#x.all.mis1 <- as.matrix(cbind(x.test.all.mis1[,i1],x.covar.mis1))
#colnames(x.all.mis1)[1] <- "gene"


Heter.result.onco = EMmvpoly(y.pheno.mis2,baselineonly = NULL,additive = x.covar.mis2,pairwise.interaction = NULL,saturated = NULL,missingTumorIndicator = 888)

  z.standard <- Heter.result.onco[[12]]
delta0 <- Heter.result.onco[[1]]
M <- 23
delta0 <- c(delta0[1:M],rep(0,2),delta0[(M+1):length(delta0)])
save(delta0,file="./whole_genome_age/ONCO/ERPRHER2GRADE_fixed_baseline/result/delta0.onco.Rdata")
save(z.standard,file="./whole_genome_age/ONCO/ERPRHER2GRADE_fixed_baseline/result/z.standard.Rdata")












