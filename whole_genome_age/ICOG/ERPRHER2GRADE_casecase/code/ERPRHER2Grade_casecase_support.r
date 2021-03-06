i1 = 1
library(devtools)
library(data.table)
#install_github("andrewhaoyu/bc2", ref = "development",args = c('--library="/home/zhangh24/R/x86_64-pc-linux-gnu-library/3.4"'))
library(bc2)
setwd("/spin1/users/zhangh24/breast_cancer_data_analysis/")
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


Heter.result.Icog = EMmvpoly(y.pheno.mis1,baselineonly = NULL,additive = x.covar.mis1,pairwise.interaction = NULL,saturated = NULL,missingTumorIndicator = 888)
number.of.tumor <- 
z.standard <- Heter.result.Icog[[12]]
delta0 <- Heter.result.Icog[[1]]
M <- 23
delta0 <- c(delta0[1:M],rep(0,2),delta0[(M+1):length(delta0)])
save(delta0,file="./whole_genome_age/ICOG/ERPRHER2GRADE_fixed_baseline/result/delta0.icog.Rdata")
#save(z.standard,file="./whole_genome_age/ICOG/ERPRHER2GRADE_fixed_baseline/result/z.standard.Rdata")












