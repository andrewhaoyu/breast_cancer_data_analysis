setwd('/spin1/users/zhangh24/breast_cancer_data_analysis/')
load('./genetic_correlation/ICOG/result/ICOG.result.Rdata')
load("./genetic_correlation/result/ICOG.result.clean.completeglm.Rdata")
ICOG.result.clean.completeglm <- ICOG.result.clean
head(ICOG.result[[1]])
head(ICOG.result.clean.completeglm)
idx <- which(ICOG.result[[1]][,2]==1&
               ICOG.result[[1]][,3]==100017701)
ICOG.result[[1]][idx,]
ICOG.result[[2]][idx,]
idx <- which(ICOG.result.clean.completeglm[,6]==1&
               ICOG.result.clean.completeglm[,7]==100017701)
ICOG.result.clean.completeglm[idx,11:15]
dim(ICOG.result[[1]])
