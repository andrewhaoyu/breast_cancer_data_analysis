#-------------------------------------------------------------------
# Update Date: 11/14/2018
# Create Date: 11/13/2018
# Goal: prepare sample dataset for qctool
# Author: Haoyu Zhang
#-------------------------------------------------------------------
#prepare sample for all of the subjects
library(data.table)
subject.file <- "/spin1/users/zhangh24/test/onco_order.txt"

onco.order <- read.table(subject.file)
n <- nrow(onco.order)
ID <- matrix("c",n,1)
for(i in 1:n){
  ID[i] <- paste0("sample_",onco.order[i,1])
}
missing <- matrix(0,n,1)
case <- matrix(rbinom(n,1,0.5),n,1)
cov_1 <- matrix(rnorm(n),n,1)
#onco.order <- matrix(paste0("sample",onco.order),n,1)
ID <- rbind(0,ID)
missing <- rbind(0,missing)
case <- rbind("B",case)
cov_1 <- rbind("C",cov_1)
#onco.order <- matrix(onco.order,ncol=1)

sample.data <- data.frame(ID,
                          ID,
                          missing,
                          case,
                          cov_1,stringsAsFactors = T)
colnames(sample.data) <- c("ID_1",
                           "ID_2",
                           "missing",
                           "case",
                           "cov_1")

write.table(sample.data,file = "/spin1/users/zhangh24/test/sample.txt",
            row.names = F, quote = F,sep = " ")

#prepare sample for people keeping in test
setwd('/spin1/users/zhangh24/breast_cancer_data_analysis/')
load(paste0("./risk_prediction/result/split.id.rdata"))
#icog.test.id <- Generatetestid(subtypes.icog)
icog.train.id <- split.id[[1]]
onco.train.id <- split.id[[2]]
onco.test.id <- split.id[[3]]
n <- length(onco.test.id)
test_ID <- matrix("c",n,1)
for(i in 1:n){
  test_ID[i] <- paste0("sample_",onco.test.id[i])
}

