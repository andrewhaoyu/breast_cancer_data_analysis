#-------------------------------------------------------------------
# Update Date: 11/22/2018
# Create Date: 11/22/2018
# Goal: calculate auc
# Author: Haoyu Zhang
#-------------------------------------------------------------------
#---------------------------------------#---------------------------------------
setwd('/spin1/users/zhangh24/breast_cancer_data_analysis/')
library(dplyr)
library(data.table)
library(pROC)
library(bcutility)
prs <- as.data.frame(fread("/spin1/users/zhangh24/BCAC/prs_out/Lu_standard_prs7.profile",header=T))
sample.data <- as.data.frame(fread("/spin1/users/zhangh24/test/sample.txt"))

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
library(dplyr)
names.subtypes <-  c("Luminal_A","Luminal_B",
                     "Luminal_B_HER2Neg",
                     "HER2Enriched",
                     "TripleNeg")
#test sample Luminal A
#for(i in 1:length(names.subtypes)){
i <- 1 
 test.sample = sample.data%>%
    filter(
      (ID_1%in%test_ID)&
        ((subtypes=="control"|
            subtypes==names.subtypes[i]))
    ) %>%
    select(ID_1,case)
#}

colnames(test.sample) <- c("IID","case")
test.sample <- inner_join(test.sample,prs)
prs.standard <- test.sample[,"SCORE"]
y.test <- test.sample[,"case"]
roc.standard <- calibration(y.test,prs.standard)
roc.standard <- roc(y.test,as.vector(prs.standard),ci=T,plot=F)
roc.standard
pla <- 4
auc.result[1] <- round(as.numeric(roc.standard$auc),pla)*100

auc.95[1] <- paste0(round(as.numeric(roc.standard$ci)[1],pla)*100,
                    "-",
                    round(as.numeric(roc.standard$ci)[3],pla)*100)
cal.standard <- calibration(y.test,prs.standard)







  