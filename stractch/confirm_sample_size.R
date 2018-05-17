setwd("/spin1/users/zhangh24/breast_cancer_data_analysis/")
library(data.table)
data2 <- fread("./data/Onco_euro_v10_10232017.csv",header=T)
data2 <- as.data.frame(data2)

data1 <- fread("./data/iCOGS_euro_v10_10232017.csv",header=T)
data1 <- as.data.frame(data1)

new.data1 <- cbind(data1$Behaviour1,data1$ER_status1,data1$PR_status1,data1$HER2_status1,data1$Grade1,data1$age)
new.data2 <- cbind(data2$Behaviour1,data2$ER_status1,data2$PR_status1,data2$HER2_status1,data2$Grade1,data2$age)
new.data <- data.frame(rbind(new.data1,new.data2))
subtypes1 <- GenerateIntrinsic_new(new.data[,1],
                                   new.data[,2],
                                   new.data[,3],
                                   new.data[,4],
                                   new.data[,5])
subtypes <- subtypes1[[1]]
age <- new.data[-subtypes1[[2]],6]
age.cut <- cut(age,breaks=c(15,39.99,49.99,59.99,69.99,100,999))
write.csv(table(age.cut,subtypes),file="./data/intrinsic_sample_size.csv")
GenerateIntrinsic_new <- function(casecon,ER,PR,HER2,Grade){
  n <- length(ER)
  idx.LA <- which(HER2==0&(ER==1|PR==1)&Grade!=3)
  idx.LB <- which(HER2==1&(ER==1|PR==1))
  idx.LUBHER2 <- which(HER2==0&(ER==1|PR==1)&Grade==3)
  idx.HER2 <- which(HER2==1&ER==0&PR==0)
  idx.Tp <- which(HER2==0&ER==0&PR==0)
  subtypes <- rep("control",n)
  idx <- which(casecon==1)
  subtypes[idx] = "missing"
  subtypes[idx.LA] <- "Luminal_A"
  subtypes[idx.LB] <- "Luminal_B"
  subtypes[idx.LUBHER2] <- "Luminal_B_HER2Enriched"
  subtypes[idx.HER2] <- "HER2Enriched"
  subtypes[idx.Tp] <- "TripleNeg"
  idx.remove <- which(subtypes=="missing")
  subtypes <- subtypes[-idx.remove]
  subtypes <- factor(subtypes,levels=c("control",
                                       "Luminal_A",
                                       "Luminal_B",
                                       "Luminal_B_HER2Enriched",
                                       "HER2Enriched",
                                       "TripleNeg"))
  return(list(subtypes,idx.remove))
}

