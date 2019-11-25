#goal: get the auc of the 313 Nasim's SNPs

setwd('/data/zhangh24/breast_cancer_data_analysis/')
library(dplyr)
library(data.table)
library(pROC)
select.names <- c("Luminal_A",
                  "Luminal_B",
                  "Luminal_B_HER2Neg",
                  "HER2Enriched",
                  "TripleNeg")

load(paste0("./risk_prediction/result/split.id.rdata"))
#icog.test.id <- Generatetestid(subtypes.icog)
#load the sample data
icog.train.id <- split.id[[1]]
onco.train.id <- split.id[[2]]
onco.test.id <- split.id[[3]]
#onco.test.id <- split.id[[5]]
icog.vad.id <- split.id[[4]]
onco.vad.id <- split.id[[5]]


load("/data/zhangh24/breast_cancer_data_analysis/risk_prediction/Nasim_prs/result/onco.nasim.snp.rdata")
data2 <- as.data.frame(fread("/data/zhangh24/breast_cancer_data_analysis/data/sig_snp_onco_prs.csv",header=T))
data2 <- data2[,-1]

#order data2 and onco.nasim.snp to be the same order
ID_data2 <- data2[,1,drop=F]

colnames(onco.nasim.snp)
onco.nasim.snp.order <- left_join(ID_data2,onco.nasim.snp,by="ID")


onco.nasim.snpvalue <- onco.nasim.snp.order[,2:ncol(onco.nasim.snp.order)]



library(bcutility,lib.loc = "/spin1/home/linux/zhangh24/R/x86_64-pc-linux-gnu-library/3.6/")
y.pheno.mis2 <- cbind(data2$Behavior,data2$ER,data2$PR,data2$HER2,data2$Grade)
subtypes.onco <- GenerateIntrinsicmis(y.pheno.mis2[,2],
                                                  y.pheno.mis2[,3],
                                                  y.pheno.mis2[,4],
                                                  y.pheno.mis2[,5])
#load the odds ratio
load("./risk_prediction/Nasim_prs/result/313_intrinsic_subtype_logodds.Rdata")
log.odds <- final_result[,4:ncol(final_result)]
prs <- as.matrix(onco.nasim.snpvalue)%*%as.matrix(log.odds)




onco.test <- which(data2[,1]%in%onco.test.id)
y.test <- data2$Behavior[onco.test]
subtypes.test <- as.character(subtypes.onco[onco.test])
prs.test <- prs[onco.test,]
total <- length(select.names)
auc <- rep(0,total)
auc95 <- rep("c",total)
subtypes <- rep("c",total)
ind <- 1

for(i in 1:5){
  target.subtype <- select.names[i]
  idx <- which(subtypes.test==target.subtype|
                 subtypes.test=="control")
  y.test.tar <- y.test[idx]
  prs.test.tar <- prs.test[idx,i]
  roc.standard <- roc(y.test.tar,as.vector(prs.test.tar),ci=T,plot=F)
  pla <- 4
  auc[ind] <- round(as.numeric(roc.standard$auc),pla)*100
  
  auc95[ind] <- paste0(round(as.numeric(roc.standard$auc)[1],pla)*100,
                       "(",
                       round(as.numeric(roc.standard$ci)[1],pla)*100,
                       "-",
                       round(as.numeric(roc.standard$ci)[3],pla)*100
                       ,")")
  subtypes[ind] <- target.subtype
  ind <- ind+1
}

auc.result <- data.frame(auc,
                         auc95,
                         subtypes)

