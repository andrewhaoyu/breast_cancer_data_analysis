#calculate the risk ratio for five subtypes based on overall and ER specific PRS
library(devtools)
library(withr)
setwd("/data/zhangh24/breast_cancer_data_analysis/")
#filter all the discovery SNPs +-500kb of the 313 SNPs or r2>=0.1
#read in nasim_snp_infor
snp.nasim <- read.csv("./data/Nasim_313_SNPs_infor.csv",header=T)
library(MASS)
colnames(snp.nasim)[c(2,3)] <- c("CHR","position")
#load in nasim's oncoarray genotype
load("/data/zhangh24/breast_cancer_data_analysis/risk_prediction/Nasim_prs/result/onco.nasim.snp.rdata")
onco.nasim.genotype <- onco.nasim.snp[,2:ncol(onco.nasim.snp)]

load(paste0("./risk_prediction/Nasim_prs/result/313_overall_logodds.Rdata"))
overall_result <- final_result


library(bc2, lib.loc ="/home/zhangh24/R/x86_64-pc-linux-gnu-library/3.6/")


snpvalue.result <- onco.nasim.genotype



library(dplyr)
logodds = overall_result %>% 
  select(logodds_overall,logodds_erpos,logodds_erneg)
all.equal(colnames(snpvalue.result),overall_result[,3])
prs_all <- as.matrix(snpvalue.result)%*%as.matrix(logodds)

#subtypes names in the files, there is a difference, one is called HER2_Enriched, the other one is called HER2Enriched
select.names <- c("Luminal_A",
                  "Luminal_B",
                  "Luminal_B_HER2Neg",
                  "HER2Enriched",
                  "TN")
#load all the datasets
#the subtypes names in sample data
names.subtypes <- c("Luminal_A",
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
icog.vad.id <- split.id[[4]]
onco.vad.id <- split.id[[5]]

#combine test and validation
onco.test.id <- c(split.id[[3]],split.id[[5]])

sample.data <- read.table("/data/zhangh24/BCAC/impute_onco/sample.txt",header=T,stringsAsFactors = F)
sample.data <- sample.data[-1,,drop=F]

onco.test.id <- as.vector(onco.test.id)
n <- length(onco.test.id)
test_ID <- rep("c",n)
for(i in 1:n){
  #fix the issue that prs file will code 100000 as 1e+05
  if(onco.test.id[i]==100000){
    test_ID[i] <- paste0("sample_1e+05")  
  }else{
    test_ID[i] <- paste0("sample_",as.numeric(onco.test.id[i]))
  }
  
}
#Get the odds ratio for one subtypes based on the sample.data and prs file
GetOddsRatio <- function(test.sample,quantile.thres){
  odds.ratio <- rep(0,length(quantile.thres)-1)
  
  prs.standard <- test.sample[,"prs"]
  y.test <- test.sample[,"case"]
  quan.set <-quantile(prs.standard,quantile.thres)  
  temp <- 1
  odds <- rep(0,length(quan.set)-1)
  for(k in 1:(length(quan.set)-1)){
    idx <- which(prs.standard>=quan.set[temp]&prs.standard<=quan.set[temp+1])
    
    p <- sum(as.numeric(y.test[idx]))/length(idx)
    
    
    odds[temp] <- p/(1-p)
    
    temp <- temp+1
  }  
  for(k in 1:(length(quan.set)-1)){
    odds.ratio[k] <- odds[k]/odds[6]
  }
  return(odds.ratio)
}






quantile.thres <- c(0,0.01,0.05,0.10,
                    0.20,0.40,0.60,0.80,0.90,0.95,0.99,1) 
odds.ratio.overall <- matrix(0,length(quantile.thres)-1,
                     length(select.names) )


total <- length(select.names)
for(j in 1:length(select.names)){
  
  #first for overall prs results
  prs <- prs_all[,'logodds_overall']
  #prs[,4] <- prs.la.temp
  temp <- j%%5
  if(temp==0){temp=5}
  sample.data.all <- cbind(sample.data,prs)
  #select the sample in the test data
  #also select sample for control and target subtypes
  test.sample = sample.data.all%>%
    filter(
      (ID_1%in%test_ID)&
        ((subtypes=="control"|
            subtypes==names.subtypes[temp])))%>%select(ID_1,case,prs)
  colnames(test.sample) <- c("IID","case","prs")
  odds.ratio.overall[,j] <- GetOddsRatio(test.sample,quantile.thres)
  
 
}
colnames(odds.ratio.overall) <- select.names
write.csv(odds.ratio.overall,file = paste0("./risk_prediction/Nasim_prs/result/odds.ratio.overall.csv"))




odds.ratio.er <- matrix(0,length(quantile.thres)-1,
                             length(select.names) )


total <- length(select.names)
for(j in 1:length(select.names)){
  
  #first three subtypes are er positive subtypes
  if(j<=3){
    prs <- prs_all[,'logodds_erpos']
    
  }else{
    #the other two subtypes are er negative subtypes
    prs <- prs_all[,'logodds_erneg']
    
  }
  #prs[,4] <- prs.la.temp
  temp <- j%%5
  if(temp==0){temp=5}
  sample.data.all <- cbind(sample.data,prs)
  #select the sample in the test data
  #also select sample for control and target subtypes
  test.sample = sample.data.all%>%
    filter(
      (ID_1%in%test_ID)&
        ((subtypes=="control"|
            subtypes==names.subtypes[temp])))%>%select(ID_1,case,prs)
  colnames(test.sample) <- c("IID","case","prs")
  odds.ratio.er[,j] <- GetOddsRatio(test.sample,quantile.thres)
  
  
}
colnames(odds.ratio.er) <- select.names
write.csv(odds.ratio.er,file = paste0("./risk_prediction/Nasim_prs/result/odds.ratio.er.csv"))


