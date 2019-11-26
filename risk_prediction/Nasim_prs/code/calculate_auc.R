#calculate auc based on nasim 313 SNPs
setwd("/data/zhangh24/breast_cancer_data_analysis/")
load(paste0("./risk_prediction/Nasim_prs/result/313_intrinsic_subtype_logodds.Rdata"))
intrinsic_subtypes_result <- final_result
load(paste0("./risk_prediction/Nasim_prs/result/313_overall_logodds.Rdata"))
overall_result <- final_result
all_result <- cbind(intrinsic_subtypes_result,
                    overall_result[,4:6])

library(bc2, lib.loc ="/home/zhangh24/R/x86_64-pc-linux-gnu-library/3.6/")

load("/data/zhangh24/breast_cancer_data_analysis/risk_prediction/Nasim_prs/result/onco.nasim.snp.rdata")
snpvalue.result <- onco.nasim.snp[,2:ncol(onco.nasim.snp)]

library(dplyr)
logodds = all_result %>% 
  select(Luminal_A,Luminal_B,Luminal_B_HER2Neg,HER2_Enriched,TN,logodds_overall,logodds_erpos,logodds_erneg)
all.equal(colnames(snpvalue.result),all_result[,3])
prs_all <- as.matrix(snpvalue.result)%*%as.matrix(logodds)

#use overall analysis to predict auc results
library(dplyr)
library(data.table)
library(pROC)
#library(bcutility)
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


total <- ncol(prs_all)*length(select.names)
auc <- rep(0,total)
auc95 <- rep("c",total)
ind <- 1
method <- rep("c",total)
subtypes <- rep("c",total)
for(j in 1:length(select.names)){
  for(i in 1:ncol(prs_all)){
 
    prs <- prs_all[,i]
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
              subtypes==names.subtypes[temp])))%>%
      select(ID_1,case,prs)
    
    
    
    
    
    
    
    colnames(test.sample) <- c("IID","case","prs")
    
    test_test.sample <- test.sample %>% 
      group_by(case) %>% 
      mutate(meanprs = mean(prs))
    
    
    prs.standard <- test.sample[,"prs"]
    y.test <- test.sample[,"case"]
    
    roc.standard <- roc(y.test,as.vector(prs.standard),ci=T,plot=F)
    
    
    pla <- 4
    auc[ind] <- round(as.numeric(roc.standard$auc),pla)*100
    
    auc95[ind] <- paste0(round(as.numeric(roc.standard$auc)[1],pla)*100,
                         "(",
                         round(as.numeric(roc.standard$ci)[1],pla)*100,
                         "-",
                         round(as.numeric(roc.standard$ci)[3],pla)*100
                         ,")")
    
    subtypes[ind] <- names.subtypes[temp]
    method[ind] <- colnames(prs_all)[i]
    ind = ind + 1
    
  }
}

auc.result <- data.frame(auc,
                         auc95,
                         subtypes,
                         method)
write.csv(auc.result,file = paste0("./risk_prediction/Nasim_prs/result/nasim_snp_auc_result_test_data.csv"))













#use validation dataset to verify
onco.test.id <- split.id[[5]]



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


total <- ncol(prs_all)*length(select.names)
auc <- rep(0,total)
auc95 <- rep("c",total)
ind <- 1
method <- rep("c",total)
subtypes <- rep("c",total)
for(j in 1:length(select.names)){
  for(i in 1:ncol(prs_all)){
    
    prs <- prs_all[,i]
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
              subtypes==names.subtypes[temp])))%>%
      select(ID_1,case,prs)
    
    vad_test.sample <- test.sample %>% 
      group_by(case) %>% 
      mutate(meanprs = mean(prs))
    
    
    
    
    
    
    colnames(test.sample) <- c("IID","case","prs")
    
    
    
    prs.standard <- test.sample[,"prs"]
    y.test <- test.sample[,"case"]
    
    roc.standard <- roc(y.test,as.vector(prs.standard),ci=T,plot=F)
    
    
    pla <- 4
    auc[ind] <- round(as.numeric(roc.standard$auc),pla)*100
    
    auc95[ind] <- paste0(round(as.numeric(roc.standard$auc)[1],pla)*100,
                         "(",
                         round(as.numeric(roc.standard$ci)[1],pla)*100,
                         "-",
                         round(as.numeric(roc.standard$ci)[3],pla)*100
                         ,")")
    
    subtypes[ind] <- names.subtypes[temp]
    method[ind] <- colnames(prs_all)[i]
    ind = ind + 1
    
  }
}

auc.result <- data.frame(auc,
                         auc95,
                         subtypes,
                         method)
write.csv(auc.result,file = paste0("./risk_prediction/Nasim_prs/result/nasim_snp_auc_result_vad_data.csv"))


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


total <- ncol(prs_all)*length(select.names)
auc <- rep(0,total)
auc95 <- rep("c",total)
ind <- 1
method <- rep("c",total)
subtypes <- rep("c",total)
for(j in 1:length(select.names)){
  for(i in 1:ncol(prs_all)){
    
    prs <- prs_all[,i]
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
              subtypes==names.subtypes[temp])))%>%
      select(ID_1,case,prs)
    
    
    
    
    
    
    
    colnames(test.sample) <- c("IID","case","prs")
    
    testvad_test.sample <- test.sample %>% 
      group_by(case) %>% 
      mutate(meanprs = mean(prs))
    
    
    
    prs.standard <- test.sample[,"prs"]
    y.test <- test.sample[,"case"]
    
    roc.standard <- roc(y.test,as.vector(prs.standard),ci=T,plot=F)
    
    
    pla <- 4
    auc[ind] <- round(as.numeric(roc.standard$auc),pla)*100
    
    auc95[ind] <- paste0(round(as.numeric(roc.standard$auc)[1],pla)*100,
                         "(",
                         round(as.numeric(roc.standard$ci)[1],pla)*100,
                         "-",
                         round(as.numeric(roc.standard$ci)[3],pla)*100
                         ,")")
    
    subtypes[ind] <- names.subtypes[temp]
    method[ind] <- colnames(prs_all)[i]
    ind = ind + 1
    
  }
}

auc.result <- data.frame(auc,
                         auc95,
                         subtypes,
                         method)
write.csv(auc.result,file = paste0("./risk_prediction/Nasim_prs/result/nasim_snp_auc_result_testvad_data.csv"))



vad_test.sample[idx,]
head(vad_test.sample)
head(test_test.sample)
tail(test_test.sample)
head(testvad_test.sample)
tail(testvad_test.sample)
