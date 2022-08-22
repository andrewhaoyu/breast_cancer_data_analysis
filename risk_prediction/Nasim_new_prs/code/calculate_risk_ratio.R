#calculate the risk ratio for five subtypes 
setwd("/data/zhangh24/breast_cancer_data_analysis/")
#filter all the discovery SNPs +-500kb of the 313 SNPs or r2>=0.1
#read in nasim_snp_infor
snp.nasim <- read.csv("./data/Nasim_313_SNPs_infor.csv",header=T)
library(MASS)
colnames(snp.nasim)[c(2,3)] <- c("CHR","position")
#load in nasim's oncoarray genotype
load("/data/zhangh24/breast_cancer_data_analysis/risk_prediction/Nasim_prs/result/onco.nasim.snp.rdata")
onco.nasim.genotype <- onco.nasim.snp[,2:ncol(onco.nasim.snp)]
#read in discovery snp infor
load(paste0("./risk_prediction/Nasim_new_prs/result/discover_snp_id.rdata"))
snp.dis <- snp_id
#load dis snp's oncoarray genotype
load("/data/zhangh24/breast_cancer_data_analysis/risk_prediction/Nasim_new_prs/result/onco.dis.snp.rdata")
onco.dis.genotype <- onco.dis.snp[,2:ncol(onco.dis.snp)]
#a bianry vector to indicate whether this snp should be kept in the analysis
idx.cover <- rep(1,32)


for(k in 1:32){
  #chech the position of the SNP
  chr.temp <- snp.dis$CHR[k]
  position.temp <- snp.dis$position[k]
  distance <- 500000
  idx <- which(snp.nasim$CHR==chr.temp&
                 (snp.nasim$position-distance<=position.temp)&
                 ((snp.nasim$position+distance>=position.temp)))
  if(length(idx)!=0){
    idx.cover[k] <- 0
  }
  #check the LD of the snp
  geno.temp <- onco.dis.genotype[,k]
  r2 = cor(geno.temp,onco.nasim.genotype)^2
  idx <- which(r2>=0.1)
  if(length(idx)!=0){
    idx.cover[k] <- 0
  }
}

idx <- which(idx.cover==1)
save(idx,file = "./risk_prediction/Nasim_prs/result/independent_SNPs_idx_in_32_novel_SNP.rdata")

load(paste0("./risk_prediction/Nasim_prs/result/313_intrinsic_subtype_logodds.Rdata"))
nasim_intrinsic_subtypes_result <- final_result
load("./risk_prediction/Nasim_new_prs/result/32_intrinsic_subtype_logodds.Rdata")
dis_intrinsic_subtypes_result <- final_result

library(bc2, lib.loc ="/home/zhangh24/R/x86_64-pc-linux-gnu-library/3.6/")


#merge 313 Nasim SNPs and new selected independent SNPs
intrinsic_subtypes_result <- rbind(nasim_intrinsic_subtypes_result,dis_intrinsic_subtypes_result[idx,])

snpvalue.result <- cbind(onco.nasim.genotype,
                         onco.dis.genotype[,idx])



library(dplyr)
logodds = intrinsic_subtypes_result %>% 
  select(Luminal_A,Luminal_B,Luminal_B_HER2Neg,HER2_Enriched,TN)
all.equal(colnames(snpvalue.result),intrinsic_subtypes_result[,3])
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

quantile.thres <- c(0,0.01,0.05,0.10,
                    0.20,0.40,0.60,0.80,0.90,0.95,0.99,1) 


total <- length(select.names)
odds.ratio <- matrix(0,length(quantile.thres)-1,
                     length(select.names) )
for(j in 1:length(select.names)){

    
    prs <- prs_all[,j]
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
  odds.ratio[k,j] <- odds[k]/odds[6]
}
}
colnames(odds.ratio) <- select.names
write.csv(odds.ratio,file = paste0("./risk_prediction/Nasim_new_prs/result/odds.ratio.csv"))
