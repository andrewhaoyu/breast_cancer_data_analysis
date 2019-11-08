#-------------------------------------------------------------------
# Update Date: 11/22/2018
# Create Date: 11/22/2018
# Goal: calculate auc
# Author: Haoyu Zhang
#-------------------------------------------------------------------
#---------------------------------------#---------------------------------------
setwd('/data/zhangh24/breast_cancer_data_analysis/')
library(dplyr)
library(data.table)
library(pROC)
#library(bcutility)
#subtypes names in the files, there is a difference, one is called HER2_Enriched, the other one is called HER2Enriched
select.names <- c("Luminal_A",
              "Luminal_B",
              "Luminal_B_HER2Neg",
              "HER2_Enriched",
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
onco.vad.id <- split.id[[4]]

sample.data <- read.table("/data/zhangh24/BCAC/impute_onco/sample.txt",header=T)
sample.data <- sample.data[-1,,drop=F]

onco.data <- as.data.frame(fread("/data/zhangh24/breast_cancer_data_analysis/data/sig_snp_onco_prs.csv",header=T))
onco.data <- onco.data[,-1]




onco.test.id <- data.frame(ID = onco.test.id)
colnames(onco.test.id) = "ID"

new.onco.data <- left_join(onco.test.id,
                           onco.data
                           ,by="ID")


temp <- as.matrix(cbind(
            2-new.onco.data[,c(20,30),drop=F],
            new.onco.data[,c(23),drop=F]))%*%
            as.vector(beta)/6
new.temp <- cbind(new.onco.data[,1],temp,2-new.onco.data[,c(20,30,23)])
head(new.temp)
idx <- which(new.onco.data[,1]==39183 )
new.onco.data[idx,]
n <- length(onco.test.id)
test_ID <- matrix("c",n,1)
for(i in 1:n){
  #fix the issue that prs file will code 100000 as 1e+05
  if(onco.test.id[i]==100000){
    test_ID[i] <- paste0("sample_1e+05")  
  }else{
    test_ID[i] <- paste0("sample_",as.numeric(onco.test.id[i]))
  }
  
}
#test_ID
#select.names <- c(subtypes,paste0("eb_",subtypes))
#select.names <- c(rep("standard",length(subtypes)),
                 # select.names)
pthres <- c(5E-08,1E-07,5E-07,1E-06,5E-06,1E-05,5E-05,
            1E-04,5E-04,1E-03,5E-03,1E-02)
#create the file for standard logisic regression
n.pthres <- length(pthres)

total <- length(select.names)*n.pthres
auc <- rep(0,total)
auc95 <- rep("c",total)
method <- rep("c",total)
p <- rep(0,total)
subtypes <- rep("c",total)
n.snp <- rep(0,total)
ind <- 1
for(j in 1:length(select.names)){
  for(i in 1:n.pthres){
    #read in result
    #load the PRS information for the subtypes
    #the prs files contain all the subjects in the genotyped data
    #the genotype data is larger than the phenotype data
    #we need to select the subset for testdata
    prs <- as.data.frame(fread(paste0("/data/zhangh24/breast_cancer_data_analysis/risk_prediction/subtypes_prs/result/",select.names[j],"_prs_",i,"_out.profile"),header=T))
   
      temp <- j%%5
    if(temp==0){temp=5}
    
      #select the sample in the test data
      #also select sample for control and target subtypes
    
      
      test.sample = sample.data%>%
      filter(
        (ID_1%in%test_ID))
        
        
        
        &
          ((subtypes=="control"|
              subtypes==names.subtypes[temp])))%>%
      select(ID_1,case)
    
    colnames(test.sample) <- c("IID","case")
    colnames(test.sample)[1] <- "IID"
    #match the test sample
    #get the prs from the prs file of all the genotype data
    test.prs <- read.table("/data/zhangh24/breast_cancer_data_analysis/risk_prediction/subtypes_prs/result/test_out.profile",header=T)
    prs <- test.prs
    
    test.sample_new <- left_join(test.sample,prs)
    idx <- which(test.sample_new[,1]==test_ID[2])
    test.sample_new[idx,]
    idx <- which(test.sample_new$case==1)
    idx2<- which(test.sample_new$case==0)
    mean(test.sample_new$SCORE[idx])
    mean(test.sample_new$SCORE[idx2])
    head(test.sample_new[idx,])    
    head(test.sample_new[idx2,])
    tail(test.sample_new[idx,])    
    tail(test.sample_new[idx2,])
    prs.standard <- test.sample_new[,"SCORE"]
    y.test <- test.sample_new[,"case"]
    #roc.standard <- calibration(y.test,prs.standard)
    
    
    
    
    roc.standard <- roc(y.test,as.vector(prs.standard),ci=T,plot=T)
    plot.roc(roc.standard)
    roc.standard
    pla <- 4
    auc[ind] <- round(as.numeric(roc.standard$auc),pla)*100
    
    auc95[ind] <- paste0(round(as.numeric(roc.standard$auc)[1],pla)*100,
                  "(",
      round(as.numeric(roc.standard$ci)[1],pla)*100,
                        "-",
                        round(as.numeric(roc.standard$ci)[3],pla)*100
      ,")")
   # if(j<=5){
    #  method.temp <- "standard"  
    #}else if(j<=10){
     
     method.temp <- "two-stage"
    #}else{method.temp = "eb"}
    subtypes[ind] <- names.subtypes[temp]
    method[ind] <- method.temp
    p[ind] <- pthres[i]
    #read in SNP file information to calculate PRS
    
    code <- paste0("wc -l /data/zhangh24/breast_cancer_data_analysis/risk_prediction/subtypes_prs/result/",select.names[j],"_prs_pvaluecut_",i,".file")
    temp.out <- system(code,intern = T)
    n.snp[ind] <- as.numeric(gsub(paste0("/data/zhangh24/breast_cancer_data_analysis/risk_prediction/subtypes_prs/result/",select.names[j],"_prs_pvaluecut_",i,".file"),"",temp.out))
    
    ind = ind + 1
    
#    cal.standard <- calibration(y.test,prs.standard)
    
    
    
    
    
  }
}

auc.result <- data.frame(auc,
                         auc95,
                         method,
                         p,
                         subtypes,
                         n.snp)
save(auc.result,file = "/data/zhangh24/breast_cancer_data_analysis/risk_prediction/LD_clumping/result/auc.result.rdata")
