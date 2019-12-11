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
#load all the datasets
#the subtypes names in sample data
select.names <- c("Luminal_A",
                  "Luminal_B",
                  "Luminal_B_HER2Neg",
                  "HER2Enriched",
                  "TripleNeg")

names.subtypes <- c("stan_logodds")
load(paste0("./risk_prediction/result/split.id.rdata"))
#icog.test.id <- Generatetestid(subtypes.icog)
#load the sample data
icog.train.id <- split.id[[1]]
onco.train.id <- split.id[[2]]
onco.test.id <- split.id[[3]]
#onco.test.id <- split.id[[5]]
icog.vad.id <- split.id[[4]]
onco.vad.id <- split.id[[5]]

#sample.data <- read.table("/data/zhangh24/BCAC/impute_onco/sample.txt",header=T,stringsAsFactors = F)
#sample.data <- sample.data[-1,,drop=F]


onco.data <- as.data.frame(fread("/data/zhangh24/breast_cancer_data_analysis/data/sig_snp_onco_prs.csv",header=T))
onco.data <- onco.data[,-1]

y.pheno.mis2 <- cbind(onco.data$Behavior,
                      onco.data$ER,
                      onco.data$PR,
                      onco.data$HER2,
                      onco.data$Grade)

all.data <- data.frame(onco.data[,1],
                       onco.data$Behavior,
                       stringsAsFactors = F)
idx.fil <- which(all.data[,1]%in%onco.test.id)
idx.match <- match(onco.test.id,all.data[idx.fil,1])
test.data <- all.data[idx.fil[idx.match],]
all.equal(test.data[,1],onco.test.id)




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
    prs <- as.data.frame(fread(paste0("/data/zhangh24/breast_cancer_data_analysis/risk_prediction/overall_prs/result/",select.names[j],"_prs_",i,"_out.profile"),header=T))
    #prs[,4] <- prs.la.temp
    temp <- j%%5
    if(temp==0){temp=5}
    
    #select the sample in the test data
    #also select sample for control and target subtypes
    ID_all <- as.numeric(gsub("sample_","",prs$IID))
    idx.fil <- which(ID_all%in%onco.test.id)
    idx.match <- match(onco.test.id,ID_all[idx.fil])
    ID_test = ID_all[idx.fil[idx.match]]     
    all.equal(ID_test,onco.test.id)    
    score <- prs[idx.fil[idx.match],4]
    
    #temp.log.odds <- read.csv("/data/zhangh24/breast_cancer_data_analysis/risk_prediction/two_stage_model/result/temp_logodds_result.csv")
    
    test.sample <- cbind(test.data,score)
    
    colnames(test.sample) <- c("IID","case","SCORE")
    colnames(test.sample)[1] <- "IID"
    
    
    # idx.fil <- which(prs[,1]%in%test.sample[,1])
    # idx.match <- match(prs[idx.fil,1],test.sample[,1])
    #    test.sample_new2 <- prs[idx.fil[idx.match],] 
    
    
    
    prs.standard <- test.sample[,"SCORE"]
    y.test <- test.sample[,"case"]
    
    
    roc.standard <- roc(y.test,as.vector(prs.standard),ci=T,plot=F)
    # plot.roc(roc.standard)
    #roc.standard
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
    
    code <- paste0("wc -l /data/zhangh24/breast_cancer_data_analysis/risk_prediction/overall_prs/result/",select.names[j],"_prs_pvaluecut_",i,".file")
    temp.out <- system(code,intern = T)
    n.snp[ind] <- as.numeric(gsub(paste0("/data/zhangh24/breast_cancer_data_analysis/risk_prediction/overall_prs/result/",select.names[j],"_prs_pvaluecut_",i,".file"),"",temp.out))
    
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
