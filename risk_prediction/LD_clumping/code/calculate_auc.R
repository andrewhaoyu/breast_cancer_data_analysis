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
library(bcutility)
subtypes <- c("Luminial_A",
              "Luminal_B",
              "Luminal_B_HER2Neg",
              "HER2_Enriched",
              "TN")
#load all the datasets
sample.data <- as.data.frame(fread("/data/zhangh24/test/sample.txt"))
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
n <- length(onco.test.id)
test_ID <- matrix("c",n,1)
for(i in 1:n){
  test_ID[i] <- paste0("sample_",onco.test.id[i])
}
select.names <- c(subtypes,paste0("eb_",subtypes))
select.names <- c(rep("standard",length(subtypes)),
                  select.names)
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
    
    prs <- as.data.frame(fread(paste0("/data/zhangh24/BCAC/prs_out/",select.names[j],"_prs_",i,"_out.profile"),header=T))
    
      temp <- j%%5
    if(temp==0){temp=5}
    
    test.sample = sample.data%>%
      filter(
        (ID_1%in%test_ID)&
          ((subtypes=="control"|
              subtypes==names.subtypes[temp]))
      ) %>%
      select(ID_1,case)
    
    colnames(test.sample) <- c("IID","case")
   #match the test sample
    test.sample <- inner_join(test.sample,prs)
    prs.standard <- test.sample[,"SCORE"]
    y.test <- test.sample[,"case"]
    roc.standard <- calibration(y.test,prs.standard)
    roc.standard <- roc(y.test,as.vector(prs.standard),ci=T,plot=F)
    roc.standard
    pla <- 4
    auc[ind] <- round(as.numeric(roc.standard$auc),pla)*100
    
    auc95[ind] <- paste0(round(as.numeric(roc.standard$auc)[1],pla)*100,
                  "(",
      round(as.numeric(roc.standard$ci)[1],pla)*100,
                        "-",
                        round(as.numeric(roc.standard$ci)[3],pla)*100
      ,")")
    if(j<=5){
      method.temp <- "standard"  
    }else if(j<=10){
      method.temp <- "two-stage"
    }else{method.temp = "eb"}
    subtypes[ind] <- names.subtypes[temp]
    method[ind] <- method.temp
    p[ind] <- pthres[i]
    code <- paste0("wc -l /data/zhangh24/BCAC/prs_file/",select.names[j],"_prs_",i,".file")
    temp.out <- system(code,intern = T)
    n.snp[ind] <- as.numeric(gsub(paste0("/data/zhangh24/BCAC/prs_file/",select.names[j],"_prs_",i,".file"),"",temp.out))
    
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
