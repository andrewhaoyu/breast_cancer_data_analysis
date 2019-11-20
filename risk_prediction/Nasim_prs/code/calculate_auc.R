#calculate auc based on nasim 313 SNPs
setwd("/data/zhangh24/breast_cancer_data_analysis/")
load(paste0("./risk_prediction/Nasim_prs/result/nasim_snp_odds.rdata"))


sample.data <- read.table("/data/zhangh24/BCAC/impute_onco/sample.txt",header=T,stringsAsFactors = F)
sample.data <- sample.data[-1,,drop=F]
n <- nrow(sample.data)  
extract.num <- nrow(snp_odds)
snpid.result <- rep("c",extract.num)
n.sub <- nrow(sample.data)
snpvalue.result <- matrix(0,n.sub,extract.num)


total <- 0

for(i in 1:567){
  
  print(i)  
  geno.file <- paste0("/data/zhangh24/breast_cancer_data_analysis/risk_prediction/Nasim_prs/result/extracted_snp_313_",i,".txt")
  num <- as.numeric(system(paste0('cat ',geno.file,' | wc -l '),intern=T))
  
  if(num!=0){
    
    
    con <- file(geno.file)
    temp <- 0
    open(con)
    for(i in 1:num){
      
      oneLine <- readLines(con,n=1)
      myVector <- strsplit(oneLine," ")
      snpid <- as.character(myVector[[1]][3])
      
      
      temp <- temp+1
      snpid.result[temp+total] <- snpid
      snpvalue <- rep(0,n)
      
      
      snppro <- as.numeric(unlist(myVector)[7:length(myVector[[1]])])
      
      snpvalue <- convert(snppro,n.raw)
      #snpvalue <- snpvalue[idx.fil][idx.match]
      snpvalue.result[,temp+total] <- snpvalue
      
      
    }
    close(con)
    
    total <- total+num
    
    
  }
  
  
}

snpid.result <- snpid.result[1:total]
snpvalue.result <- snpvalue.result[,1:total]
#three snps were not included due to close to 0.01 MAF
#need to be fixed
snp_odds <- snp_odds[complete.cases(snp_odds),]

extract.list.shared <- snp_odds[,4]
idx.match <- match(extract.list.shared,snpid.result)
snpid.result <- snpid.result[idx.match]
all.equal(snpid.result,extract.list.shared)
snpvalue.result <- snpvalue.result[,idx.match]

logodds = snp_odds %>% 
  select(Overall.Breast.Cancerd,ER.positivee,ER.negativef,Luminal_A,Luminal_B,Luminal_B_HER2Neg,HER2_Enriched,TN)

prs_all <- snpvalue.result%*%as.matrix(logodds)

#use overall analysis to predict auc results
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
#onco.test.id <- split.id[[5]]
icog.vad.id <- split.id[[4]]
onco.vad.id <- split.id[[4]]

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
write.csv(auc.result,file = paste0("./risk_prediction/Nasim_prs/result/nasim_snp_auc_result.csv"))
