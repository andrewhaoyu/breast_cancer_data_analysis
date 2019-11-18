#-------------------------------------------------------------------
# Update Date: 11/15/2018
# Create Date: 11/13/2018
# Goal: prepare sample dataset for all onco
# Goal: prepare sample dataset for subtyeps test 
# Author: Haoyu Zhang
#-------------------------------------------------------------------
#prepare sample for all of the subjects
library(withr)
#with_libpaths(new = "/home/zhangh24/R/x86_64-pc-linux-gnu-library/3.6", install_github("andrewhaoyu/bcutility"))
#install_github("andrewhaoyu/bcutility")
library(bcutility,lib.loc = "/spin1/home/linux/zhangh24/R/x86_64-pc-linux-gnu-library/3.6/")
#read in the onco array genotype data file order
subject.file <- "/data/zhangh24/test/onco_order.txt"
library(data.table)
onco.data <- as.data.frame(fread("/data/zhangh24/breast_cancer_data_analysis/data/sig_snp_onco_prs.csv",header=T))
onco.data <- onco.data[,-1]

idx <- which(onco.data[,1]==39177)
library(tidyverse)
y.pheno.mis2 <- select(onco.data,Behavior,ER,PR,HER2,Grade)
############put onco array unknown cases as 1
idx.insi.unknown.onco <- which(y.pheno.mis2[,1]==888|
                                 y.pheno.mis2[,1]==2) 
y.pheno.mis2[idx.insi.unknown.onco,1] <- 1
subtypes.onco <- GenerateIntrinsicmis(y.pheno.mis2[,2],
                                      y.pheno.mis2[,3],
                                      y.pheno.mis2[,4],
                                      y.pheno.mis2[,5])
subtypes.onco <- GenerateIntrinsicmis(y.pheno.mis2[,2],y.pheno.mis2[,3],
                                      y.pheno.mis2[,4],y.pheno.mis2[,5])
disease_onco = select(onco.data,ID,Behavior) %>% 
              cbind(subtypes.onco)
onco.order <- read.table(subject.file)
colnames(onco.order) <- "ID"

onco.order.new <- left_join(onco.order,disease_onco,
                          by= "ID")

#head(disease_onco)

idx <- which(disease_onco[,1]%in%onco.order$ID==T)


#idx <- which(onco.order.new$ID==999999357290)
all.equal(onco.order.new[,1],onco.order[,1])
subtypes.onco <- data.frame(
                    as.character(onco.order.new[,3]),
                            stringsAsFactors = F)

# 
# onco.order.new <- merge(onco.order,disease_onco,
#                         by.x= "ID",
#                         by.y = "ID",
#                         all.x = T
#                         )
# idx.match <- match(as.numeric(onco.order[,1]),
#                    as.numeric(onco.order.new[,1]))
# onco.order.new <- onco.order.new[idx.match,]
# 

n <- nrow(onco.order)
ID <- matrix("c",n,1)
for(i in 1:n){
  #avoid the isseu of coding 100000 as 1e+05 since it will cause trouble in future sample match as character

    ID[i] <- paste0("sample_",as.numeric(onco.order[i,1]))
  
  
}
missing <- matrix(0,n,1)
case <- onco.order.new[,2,drop=F]
#fam file have 6 different columns
#first column is family ID
#second colunn is individual ID
# 3rd column Paternal ID
# 4th column Maternal ID
#5th column (1=male; 2=female; other=unknown)
#sixth column is the case with 2 as cases, 1 as controls, -9 as missing
plink.case <- case+1
plink.case[is.na(case)]= -9
fam <- data.frame(ID,ID,
                  rep(0,n),
                  rep(0,n),
                  rep(2,n),
                  case+1)
write.table(fam,file = "/data/zhangh24/BCAC/impute_onco/onco_plink.fam",row.names = F,col.names = F,quote=F)



cov_1 <- matrix(rnorm(n),n,1)




#onco.order <- matrix(paste0("sample",onco.order),n,1)
ID <- rbind(0,ID)
missing <- rbind(0,missing)
case <- rbind("B",case)
cov_1 <- rbind("C",cov_1)
subtypes.onco <- rbind("D",
                       subtypes.onco)
#onco.order <- matrix(onco.order,ncol=1)

sample.data <- data.frame(ID,
                          ID,
                          missing,
                          case,
                          cov_1,
                          subtypes.onco,
                          stringsAsFactors = F)
colnames(sample.data) <- c("ID_1",
                           "ID_2",
                           "missing",
                           "case",
                           "cov_1",
                           "subtypes")
write.table(sample.data,file = "/data/zhangh24/BCAC/impute_onco/sample.txt",
            row.names = F, quote = F,sep = " ")












#prepare sample for people keeping in test
setwd('/data/zhangh24/breast_cancer_data_analysis/')
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
for(i in 1:length(names.subtypes)){
  test.sample = sample.data%>%
                filter(
                  (ID_1%in%test_ID)&
              ((subtypes=="control"|
                    subtypes==names.subtypes[i]))
                ) %>%
                select(ID_1)
  write.table(test.sample,file = 
                paste0("/data/zhangh24/BCAC/test_sample_",
                        names.subtypes[i]
                        ,".txt"),
              row.names = F, quote = F,sep = " ")
}






write.table(test_ID,file = "/data/zhangh24/BCAC/sample_onco_test.txt",
            row.names = F, quote = F,sep = " ")
