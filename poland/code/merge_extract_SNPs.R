
library(R.utils)
library(data.table)
setwd("/data/zhangh24/breast_cancer_data_analysis/")

subject.file <- "/gpfs/gsfs4/users/NC_BW/icogs_onco/genotype/imputed2/onco_order.txt.gz"
onco.order <- read.table(gzfile(subject.file))
setwd("/data/zhangh24/breast_cancer_data_analysis/")
data2 <- fread("./data/Onco_euro_v10_10232017.csv",header=T)
data2 <- as.data.frame(data2)
#idx <- which(data2$StudyCountry=="Poland")
#data2 <- data2[idx,]

y.pheno.mis2 <- cbind(data2$Behaviour1,data2$ER_status1,data2$PR_status1,data2$HER2_status1,data2$Grade1)
#y.pheno.mis2 <- cbind(data2$Behaviour1,data2$PR_status1,data2$ER_status1,data2$HER2_status1)
colnames(y.pheno.mis2) = c("Behaviour","ER",
                           "PR","HER2","Grade")

x.covar.mis2 <- data2[,c(5:8,204)]
Onc_ID <- data2$Onc_ID



idx.fil <- onco.order[,1]%in%Onc_ID
n <- length(idx.fil)
snpvalue <- rep(0,n)
idx.match <- match(Onc_ID,onco.order[idx.fil,1])
#Icog.order.match <- Icog.order[idx.fil,1][idx.match]
library(bc2)
#load("./whole_genome/ONCO/ERPRHER2GRADE_fixed_baseline/result/score.test.support.onco.ERPRHER2Grade.Rdata")

top_signal_in_poland <- matrix(0,nrow(data2),2)
library(data.table)
result <- NULL
length.file <- rep(0,567)
file.id <- c(342,384)
for(i in 1:2){
  
  file <- paste0("/data/zhangh24/breast_cancer_data_analysis/poland/result/whole_genome/top_signal_extract",file.id[i],".txt")
  oneLine <- readLines(file)
  myVector <- strsplit(oneLine," ")
  snpid <- as.character(myVector[[1]][3])
  #snpid_result[temp] <- snpid
  snpvalue <- rep(0,n)
  
  
  snppro <- as.numeric(unlist(myVector)[7:length(myVector[[1]])])
  snpvalue <- convert(snppro,n)
  snpvalue <- snpvalue[idx.fil][idx.match]
  top_signal_in_poland[,i] <- snpvalue
}

colnames(top_signal_in_poland) <- c("rs11200014","rs78540526")
write.csv(top_signal_in_poland,file="/data/zhangh24/breast_cancer_data_analysis/poland/result/whole_genome/top_signal_in_poland.csv")  
