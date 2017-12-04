library(R.utils)
library(data.table)
setwd("/spin1/users/zhangh24/breast_cancer_data_analysis/")

n <- 109713
snpvalue <- rep(0,n)
subject.file <- "/gpfs/gsfs4/users/NC_BW/icogs_onco/genotype/imputed2/icogs_order.txt.gz"
Icog.order <- read.table(gzfile(subject.file))

setwd("/spin1/users/zhangh24/breast_cancer_data_analysis/")
data1 <- fread("./data/iCOGS_euro_v10_10232017.csv",header=T)
data1 <- as.data.frame(data1)
y.pheno.mis1 <- cbind(data1$Behaviour1,data1$ER_status1,data1$PR_status1,data1$HER2_status1,data1$Grade1)
colnames(y.pheno.mis1) = c("Behavior","ER","PR","HER2","Grade")
#x.test.all.mis1 <- data1[,c(27:206)]
SG_ID <- data1$SG_ID
x.covar.mis1 <- data1[,c(5:14,204)]
age <- data1[,204]
idx.complete <- which(age!=888)
y.pheno.mis1 <- y.pheno.mis1[idx.complete,]
x.covar.mis1 <- x.covar.mis1[idx.complete,]
SG_ID <- SG_ID[idx.complete]




idx.fil <- Icog.order[,1]%in%SG_ID
idx.match <- match(SG_ID,Icog.order[idx.fil,1])
#Icog.order.match <- Icog.order[idx.fil,1][idx.match]
library(bc2)
#load("./whole_genome_age/ICOG/ERPRHER2_fixed/result/score.test.support.icog.ERPRHER2.Rdata")



Filesdir <- "/gpfs/gsfs4/users/NC_BW/icogs_onco/genotype/imputed2/icogs_imputed/"
Files <- dir(Filesdir,pattern="icogs_merged_b1_12.",full.names=T)
Filesex <- dir(Filesdir,pattern="icogs_merged_b1_12.chr23",full.names=T)
idx.sex <- Files%in%Filesex
Files <- Files[!idx.sex]
library(gtools)
Files <- mixedsort(Files)



qctool.command <- rep("c",564)
qctool.command <- data.frame(qctool.command,stringsAsFactors=F)


for(i in 1:564){
  geno.file <- Files[i]
  temp <- paste0("/spin1/users/zhangh24/qctool_v1.4-linux-x86_64/qctool -g ",Files[i]," -incl-rsids /spin1/users/zhangh24/breast_cancer_data_analysis/discovery_SNP/conditional_analysis/result/all.conditional.snps.icogs.txt -og /spin1/users/zhangh24/breast_cancer_data_analysis/discovery_SNP/conditional_analysis/result/conditional_extract_icog",i,".txt")
  qctool.command[i,1] <- temp
  
}


write.table(qctool.command,file = paste0("/spin1/users/zhangh24/breast_cancer_data_analysis/discovery_SNP/conditional_analysis/code/qc_extract_conditional_icog.sh"),col.names = F,row.names = F,quote=F)















setwd("/spin1/users/zhangh24/breast_cancer_data_analysis/")

subject.file <- "/gpfs/gsfs4/users/NC_BW/icogs_onco/genotype/imputed2/onco_order.txt.gz"
onco.order <- read.table(gzfile(subject.file))

library(data.table)
data2 <- fread("./data/Onco_euro_v10_10232017.csv",header=T)
data2 <- as.data.frame(data2)

#Icog.order.match <- Icog.order[idx.fil,1][idx.match]
library(bc2)
#load("./whole_genome/ONCO/ERPRHER2_fixed/result/score.test.support.onco.ERPRHER2.Rdata")
load("/spin1/users/zhangh24/breast_cancer_data_analysis/whole_genome_age/ICOG/ERPRHER2GRADE_fixed_baseline/result/extract_list.Rdata")

#extract.list <- extract.list[-c(1721,1722,1724,1726,1727,1728),]
#save(extract.list,file="/spin1/users/zhangh24/breast_cancer_data_analysis/whole_genome_age/ICOG/ERPRHER2GRADE_fixed_baseline/result/extract_list.Rdata")


Filesdir <- "/gpfs/gsfs4/users/NC_BW/icogs_onco/genotype/imputed2/onco_imputed"
Files <- dir(Filesdir,pattern="OncoArray_european_merged_b1_15.",full.names=T)
Filesex <- dir(Filesdir,pattern="OncoArray_european_merged_b1_15.chr23",full.names=T)
idx.sex <- Files%in%Filesex
Files <- Files[!idx.sex]
library(gtools)
Files <- mixedsort(Files)


qctool.command <- rep("c",567)
qctool.command <- data.frame(qctool.command,stringsAsFactors=F)


for(i in 1:567){
  geno.file <- Files[i]
  temp <- paste0("/spin1/users/zhangh24/qctool_v1.4-linux-x86_64/qctool -g ",Files[i]," -incl-rsids /spin1/users/zhangh24/breast_cancer_data_analysis/discovery_SNP/conditional_analysis/result/all.conditional.snps.onco.txt -og /spin1/users/zhangh24/breast_cancer_data_analysis/discovery_SNP/conditional_analysis/result/conditional_extract_onco",i,".txt")
  qctool.command[i,1] <- temp
  
}


write.table(qctool.command,file = paste0("/spin1/users/zhangh24/breast_cancer_data_analysis/discovery_SNP/conditional_analysis/code/qc_extract_conditional_onco.sh"),col.names = F,row.names = F,quote=F)






