
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
load("./genetic_correlation/result/hapmap3list.Rdata")


snp.icogs.extract.id <- shared.data[,7,drop=F]
write.table(snp.icogs.extract.id,file = paste0("./genetic_correlation/result/extract_id_icog.txt"),quote = F,row.names=F)

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
  temp <- paste0("/spin1/users/zhangh24/qctool_v1.4-linux-x86_64/qctool -g ",Files[i]," -incl-rsids /spin1/users/zhangh24/breast_cancer_data_analysis/genetic_correlation/result/extract_id_icog.txt -og /spin1/users/zhangh24/breast_cancer_data_analysis/genetic_correlation/ICOG/result/hapmap_icog",i,".txt")
  qctool.command[i,1] <- temp
  
}


write.table(qctool.command,file = paste0("./genetic_correlation/ICOG/code/qc_extract_snp.sh"),col.names = F,row.names = F,quote=F)


