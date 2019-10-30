rm(list=ls())
#args=(commandArgs(TRUE))
#for(p in 1:length(args)){
#       eval(parse(text=args[[p]]))
#  }
#print(i)
#i1 <- i
#commandarg <- commandArgs(trailingOnly=F)
#myarg <- commandarg[length(commandarg)]
#myarg <- sub("-","",myarg)
#i <- as.numeric(myarg)
#print(i)
#pheno is ICOGS,data2 is onco_array
arg <- commandArgs(trailingOnly=T)
i <- as.numeric(arg[[1]])
i1 <- i
print(i)
library(R.utils)
setwd("/data/zhangh24/breast_cancer_data_analysis/")

n <- 109713
snpvalue <- rep(0,n)
subject.file <- "/gpfs/gsfs4/users/NC_BW/icogs_onco/genotype/imputed2/icogs_order.txt.gz"
Icog.order <- read.table(gzfile(subject.file))
pheno.file <- "./data/pheno.Icog"
load(pheno.file)
n.sub = nrow(pheno)
y.pheno.mis1 <- cbind(pheno$Behaviour1,pheno$PR_status1,pheno$ER_status1,pheno$HER2_status1)
colnames(y.pheno.mis1) = c("Behaviour1","PR_status1",
                           "ER_status1","HER2_status1")
idx.fil <- Icog.order[,1]%in%pheno$SG_ID
idx.match <- match(pheno$SG_ID,Icog.order[idx.fil,1])
#Icog.order.match <- Icog.order[idx.fil,1][idx.match]
library(bc2)
#load("./whole_genome/ICOG/ERPRHER2_fixed/result/score.test.support.icog.ERPRHER2.Rdata")
load("/data/zhangh24/breast_cancer_data_analysis/whole_genome/ICOG/ERPRHER2GRADE_fixed_baseline/result/extract_list.Rdata")

snp.icogs.extract.id <- extract.list[,(ncol(extract.list)-2),drop=F]
write.table(snp.icogs.extract.id,file = paste0("./whole_genome/ICOG/ERPRHER2GRADE_fixed_baseline/result/extract_id_icog.txt"),quote = F,row.names=F)
Filesdir <- "/gpfs/gsfs4/users/NC_BW/icogs_onco/genotype/imputed2/icogs_imputed/"

Files <- dir(Filesdir,pattern="icogs_merged_b1_12.",full.names=T)
Filesex <- dir(Filesdir,pattern="icogs_merged_b1_12.chr23",full.names=T)
idx.sex <- Files%in%Filesex
Files <- Files[!idx.sex]
geno.file <- Files[i1]


qctool.command <- rep("c",564)
qctool.command <- data.frame(qctool.command,stringsAsFactors=F)


for(i in 1:564){
  geno.file <- Files[i]
  temp <- paste0("/data/zhangh24/qctool_v1.4-linux-x86_64/qctool -g ",Files[i]," -incl-rsids /data/zhangh24/breast_cancer_data_analysis/whole_genome/ICOG/ERPRHER2GRADE_fixed_baseline/result/extract_id_icog.txt -og /data/zhangh24/breast_cancer_data_analysis/whole_genome/ICOG/ERPRHER2GRADE_fixed_baseline/result/ERPRHER2_extract",i,".txt")
  qctool.command[i,1] <- temp
  
}


write.table(qctool.command,file = paste0("/data/zhangh24/breast_cancer_data_analysis/whole_genome/ICOG/ERPRHER2GRADE_fixed_baseline/code/qc_extract_snp.sh"),col.names = F,row.names = F,quote=F)
#qctool.command <- as.data.frame(qctool.command,stringsAsFactors=F)





