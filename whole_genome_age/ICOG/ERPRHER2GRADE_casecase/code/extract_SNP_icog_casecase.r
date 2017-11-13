
library(R.utils)
library(data.table)
setwd("/spin1/users/zhangh24/breast_cancer_data_analysis/")

n <- 109713
snpvalue <- rep(0,n)
subject.file <- "/gpfs/gsfs4/users/NC_BW/icogs_onco/genotype/imputed2/icogs_order.txt.gz"
Icog.order <- read.table(gzfile(subject.file))

setwd("/spin1/users/zhangh24/breast_cancer_data_analysis/")
#Icog.order.match <- Icog.order[idx.fil,1][idx.match]
library(bc2)
#load("./whole_genome_age/ICOG/ERPRHER2_fixed/result/score.test.support.icog.ERPRHER2.Rdata")
load("/spin1/users/zhangh24/breast_cancer_data_analysis/whole_genome_age/ICOG/ERPRHER2GRADE_casecase/result/extract_list.Rdata")


snp.icogs.extract.id <- extract.list[,(ncol(extract.list)-2),drop=F]
write.table(snp.icogs.extract.id,file = paste0("./whole_genome_age/ICOG/ERPRHER2GRADE_casecase/result/extract_id_icog.txt"),quote = F,row.names=F)

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
  temp <- paste0("/spin1/users/zhangh24/qctool_v1.4-linux-x86_64/qctool -g ",Files[i]," -incl-rsids /spin1/users/zhangh24/breast_cancer_data_analysis/whole_genome_age/ICOG/ERPRHER2GRADE_casecase/result/extract_id_icog.txt -og /spin1/users/zhangh24/breast_cancer_data_analysis/whole_genome_age/ICOG/ERPRHER2GRADE_casecase/result/ERPRHER2_extract",i,".txt")
  qctool.command[i,1] <- temp
  
}


write.table(qctool.command,file = paste0("/spin1/users/zhangh24/breast_cancer_data_analysis/whole_genome_age/ICOG/ERPRHER2GRADE_casecase/code/qc_extract_snp.sh"),col.names = F,row.names = F,quote=F)


