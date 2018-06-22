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

extract_list <- 

#extract.list <- extract.list[-c(1721,1722,1724,1726,1727,1728),]
#save(extract.list,file="/spin1/users/zhangh24/breast_cancer_data_analysis/whole_genome_age/ICOG/ERPRHER2GRADE_fixed_baseline/result/extract_list.Rdata")

snp.onco.extract.id <- extract.list[,(ncol(extract.list)-1),drop=F]
write.table(snp.onco.extract.id,file = paste0("./whole_genome_age/ONCO/ERPRHER2GRADE_fixed_baseline/result/extract_id_onco.txt"),quote = F,row.names=F)



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
  temp <- paste0("/spin1/users/zhangh24/qctool_v1.4-linux-x86_64/qctool -g ",Files[i]," -incl-rsids /spin1/users/zhangh24/breast_cancer_data_analysis/whole_genome_age/ONCO/ERPRHER2GRADE_fixed_baseline/result/extract_id_onco.txt -og /spin1/users/zhangh24/breast_cancer_data_analysis/whole_genome_age/ONCO/ERPRHER2GRADE_fixed_baseline/result/ERPRHER2_extract",i,".txt")
  qctool.command[i,1] <- temp
  
}


write.table(qctool.command,file = paste0("/spin1/users/zhangh24/breast_cancer_data_analysis/whole_genome_age/ONCO/ERPRHER2GRADE_fixed_baseline/code/qc_extract_snp.sh"),col.names = F,row.names = F,quote=F)





