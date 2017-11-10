
  load("/spin1/users/zhangh24/breast_cancer_data_analysis/whole_genome_age/ICOG/ERPRHER2GRADE_fixed_baseline/result/meta_result_shared_1p.Rdata")
  new_filter <- read.csv("/spin1/users/zhangh24/breast_cancer_data_analysis/whole_genome/ICOG/ERPRHER2_fixed/result/Filter_based_on_Montse.csv",header=T,stringsAsFactors = F)
  new_filter[,2] <- as.numeric(gsub(",","",new_filter[,2]))
  


idx <- which((meta_result_shared_1p$position%in%new_filter$position)&
(meta_result_shared_1p$CHR%in%new_filter$CHR))

Julie_snp <- meta_result_shared_1p[idx,]


setwd("/spin1/users/zhangh24/breast_cancer_data_analysis/")
snp.icogs.extract.id <- Julie_snp$SNP.ICOGS
write.table(snp.icogs.extract.id,file = paste0("./whole_genome_age/ICOG/ERPRHER2GRADE_fixed_baseline/result/extract_id_icog_Julie.txt"),quote = F,row.names=F)
snp.onco.extract.id <- Julie_snp$SNP.ONCO
write.table(snp.icogs.extract.id,file = paste0("./whole_genome_age/ONCO/ERPRHER2GRADE_fixed_baseline/result/extract_id_onco_Julie.txt"),quote = F,row.names=F)




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
  temp <- paste0("/spin1/users/zhangh24/qctool_v1.4-linux-x86_64/qctool -g ",Files[i]," -incl-rsids /spin1/users/zhangh24/breast_cancer_data_analysis/whole_genome_age/ICOG/ERPRHER2GRADE_fixed_baseline/result/extract_id_icog.txt -og /spin1/users/zhangh24/breast_cancer_data_analysis/whole_genome_age/ICOG/ERPRHER2GRADE_fixed_baseline/result/ERPRHER2_extract",i,".txt")
  qctool.command[i,1] <- temp
  
}


write.table(qctool.command,file = paste0("/spin1/users/zhangh24/breast_cancer_data_analysis/whole_genome_age/ICOG/ERPRHER2GRADE_fixed_baseline/code/qc_extract_snp_julie.sh"),col.names = F,row.names = F,quote=F)

