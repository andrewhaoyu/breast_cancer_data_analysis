load("/data/zhangh24/breast_cancer_data_analysis/whole_genome_age/ICOG/ERPRHER2GRADE_fixed_baseline/result/meta_result_shared_1p.Rdata")
new_filter <- read.csv("/data/zhangh24/breast_cancer_data_analysis/data/Filter_based_on_CIMBA.csv",header=T,stringsAsFactors = F)
colnames(new_filter) <- c("SNP","position","CHR")
new_filter[,2] <- as.numeric(gsub(",","",new_filter[,2]))



idx <- which((meta_result_shared_1p$position>=(new_filter$position-500000)&
                (meta_result_shared_1p$position<=(new_filter$position+500000)&(meta_result_shared_1p$CHR%in%new_filter$CHR))))
CIMBA_snp <- meta_result_shared_1p[idx,]
save(CIMBA_snp,file = paste0("./discovery_SNP/CIMBA_conditional_analysis/result/CIMBA_snp_name_match.Rdata"))
# save(Julie_snp,file=paste0("./whole_genome_age/ICOG/ERPRHER2GRADE_fixed_baseline/result/Julie_snp_name_match.Rdata"))

setwd("/data/zhangh24/breast_cancer_data_analysis/")
snp.icogs.extract.id <- CIMBA_snp$SNP.ICOGS
write.table(snp.icogs.extract.id,file = paste0("./discovery_SNP/CIMBA_conditional_analysis/result/extract_id_icog_CIMBA.txt"),quote = F,row.names=F)
snp.onco.extract.id <- CIMBA_snp$SNP.ONCO
write.table(snp.onco.extract.id,file = paste0("./discovery_SNP/CIMBA_conditional_analysis/result/extract_id_onco_CIMBA.txt"),quote = F,row.names=F)


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
  temp <- paste0("/data/zhangh24/qctool_v1.4-linux-x86_64/qctool -g ",Files[i]," -incl-rsids /data/zhangh24/breast_cancer_data_analysis/discovery_SNP/CIMBA_conditional_analysis/result/extract_id_icog_CIMBA.txt -og /data/zhangh24/breast_cancer_data_analysis/discovery_SNP/CIMBA_conditional_analysis/result/ICOG/CIMBA_Icog",i,".txt")
  qctool.command[i,1] <- temp
  
}


write.table(qctool.command,file = paste0("/data/zhangh24/breast_cancer_data_analysis/discovery_SNP/CIMBA_conditional_analysis/code/qc_extract_snp_CIMBA_icog.sh"),col.names = F,row.names = F,quote=F)












#################extract the SNPs from ONCO
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
  temp <- paste0("/data/zhangh24/qctool_v1.4-linux-x86_64/qctool -g ",Files[i]," -incl-rsids /data/zhangh24/breast_cancer_data_analysis/discovery_SNP/CIMBA_conditional_analysis/result/extract_id_onco_CIMBA.txt -og /data/zhangh24/breast_cancer_data_analysis/discovery_SNP/CIMBA_conditional_analysis/result/ONCO/CIMBA_onco",i,".txt")
  qctool.command[i,1] <- temp
  
}


write.table(qctool.command,file = paste0("/data/zhangh24/breast_cancer_data_analysis/discovery_SNP/CIMBA_conditional_analysis/code/qc_extract_snp_CIMBA_onco.sh"),col.names = F,row.names = F,quote=F)








