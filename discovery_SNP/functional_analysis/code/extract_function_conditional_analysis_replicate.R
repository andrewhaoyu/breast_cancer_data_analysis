load("/data/zhangh24/breast_cancer_data_analysis/whole_genome_age/ICOG/ERPRHER2GRADE_fixed_baseline/result/meta_result_shared_1p.Rdata")
library(data.table)
new_filter <- fread("/data/zhangh24/breast_cancer_data_analysis/data/two_new_SNPs_for_replicate.csv",header=T,stringsAsFactors = F)
colnames(new_filter) <- c("SNP","position","CHR")
#new_filter[,2] <- as.numeric(gsub(",","",new_filter[,2]))
new_filter <- as.data.frame(new_filter)
idx <- NULL
for(i in 1:nrow(new_filter)){
  idx.temp <- which((meta_result_shared_1p$position>=(new_filter$position[i]-500000)&
                       (meta_result_shared_1p$position<=(new_filter$position[i]+500000)&(meta_result_shared_1p$CHR==new_filter$CHR[i]))))
  print(length(idx.temp))
  idx <- c(idx,idx.temp)
  
}



setwd('/data/zhangh24/breast_cancer_data_analysis/')
functional_snp_conditional <- meta_result_shared_1p[idx,]
save(functional_snp_conditional,file = paste0("./discovery_SNP/functional_analysis/result/functional_snp_conditional_replicate.Rdata"))
# save(Julie_snp,file=paste0("./whole_genome_age/ICOG/ERPRHER2GRADE_fixed_baseline/result/Julie_snp_name_match.Rdata"))

setwd("/data/zhangh24/breast_cancer_data_analysis/")
snp.icogs.extract.id <- functional_snp_conditional$SNP.ICOGS
write.table(snp.icogs.extract.id,file = paste0("./discovery_SNP/functional_analysis/result/extract_id_icog_functional_conditional_replicate.txt"),quote = F,row.names=F)
snp.onco.extract.id <- functional_snp_conditional$SNP.ONCO
write.table(snp.onco.extract.id,file = paste0("./discovery_SNP/functional_analysis/result/extract_id_onco_functional_conditional_replicate.txt"),quote = F,row.names=F)


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
  temp <- paste0("/data/zhangh24/qctool_v1.4-linux-x86_64/qctool -g ",Files[i]," -incl-rsids /data/zhangh24/breast_cancer_data_analysis/discovery_SNP/functional_analysis/result/extract_id_icog_functional_conditional_replicate.txt -og /data/zhangh24/breast_cancer_data_analysis/discovery_SNP/functional_analysis/result/ICOG/funcitonal_conditional_icog_replicate",i,".txt")
  qctool.command[i,1] <- temp
  
}


write.table(qctool.command,file = paste0("/data/zhangh24/breast_cancer_data_analysis/discovery_SNP/functional_analysis/code/qc_extract_snp_funcitonal_conditional_icog_replicate.sh"),col.names = F,row.names = F,quote=F)












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
  temp <- paste0("/data/zhangh24/qctool_v1.4-linux-x86_64/qctool -g ",Files[i]," -incl-rsids /data/zhangh24/breast_cancer_data_analysis/discovery_SNP/functional_analysis/result/extract_id_onco_functional_conditional_replicate.txt -og /data/zhangh24/breast_cancer_data_analysis/discovery_SNP/functional_analysis/result/ONCO/functional_conditional_onco_replicate",i,".txt")
  qctool.command[i,1] <- temp
  
}


write.table(qctool.command,file = paste0("/data/zhangh24/breast_cancer_data_analysis/discovery_SNP/functional_analysis/code/qc_extract_snp_functional_conditional_onco_replicate.sh"),col.names = F,row.names = F,quote=F)










