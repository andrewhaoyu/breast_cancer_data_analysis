##############generate the discovery snps list
load("/spin1/users/zhangh24/breast_cancer_data_analysis/whole_genome_age/ICOG/ERPRHER2GRADE_casecase/result/extract_result_shared.Rdata")
# load("/spin1/users/zhangh24/breast_cancer_data_analysis/whole_genome_age/ICOG/ERPRHER2GRADE_fixed_baseline/result/LD_pruning.result")
# load("/spin1/users/zhangh24/breast_cancer_data_analysis/whole_genome_age/ICOG/ERPRHER2GRADE_fixed_baseline/result/extract_list.Rdata")
load("/spin1/users/zhangh24/breast_cancer_data_analysis/whole_genome_age/ICOG/ERPRHER2GRADE_fixed_baseline/result/meta_result_shared_1p.Rdata")
discovery_snp <- read.csv("/spin1/users/zhangh24/breast_cancer_data_analysis/data/SNP_infor_conditional_check.csv",header=T)
discovery_snp <- discovery_snp[1:22,]

###############find the potential SNPs 
chr.pos.dis <- paste0(discovery_snp$CHR,":",
                       discovery_snp$Pos)
 
 
 
 chr.pos.all <- paste0(meta_result_shared_1p$CHR,":",
                       meta_result_shared_1p$position)
 idx.fil <- which(chr.pos.all%in%chr.pos.dis)
 length(chr.pos.all[idx.fil])
 discovery_snp<-  meta_result_shared_1p[idx.fil,]
# 



 setwd("/spin1/users/zhangh24/breast_cancer_data_analysis/")
 snp.icogs.extract.id <- as.character(discovery_snp$SNP.ICOGS)
 write.table(snp.icogs.extract.id,file = paste0("./discovery_SNP/result/extract_id_icog_discovery.txt"),quote = F,row.names=F)
 snp.onco.extract.id <- as.character(discovery_snp$SNP.ONCO)
 write.table(snp.onco.extract.id,file = paste0("./discovery_SNP/result/extract_id_onco_discovery.txt"),quote = F,row.names=F)
 
 
 # /spin1/users/zhangh24/breast_cancer_data_analysis/discovery_SNP/result/extract_id_icog_discovery.txt
 
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
   temp <- paste0("/spin1/users/zhangh24/qctool_v1.4-linux-x86_64/qctool -g ",Files[i]," -incl-rsids /spin1/users/zhangh24/breast_cancer_data_analysis/discovery_SNP/result/extract_id_icog_discovery.txt -og /spin1/users/zhangh24/breast_cancer_data_analysis/discovery_SNP/result/discovery_icog",i,".txt")
   qctool.command[i,1] <- temp
   
 }
 
 
 write.table(qctool.command,file = paste0("/spin1/users/zhangh24/breast_cancer_data_analysis/discovery_SNP/code/qc_extract_discovery_icog.sh"),col.names = F,row.names = F,quote=F)
 
 
 
 
 
 Filesdir <- "/gpfs/gsfs4/users/NC_BW/icogs_onco/genotype/imputed2/onco_imputed"
 Files <- dir(Filesdir,pattern="OncoArray_european_merged_b1_15.",full.names=T)
 Filesex <- dir(Filesdir,pattern="OncoArray_european_merged_b1_15.chr23",full.names=T)
 idx.sex <- Files%in%Filesex
 Files <- Files[!idx.sex]
 Files <- mixedsort(Files)
 
 
 qctool.command <- rep("c",567)
 qctool.command <- data.frame(qctool.command,stringsAsFactors=F)
 
 
 for(i in 1:567){
   geno.file <- Files[i]
   temp <- paste0("/spin1/users/zhangh24/qctool_v1.4-linux-x86_64/qctool -g ",Files[i]," -incl-rsids /spin1/users/zhangh24/breast_cancer_data_analysis/discovery_SNP/result/extract_id_onco_discovery.txt -og /spin1/users/zhangh24/breast_cancer_data_analysis/discovery_SNP/result/discovery_onco",i,".txt")
   qctool.command[i,1] <- temp
   
 }
 
 
 write.table(qctool.command,file = paste0("/spin1/users/zhangh24/breast_cancer_data_analysis/discovery_SNP/code/qc_extract_discovery_onco.sh"),col.names = F,row.names = F,quote=F)
 
 
 
 
 
 
 







