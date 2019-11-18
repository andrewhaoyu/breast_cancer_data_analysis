#goal: extract nasim snps from oncoarray data
setwd("/data/zhangh24/breast_cancer_data_analysis/")
load(paste0("./risk_prediction/Nasim_prs/result/nasim_snp_odds.rdata"))

extract_snp <-  snp_odds %>%
 select(SNP.ONCO)

write.table(extract_snp,file = paste0("./risk_prediction/Nasim_prs/result/extract_snp_313.txt"),row.names=F,col.names=T,quote=F)



Filesdir <- "/gpfs/gsfs4/users/NC_BW/icogs_onco/genotype/imputed2/onco_imputed"
Files <- dir(Filesdir,pattern="OncoArray_european_merged_b1_15.",full.names=T)
Filesex <- dir(Filesdir,pattern="OncoArray_european_merged_b1_15.chr23",full.names=T)
idx.sex <- Files%in%Filesex
Files <- Files[!idx.sex]


qctool.command <- rep("c",567)
qctool.command <- data.frame(qctool.command,stringsAsFactors=F)


for(i in 1:567){
  geno.file <- Files[i]
  
  temp <- paste0("/data/zhangh24/qctool_v1.4-linux-x86_64/qctool -g ",Files[i]," -incl-rsids /data/zhangh24/breast_cancer_data_analysis/risk_prediction/Nasim_prs/result/extract_snp_313.txt -og /data/zhangh24/breast_cancer_data_analysis/risk_prediction/Nasim_prs/result/extracted_snp_313_",i,".txt")
  qctool.command[i,1] <- temp
  
}


write.table(qctool.command,file = paste0("./risk_prediction/Nasim_prs/code/qc_extract_snp.sh"),col.names = F,row.names = F,quote=F)

