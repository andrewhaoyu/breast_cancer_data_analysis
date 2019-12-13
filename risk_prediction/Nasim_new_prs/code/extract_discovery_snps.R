#read in discovery snps id information
setwd("/data/zhangh24/breast_cancer_data_analysis/")
snp <- read.csv("./data/discovery_snp_paper_order.csv")
library(dplyr)
snp = snp %>% mutate(chr.pos = paste0(CHR,":",position))
#load whole genome information to get the SNP id in OncoArray and iCOGs
load("./risk_prediction/intrinsic_subtypes_whole_genome/ICOG/result/whole_gonome.rdata")
whole_genome = whole_genome %>% 
  mutate(chr.pos = paste0(CHR,":",position))
snp_id = left_join(snp,whole_genome,
                   by="chr.pos") %>% 
  select(var_name,SNP.ICOGS,SNP.ONCO)

#goal: extract discovery snps from oncoarray data
setwd("/data/zhangh24/breast_cancer_data_analysis/")
load(paste0("./risk_prediction/Nasim_prs/result/nasim_snp_odds.rdata"))
load(paste0("./risk_prediction/Nasim_prs/result/nasim_snp_id.rdata"))
#extract oncoarray data
extract_snp <-  snp_id %>%
  select(SNP.ONCO)

write.table(extract_snp,file = paste0("./risk_prediction/Nasim_prs/result/extract_snp_313_onco.txt"),row.names=F,col.names=T,quote=F)



Filesdir <- "/gpfs/gsfs4/users/NC_BW/icogs_onco/genotype/imputed2/onco_imputed"
Files <- dir(Filesdir,pattern="OncoArray_european_merged_b1_15.",full.names=T)
Filesex <- dir(Filesdir,pattern="OncoArray_european_merged_b1_15.chr23",full.names=T)
idx.sex <- Files%in%Filesex
Files <- Files[!idx.sex]


qctool.command <- rep("c",567)
qctool.command <- data.frame(qctool.command,stringsAsFactors=F)


for(i in 1:567){
  geno.file <- Files[i]
  
  temp <- paste0("/data/zhangh24/qctool_v1.4-linux-x86_64/qctool -g ",Files[i]," -incl-rsids /data/zhangh24/breast_cancer_data_analysis/risk_prediction/Nasim_prs/result/extract_snp_313_onco.txt -og /data/zhangh24/breast_cancer_data_analysis/risk_prediction/Nasim_prs/result/extracted_snp_onco_313_",i,".txt")
  qctool.command[i,1] <- temp
  
}

write.table(qctool.command,file = paste0("./risk_prediction/Nasim_prs/code/qc_extract_snp_onco.sh"),col.names = F,row.names = F,quote=F)



#extract the genotype from icogs data
extract_snp <-  snp_id %>%
  select(SNP.ICOGS)

write.table(extract_snp,file = paste0("./risk_prediction/Nasim_prs/result/extract_snp_313_icog.txt"),row.names=F,col.names=T,quote=F)

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
  
  
  
  temp <- paste0("/data/zhangh24/qctool_v1.4-linux-x86_64/qctool -g ",Files[i]," -incl-rsids /data/zhangh24/breast_cancer_data_analysis/risk_prediction/Nasim_prs/result/extract_snp_313_icog.txt -og /data/zhangh24/breast_cancer_data_analysis/risk_prediction/Nasim_prs/result/extracted_snp_icog_313_",i,".txt")
  qctool.command[i,1] <- temp
  
}


write.table(qctool.command,file = paste0("./risk_prediction/Nasim_prs/code/qc_extract_snp_icog.sh"),col.names = F,row.names = F,quote=F)


