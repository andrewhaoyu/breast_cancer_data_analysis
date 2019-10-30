#Goal: prepare data set into the format of LD-score regression
#load icogs results
load("/Users/zhangh24/GoogleDrive/breast_cancer_data_analysis/whole_genome_age/ICOG/Intrinsic_subtypes/result/intrinsic_subtype_triple_negative_results.Rdata")
load("/Users/zhangh24/GoogleDrive/breast_cancer_data_analysis/genetic_correlation/standard_analysis/result/standard_gwas_result/case_control/ICOG/icog_result_odds_sd.Rdata")
library(data.table)
icog_result <- icog_result %>% mutate(CHRpos=paste0(CHR,":",bp))
library(dplyr)

new_result <- merge(icog_result,intrinsic_subtype_triple_negative_results,by="rs_id")
#change the data format into the LD score regression format
icog_result <- new_result  %>% 
  mutate(or = exp(logodds)) %>% 
  mutate(p =  2*pnorm(-abs(logodds/sd)))%>% 
  mutate(ngt=0) %>% 
  select(rs_id,CHR,position.x,allele1,allele2,or,sd,p,imputation_quality,ngt,freq_a1) 
colnames(icog_result) <- c("snpid",
                           "CHR",
                           "bp",
                           "a1",
                           "a2",
                           "or",
                           "se",
                           "pval",
                           "info",
                           "ngt",
                           "freq")
KG <- fread("/Users/zhangh24/GoogleDrive/breast_cancer_data_analysis/data/KG.all.chr.bim")
KG <- KG %>% mutate(CHRpos = paste0(V1,":",V4))
new_result <- merge(icog_result,KG,by = "CHRpos")
icog_result <- new_result %>% select(V2,CHR,bp,a1,a2,or,se,pval,info,ngt,freq)
colnames(icog_result) <- c("snpid",
                           "CHR",
                           "bp",
                           "a1",
                           "a2",
                           "or",
                           "se",
                           "pval",
                           "info",
                           "ngt",
                           "freq")

write.table(icog_result,file = "/Users/zhangh24/GoogleDrive/breast_cancer_data_analysis/genetic_correlation/standard_analysis/result/standard_gwas_result/case_control/ICOG/icog_result_standard_analysis.txt",col.names = T,quote = F)
#p = 2*pnorm(-abs(-1.96))

#load onco array results
load("/Users/zhangh24/GoogleDrive/breast_cancer_data_analysis/genetic_correlation/standard_analysis/result/standard_gwas_result/case_control/ONCO/onco_result_odds_sd.Rdata")
#onco_result <- onco_result %>% mutate(CHRpos=paste0(CHR,":",bp))
new_result <- merge(onco_result,intrinsic_subtype_triple_negative_results,by="rs_id")
onco_result <- new_result  %>% 
  mutate(or = exp(logodds)) %>% 
  mutate(p =  2*pnorm(-abs(logodds/sd)))%>% 
  mutate(ngt=0) %>% 
  select(rs_id,CHR,position.x,allele1,allele2,or,sd,p,imputation_quality,ngt,freq_a1)
colnames(onco_result) <- c("snpid",
                           "CHR",
                           "bp",
                           "a1",
                           "a2",
                           "or",
                           "se",
                           "pval",
                           "info",
                           "ngt",
                           "freq")
onco_result <- onco_result %>% mutate(CHRpos=paste0(CHR,":",bp))
new_result <- merge(onco_result,KG,by = "CHRpos")
onco_result <- new_result %>% select(V2,CHR,bp,a1,a2,or,se,pval,info,ngt,freq)
colnames(onco_result) <- c("snpid",
                           "CHR",
                           "bp",
                           "a1",
                           "a2",
                           "or",
                           "se",
                           "pval",
                           "info",
                           "ngt",
                           "freq")
write.table(onco_result,file = "/Users/zhangh24/GoogleDrive/breast_cancer_data_analysis/genetic_correlation/standard_analysis/result/standard_gwas_result/case_control/ONCO/onco_result_standard_analysis.txt",col.names=T,quote=F)

#load the meta-analysis results
load("/Users/zhangh24/GoogleDrive/breast_cancer_data_analysis/genetic_correlation/standard_analysis/result/standard_gwas_result/case_control/meta/case_control_meta_result_final.Rdata")
meta_result <- case_control_meta_result_final


new_result <- merge(meta_result,intrinsic_subtype_triple_negative_results,by="rs_id")

meta_result <- new_result  %>% 
  mutate(or = exp(logodds_meta)) %>% 
  mutate(p =  2*pnorm(-abs(logodds_meta/sd_meta)))%>% 
  mutate(ngt=0) %>% 
  select(rs_id,CHR.x,position.x,allele1,allele2,or,sd_meta,p,imputation_quality,ngt,freq_a1)
colnames(meta_result) <- c("snpid",
                           "CHR",
                           "bp",
                           "a1",
                           "a2",
                           "or",
                           "se",
                           "pval",
                           "info",
                           "ngt",
                           "freq")
meta_result <- meta_result %>% mutate(CHRpos=paste0(CHR,":",bp))
new_result <- merge(meta_result,KG,by = "CHRpos")
meta_result <- new_result %>% select(V2,CHR,bp,a1,a2,or,se,pval,info,ngt,freq)
colnames(meta_result) <- c("snpid",
                           "CHR",
                           "bp",
                           "a1",
                           "a2",
                           "or",
                           "se",
                           "pval",
                           "info",
                           "ngt",
                           "freq")
write.table(meta_result,file = "/Users/zhangh24/GoogleDrive/breast_cancer_data_analysis/genetic_correlation/standard_analysis/result/standard_gwas_result/case_control/meta/meta_result_standard_analysis.txt",col.names=T,quote=F)

#load onco array results
library(data.table)
library(tidyverse)
setwd("/data/zhangh24/breast_cancer_data_analysis/")
meta_result <- read.table("./data/oncoarray_bcac_public_release_oct17.txt",header=T)
meta_result <- meta_result %>% 
  mutate(chrpos=paste0(chr,":",position_b37))

