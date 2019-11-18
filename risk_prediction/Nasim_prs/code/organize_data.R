#organize the Nasim snp list
setwd("/data/zhangh24/breast_cancer_data_analysis/")
snp <- read.csv("./data/Nasim_313_SNPs_infor.csv",header=T)
load("./intrinsic_subtypes_whole_genome/ICOG/result/whole_gonome.rdata")
library(dplyr)
colnames(snp)[1] <- "var_name"
snp_all <- left_join(snp,whole_genome,
                     by="var_name")
snp_odds = snp_all %>% 
  select(var_name,rs_id,SNP.ICOGS,SNP.ONCO,CHR,position,reference_allele,effect_allele,freq_a1_icog,freq_a1_onco,Overall.Breast.Cancerd,ER.positivee,ER.negativef,stan_logodds,Luminal_A,Luminal_B,Luminal_B_HER2Neg,HER2_Enriched,TN,info_icog,info_onco)