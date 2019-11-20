#organize the Nasim snp list
setwd("/data/zhangh24/breast_cancer_data_analysis/")
snp <- read.csv("./data/Nasim_313_SNPs_infor.csv",header=T)
load("./risk_prediction/intrinsic_subtypes_whole_genome/ICOG/result/whole_gonome.rdata")

library(dplyr)
colnames(snp)[1] <- "var_name"
snp_all <- left_join(snp,whole_genome,
                     by="var_name")
snp_odds = snp_all %>% 
  select(var_name,rs_id,SNP.ICOGS,SNP.ONCO,CHR,position,reference_allele,effect_allele,freq_a1_icog,freq_a1_onco,Overall.Breast.Cancerd,ER.positivee,ER.negativef,stan_logodds,Luminal_A,Luminal_B,Luminal_B_HER2Neg,HER2_Enriched,TN,info_icog,info_onco,p.min)
save(snp_odds,file = paste0("./risk_prediction/Nasim_prs/result/nasim_snp_odds.rdata"))


load("/data/zhangh24/match.Rdata")
snp_id = left_join(snp,data,
                   by="var_name") %>% 
  select(var_name,SNP.ICOGS,SNP.ONCO)
idx <- which(is.na(snp_id$SNP.ONCO))
save(snp_id,file = paste0("./risk_prediction/Nasim_prs/result/nasim_snp_id.rdata"))
