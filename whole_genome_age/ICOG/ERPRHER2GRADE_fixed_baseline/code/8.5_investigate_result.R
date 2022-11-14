load("/data/NC_BW/HZ_SF/data3.rdata")
#load("data.rdata")
#p.value.x is the MTOP pairwise p-value
#p.value.y is the MTOP saturated p-value
#p.value is the FTOP p-value
load("data_clean.rdata")
library(dplyr)
data_clean_select = data_clean %>% 
  select(var_name, SNP.ICOGS.x, SNP.ONCO.x, position.x, exp_freq_a1.x, 
         info.x, certainty.x, CHR.x, p.value.x, p.value.y, p.value, p.acat) %>% 
  rename(SNP.ICOGS = SNP.ICOGS.x,
         SNP.ONCO = SNP.ONCO.x,
         position = position.x,
         exp_freq_a1 = exp_freq_a1.x,
         info = info.x,
         certainty = certainty.x,
         CHR = CHR.x,
         p.pairwise = p.value.x,
         p.saturated = p.value.y,
         p.FTOP = p.value)

data_clean_select_sig = data_clean_select %>% 
  filter(p.acat<=5E-08)
tail(data_clean_select_sig)
library(data.table)
known_snp = as.data.frame(fread("/data/zhangh24/breast_cancer_data_analysis/data/210_known_discovery_snp_paper_order.csv"))
head(known_snp)


#filter SNPs within +-2Mb basepair of known SNPs
# idx <- which(data_clean_select_sig$CHR==known_snp$Chr&
#                (data_clean_select_sig$position>=known_snp$Position-2*10^6)&
#                (data_clean_select_sig$position<=known_snp$Position+2*10^6))

n.known = nrow(known_snp)
idx = NULL
for(i in 1:n.known){
  CHRi = known_snp$CHR[i]
  lb = known_snp$position[i]-2*10^6
  ub = known_snp$position[i]+2*10^6
  temp =  which((data_clean_select_sig$CHR==CHRi) & (data_clean_select_sig$position>=lb) & (data_clean_select_sig$position<=ub))
  idx = c(idx,temp)
  cat(i,length(idx),"\n")
}

idx = unique(idx)
length(idx)  
#
data_clean_select_sig_best= data_clean_select_sig[-idx,]
library(data.table)

head(data)
idx <- which(data$chr.Onco==3&data$Position.Onco==115778273)
data[idx,]
 data = fread("/data/zhangh24/breast_cancer_data_analysis/discovery_SNP/prepare_summary_level_statistics/result/icogs_onco_gwas_meta_overall_breast_cancer_summary_level_statistics.txt")
# data2 = fread("/data/zhangh24/breast_cancer_data_analysis/discovery_SNP/prepare_summary_level_statistics/result/CIMBA_BRCA1_BCAC_TN_meta_summary_level_statistics.txt")

 
 # setwd("/data/NC_BW/HZ_SF")
 # load("data_clean_select_sig_best_p3.rdata")
 # data_clean_select_sig_best1 = data_clean_select_sig_best
 # load("data_clean_select_sig_best_p4.rdata")
 # data_clean_select_sig_best2 = data_clean_select_sig_best
 # load("data_clean_select_sig_best_p5.rdata")
 # data_clean_select_sig_best3 = data_clean_select_sig_best
 # 
 
#        position = as.numeric(position)) %>% 
  # filter(CHR==22&(position-28195386)<=500000&(position-28195386)>=-500000)
CHR = data_select$CHR
position= data_select$position
idx <- which(CHR==22&(position-28195386)<=500000&(position-28195386)>=-500000)

