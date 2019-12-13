#organize Nasims and new 32 SNPs
#load Nasim's snps information
snp <- read.csv("./data/Nasim_313_SNPs_infor.csv",header=T)
library(dplyr)
snp = snp %>% mutate(chr.pos = paste0(CHR,":",position))
#load whole genome information to get the SNP id in OncoArray and iCOGs
load("./risk_prediction/intrinsic_subtypes_whole_genome/ICOG/result/whole_gonome.rdata")

#load Nasim's SNPs onco array data
load("/data/zhangh24/breast_cancer_data_analysis/risk_prediction/Nasim_prs/result/onco.nasim.snp.rdata")
#load Nasim's SNPs icogs data
load("/data/zhangh24/breast_cancer_data_analysis/risk_prediction/Nasim_prs/result/icog.nasim.snp")