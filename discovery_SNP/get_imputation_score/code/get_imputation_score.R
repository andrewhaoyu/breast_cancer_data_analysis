#get the imputation score for 32 discovery SNPs


setwd("/spin1/users/zhangh24/breast_cancer_data_analysis")
load("/data/zhangh24/breast_cancer_data_analysis/whole_genome/ONCO/ERPRHER2GRADE_fixed_baseline/result/onco_info.Rdata")
load("/data/zhangh24/breast_cancer_data_analysis/whole_genome/ICOG/ERPRHER2GRADE_fixed_baseline/result/icog_info.Rdata")
discovery_snp_info <- read.csv("./data/discovery_snp_paper_order.csv")
discovery_snp_info_sub = discovery_snp_info %>% 
  mutate(chr.pos = paste0(CHR,":",position))
library(dplyr)
load("/data/zhangh24/match.Rdata")
onco_info_sub = onco_info %>% select(rs_id,CHR,position,info,exp_freq_a1)
colnames(onco_info_sub) <- c("SNP.ONCO","CHR_onco","position_onco","info_onco","freq_a1_onco")
icog_info_sub = icog_info %>% select(rs_id,CHR,position,info,exp_freq_a1)
colnames(icog_info_sub) <- c("SNP.ICOGS","CHR_icog","position_icog","info_icog","freq_a1_icog")


icog_info_sub_temp <- full_join(icog_info_sub,
                                data,
                                by="SNP.ICOGS")

icog_onco_infor <- full_join(icog_info_sub_temp,
                             onco_info_sub,
                             by="SNP.ONCO")
save(icog_onco_infor,file = "/data/zhangh24/icog_onco_information_data.rdata")
head(icog_onco_infor)

icog_onco_info_sub = icog_onco_infor %>% 
  mutate(chr.pos=paste0(CHR_onco,":",position_onco))



discovery_snp_info_sub <- left_join(discovery_snp_info_sub,icog_onco_info_sub,by="chr.pos")
write.csv(discovery_snp_info_sub,file="")

