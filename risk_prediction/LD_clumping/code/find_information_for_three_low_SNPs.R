#-------------------------------------------------------------------
# Update Date: 01/22/2020
# Create Date: 01/22/2020
# Goal: Three of the 313 SNPs are low frequency list. Find them and add them into the whole_genome_result
# Author: Haoyu Zhang
#-------------------------------------------------------------------
setwd('/data/zhangh24/breast_cancer_data_analysis/risk_prediction/')
load("./intrinsic_subtypes_whole_genome/ICOG/result/whole_gonome.rdata")

#load("./EB_whole_genome/result/whole_gonome.rdata")


library(dplyr)
n <- nrow(whole_genome)
assoc <- whole_genome %>% mutate(p.min = pmin(stan_p,FTOP_result),
                                 TEST = rep("ADD",n),
                                 NMISS = rep(0,n),
                                 OR = exp(stan_logodds),
                                 STAT = rnorm(n)) %>%
  select(CHR,SNP.ONCO,position,
         effect_allele,
         TEST,NMISS,
         OR,STAT,p.min,
         var_name)

colnames(assoc) <- c("CHR",
                     "SNP",
                     "BP",
                     "A1",
                     "TEST",
                     "NMISS",
                     "OR",
                     "STAT",
                     "P",
                     "var_name")
#load 313 Nasim SNPs
setwd('/data/zhangh24/breast_cancer_data_analysis/')
library(dplyr)
snp <- read.csv("./data/Nasim_313_SNPs_infor.csv",header=T)
load("/data/zhangh24/match.Rdata")
colnames(snp)[1] <- "var_name"
snp.new <- left_join(snp,data, by="var_name")
idx <- which((snp.new$var_name%in%assoc$var_name)!=T)
length(idx)
#three SNPs were not in the whole_gonome.rdata summary level statistics due to low allele frequency. Add back to it manually
load(paste0("./risk_prediction/Nasim_prs/result/313_intrinsic_subtype_logodds_var.Rdata"))
intrinsic_subtypes_result <- final_result
load(paste0("./risk_prediction/Nasim_prs/result/313_overall_logodds_var.Rdata"))
overall_result <- final_result
snp.new <- cbind(snp.new,intrinsic_subtypes_result[,4:33],overall_result[,4:5])
#load the icog and onco summary level statistics to get freq_a1_icog，freq_a1_onco
load("./risk_prediction/intrinsic_subtypes_whole_genome/ICOG/result/Icog_result_intrinsic_subtype.Rdata")
load("./risk_prediction/intrinsic_subtypes_whole_genome/ONCO/result/onco_result_intrinsic_subtype.Rdata")
icog_result <- icog_result_casecase
onco_result <- onco_result_casecase

icog_result = icog_result %>% 
  select(rs_id,exp_freq_a1,
         info) %>% 
  rename(SNP.ICOGS=rs_id,
         freq_a1_icog=exp_freq_a1,
         info_icog=info)
snp.new <- left_join(snp.new,icog_result,
                     by="SNP.ICOGS")
onco_result = onco_result %>% 
  select(rs_id,exp_freq_a1,
         info) %>% 
  rename(SNP.ONCO=rs_id,
         freq_a1_onco=exp_freq_a1,
         info_onco=info)
snp.new <- left_join(snp.new,onco_result,
                     by="SNP.ONCO")
#create the information for three SNPs to save in whole_genome_analysis
snp_id <- c("---",
            "---",
            "---")
rs_id <- snp.new$SNP.ONCO[idx]
position <- snp.new$Position[idx]
exp_freq_a1 <- snp.new$freq_a1_icog[idx]
info <- snp.new$info_icog[idx]
certainty <- rep(1,3)
type <- rep(0,3)
info_type0 <- rep(-1,3)
concord_type0 <- rep(-1,3)
r2_type0 <- rep(-1,3)
CHR <- snp.new$Chromosome[idx]
var_name <- snp.new$var_name[idx]
SNP.ICOGS <-snp.new$SNP.ICOGS[idx]
SNP.ONCO <- snp.new$SNP.ONCO[idx]
stan_logodds <- snp.new$logodds_overall[idx]
stan_var <- snp.new$var_overall[idx]
stan_p <- rep(0,3)
FTOP_rsult <- rep(0,3)
intrinsic_log_var = snp.new[idx,16:45]
in_p <- rep(0,3)
reference_allele = snp.new$Reference.Allele[idx]
effect_allele = snp.new$Effect.Allele[idx]
p.min <- rep(0,3)
freq_a1_icog = snp.new$freq_a1_icog[idx]
info_icog = snp.new$info_icog[idx]
freq_a1_onco = snp.new$freq_a1_onco[idx]
info_onco = snp.new$info_onco[idx]
three.snps <- data.frame(snp_id,
                         rs_id,
                         position,
                         exp_freq_a1,
                         info,
                         certainty,
                         type,
                         info_type0,
                         concord_type0,
                         r2_type0,
                         CHR,
                         var_name,
                         SNP.ICOGS,
                         SNP.ONCO,
                         stan_logodds,
                         stan_var,
                         stan_p,
                         FTOP_rsult,
                         intrinsic_log_var,
                         in_p,
                         reference_allele,
                         effect_allele,
                         p.min,
                         freq_a1_icog,
                         info_icog,
                         freq_a1_onco,
                         info_onco)
colnames(three.snps) = colnames(whole_genome)
whole_genome = rbind(whole_genome,three.snps)
setwd('/data/zhangh24/breast_cancer_data_analysis/risk_prediction/')
save(whole_genome,file = "./intrinsic_subtypes_whole_genome/ICOG/result/whole_genome_threeadd.rdata")
setwd('/data/zhangh24/breast_cancer_data_analysis/')

save(snp.new, file = "./data/Nasim_313SNPs_complete_information.Rdata")
#load whole genome information to get the SNP id in OncoArray and iCOGs



#write.table(assoc,file="/data/zhangh24/BCAC/impute_plink_onco/LD_assoc",col.names = T, row.names=F, quote = F)

