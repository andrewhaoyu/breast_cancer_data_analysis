#Goal estimate heritability using standard log odds ratio

#prepare the data for ld score regression file
library(data.table)
library(tidyverse)
setwd("/spin1/users/zhangh24/breast_cancer_data_analysis/")
all_result <- fread("./data/oncoarray_bcac_public_release_oct17.txt",header=T)
all_result <- meta_result %>% 
  mutate(chrpos=paste0(chr,":",position_b37))
KG <- fread("./data/KG.all.chr.bim")
KG <- KG %>% mutate(chrpos = paste0(V1,":",V4))
all_result <- merge(all_result,KG,by = "chrpos")
all_result <- all_result %>% 
  mutate(meta_or = exp(as.numeric(bcac_onco_icogs_gwas_beta)),
         icog_or = exp(as.numeric(bcac_icogs2_beta)),
         onco_or = exp(as.numeric(bcac_onco2_beta)),
         gwas_or = exp(as.numeric(bcac_gwas_all_beta)))
icog_result <- all_result %>% 
  select(V2,chr,position_b37,a0,a1,icog_or,bcac_icogs2_se,
         bcac_icogs2_P1df_Wald,
         bcac_onco2_r2,
         bcac_onco_icogs_gwas_eaf_controls)
colnames(icog_result) <- c("snpid",
                           "CHR",
                           "bp",
                           "a1",
                           "a2",
                           "or",
                           "se",
                           "pval",
                           "info",
                           "freq")
write.table(icog_result,file="./heritability/result/icog_result.txt",col.names = T,quote=F)
###############LD score regression was run locally with ldsc
