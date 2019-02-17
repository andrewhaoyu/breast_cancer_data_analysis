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
###############iCOGs 46785 cases 42892 controls, effective 22377
###############Onco Array 61282 cases 45494 controls, effective 26110.39
###############GWAS 14910 cases 17588 controls, effective 8069.33
cd /Users/zhangh24/GoogleDrive/ldsc
conda env create --file environment.yml
source activate ldsc
./ldsc.py -h
./munge_sumstats.py -h
./munge_sumstats.py --sumstats icog_result.txt --N 22377 --out icog --merge-alleles w_hm3.snplist --info-min 0.3
./ldsc.py --h2 icog.sumstats.gz --ref-ld-chr eur_w_ld_chr/ --w-ld-chr eur_w_ld_chr/ --out icog
less icog.log
######onco result
onco_result <- all_result %>% 
  select(V2,chr,position_b37,a0,a1,onco_or,bcac_onco2_se,
         bcac_onco2_P1df_Wald,
         bcac_onco2_r2,
         bcac_onco_icogs_gwas_eaf_controls)
colnames(onco_result) <- c("snpid",
                           "CHR",
                           "bp",
                           "a1",
                           "a2",
                           "or",
                           "se",
                           "pval",
                           "info",
                           "freq")
write.table(onco_result,file="./heritability/result/onco_result.txt",col.names = T,quote=F)
./munge_sumstats.py --sumstats onco_result.txt --N 26110.39 --out onco --merge-alleles w_hm3.snplist --info-min 0.3
./ldsc.py --h2 onco.sumstats.gz --ref-ld-chr eur_w_ld_chr/ --w-ld-chr eur_w_ld_chr/ --out onco
less onco.log

######meta result
meta_result <- all_result %>% 
  select(V2,chr,position_b37,a0,a1,meta_or,bcac_onco_icogs_gwas_se ,
         bcac_onco_icogs_gwas_P1df,
         bcac_onco2_r2,
         bcac_onco_icogs_gwas_eaf_controls)
colnames(meta_result) <- c("snpid",
                           "CHR",
                           "bp",
                           "a1",
                           "a2",
                           "or",
                           "se",
                           "pval",
                           "info",
                           "freq")
write.table(meta_result,file="./heritability/result/meta_result.txt",col.names = T,quote=F)
./munge_sumstats.py --sumstats meta_result.txt --N 56922.07 --out meta --merge-alleles w_hm3.snplist --info-min 0.3
./ldsc.py --h2 meta.sumstats.gz --ref-ld-chr eur_w_ld_chr/ --w-ld-chr eur_w_ld_chr/ --out meta
less meta.log
