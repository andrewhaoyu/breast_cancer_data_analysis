#Goal estimate heritability using standard log odds ratio

#prepare the data for ld score regression file
library(data.table)
library(dplyr)
setwd('/dcl01/chatterj/data/hzhang1/breast_cancer_data_analysis')
standard_result <- as.data.frame(fread("./discovery_SNP/result/ResultsMeta_GWAS_iCOGs_Onco_filter_R2_MAF.txt"))
#all_result <- fread("./data/oncoarray_bcac_public_release_oct17.txt",header=T)

temp.str <- strsplit(standard_result$SNP.Onco,":")
n <- nrow(standard_result)
snpid <- rep("c",n)
for(i in 1:n){
  snpid[i] <- temp.str[[i]][1]
}
idx <- which(c(1:n)%in%grep("rs",snpid)!=T)
snpid[idx] <- paste0("chr",standard_result$chr.Onco,":",
                     standard_result$Position.Onco)[idx]
standard_result$snp.id <- snpid

#onco array data
Twometa <- function(beta1,var1,beta2,var2){
  var_meta <- 1/(1/var1+1/var2)
  beta_meta <- (var_meta)*(1/var1*beta1+
                               1/var2*beta2)
  return(list(beta_meta,var_meta))
}

standard_result = standard_result %>% 
  mutate(BCAC_meta_beta = Twometa(beta.iCOGs,SE.iCOGs^2,beta.Onco,SE.Onco^2)[[1]],
         BCAC_meta_var = Twometa(beta.iCOGs,SE.iCOGs^2,beta.Onco,SE.Onco^2)[[2]])



onco_result <- standard_result %>% mutate(
  or = exp(beta.Onco),
  sample_size = 1/(SE.Onco^2*2*EAFcontrols.Onco*(1-EAFcontrols.Onco))) %>% 
  select(snp.id,chr.Onco,Position.Onco,Effect.Onco,Baseline.Onco,or,SE.Onco,P1df_risk_chi.Onco,
         r2.Onco,
         EAFcontrols.Onco,
         sample_size)
colnames(onco_result) <- c("snpid",
                           "CHR",
                           "bp",
                           "A1",
                           "A2",
                           "or",
                           "se",
                           "P",
                           "info",
                           "freq",
                           "N")
write.table(onco_result,file="/dcl01/chatterj/data/hzhang1/ldsc/onco_result.txt",col.names = T,quote=F)
conda env create --file environment.yml
source activate ldsc
./munge_sumstats.py --sumstats onco_result.txt --out onco --merge-alleles w_hm3.snplist --info-min 0.3
./ldsc.py --h2 onco.sumstats.gz --ref-ld-chr eur_w_ld_chr/ --w-ld-chr eur_w_ld_chr/ --out onco
less onco.log
#estimate 0.4227
#write.table(meta_result,file="./heritability/result/meta_result.txt",col.names = T,quote=F)








bcac_result <- standard_result %>% mutate(
  z = BCAC_meta_beta/sqrt(BCAC_meta_var),
  sample_size = 1/(BCAC_meta_var*2*EAFcontrols.Onco*(1-EAFcontrols.Onco)),
  P = 2*pnorm(-abs(z))) %>% 
  select(snp.id,chr.Onco,Position.Onco,Effect.Onco,Baseline.Onco,z,P,
         r2.Onco,
         EAFcontrols.Onco,
         sample_size)
colnames(bcac_result) <- c("snpid",
                           "CHR",
                           "bp",
                           "A2",
                           "A1",
                           "Z",
                           "P",
                           "info",
                           "freq",
                           "N")
write.table(bcac_result,file="/dcl01/chatterj/data/hzhang1/ldsc/bcac_result.txt",col.names = T,quote=F)
library(data.table)
bcac_result <- as.data.frame(fread("/dcl01/chatterj/data/hzhang1/ldsc/bcac_result.txt"))

./munge_sumstats.py --sumstats bcac_result.txt --out bcac --merge-alleles w_hm3.snplist --info-min 0.3  --signed-sumstats Z,0 --frq freq
./ldsc.py --h2 bcac.sumstats.gz --ref-ld-chr eur_w_ld_chr/ --w-ld-chr eur_w_ld_chr/ --out bcac
less bcac.log
#estimate 0.4777

#load genetic correlation data sent to Gunaghao






######meta result based on icogs, onco and gwas
# meta_result <- standard_result %>% mutate(
#   z = Beta.meta/sqrt(var.meta),
#   sample_size = 1/(var.meta*2*EAFcontrols.Onco*(1-EAFcontrols.Onco))
# ) %>% 
#   select(snp.id,chr.Onco,Position.Onco,Effect.Meta,Baseline.Meta,z,p.meta,
#          r2.Onco,
#          EAFcontrols.Onco,
#          sample_size)
# colnames(meta_result) <- c("snpid",
#                            "CHR",
#                            "bp",
#                            "A2",
#                            "A1",
#                            "Z",
#                            "P-value",
#                            "info",
#                            "freq",
#                            "N")
# write.table(meta_result,file="./heritability/result/meta_result.txt",col.names = T,quote=F)
# conda env create --file environment.yml
# source activate ldsc
# ./munge_sumstats.py --sumstats meta_result.txt --out meta --merge-alleles w_hm3.snplist --info-min 0.3  --signed-sumstats Z,0 --frq freq
# ./ldsc.py --h2 meta.sumstats.gz --ref-ld-chr eur_w_ld_chr/ --w-ld-chr eur_w_ld_chr/ --out meta
# less meta.log
# 
# 
# 
# 
# 
# 
# all_result <- meta_result %>% 
#   mutate(chrpos=paste0(chr,":",position_b37))
# KG <- fread("./data/KG.all.chr.bim")
# KG <- KG %>% mutate(chrpos = paste0(V1,":",V4))
# all_result <- merge(all_result,KG,by = "chrpos")
# all_result <- all_result %>% 
#   mutate(meta_or = exp(as.numeric(bcac_onco_icogs_gwas_beta)),
#          icog_or = exp(as.numeric(bcac_icogs2_beta)),
#          onco_or = exp(as.numeric(bcac_onco2_beta)),
#          gwas_or = exp(as.numeric(bcac_gwas_all_beta)))
# icog_result <- all_result %>% 
#   select(V2,chr,position_b37,a0,a1,icog_or,bcac_icogs2_se,
#          bcac_icogs2_P1df_Wald,
#          bcac_onco2_r2,
#          bcac_onco_icogs_gwas_eaf_controls)
# colnames(icog_result) <- c("snpid",
#                            "CHR",
#                            "bp",
#                            "a1",
#                            "a2",
#                            "or",
#                            "se",
#                            "pval",
#                            "info",
#                            "freq")
# write.table(icog_result,file="./heritability/result/icog_result.txt",col.names = T,quote=F)
# ./munge_sumstats.py --sumstats meta_result_new.txt --out meta_new --merge-alleles w_hm3.snplist --info-min 0.3 --frq freq 
# ./ldsc.py --h2 meta_new.sumstats.gz --ref-ld-chr eur_w_ld_chr/ --w-ld-chr eur_w_ld_chr/ --out meta_new
# less meta_new.log
# #estimate 0.5088
# 
# 
# 
# ######meta result
# meta_result <- standard_result %>% mutate(
#   or = exp(Beta.meta),
#   sample_size = 1/(var.meta*2*EAFcontrols.Onco*(1-EAFcontrols.Onco))
# ) %>% 
#   select(snp.id,chr.Onco,Position.Onco,Effect.Meta,Baseline.Meta,or,sdE.meta,p.meta,
#          r2.Onco,
#          EAFcontrols.Onco,
#          sample_size)
# colnames(meta_result) <- c("snpid",
#                            "CHR",
#                            "bp",
#                            "A1",
#                            "A2",
#                            "or",
#                            "se",
#                            "P-value",
#                            "info",
#                            "freq",
#                            "N")
# write.table(meta_result,file="./heritability/result/meta_result_new.txt",col.names = T,quote=F)
# 
# 
# 
# 
# 
# 
# 
# 
# 
# ###############LD score regression was run locally with ldsc
# ###############iCOGs 46785 cases 42892 controls, effective 22377
# ###############Onco Array 61282 cases 45494 controls, effective 26110.39
# ###############GWAS 14910 cases 17588 controls, effective 8069.33
# cd /Users/zhangh24/GoogleDrive/ldsc
# conda env create --file environment.yml
# source activate ldsc
# ./ldsc.py -h
# ./munge_sumstats.py -h
# ./munge_sumstats.py --sumstats icog_result.txt --N 22377 --out icog --merge-alleles w_hm3.snplist --info-min 0.3
# ./ldsc.py --h2 icog.sumstats.gz --ref-ld-chr eur_w_ld_chr/ --w-ld-chr eur_w_ld_chr/ --out icog
# less icog.log
# ######onco result
# onco_result <- all_result %>% 
#   select(V2,chr,position_b37,a0,a1,onco_or,bcac_onco2_se,
#          bcac_onco2_P1df_Wald,
#          bcac_onco2_r2,
#          bcac_onco_icogs_gwas_eaf_controls)
# colnames(onco_result) <- c("snpid",
#                            "CHR",
#                            "bp",
#                            "a1",
#                            "a2",
#                            "or",
#                            "se",
#                            "pval",
#                            "info",
#                            "freq")
# write.table(onco_result,file="./heritability/result/onco_result.txt",col.names = T,quote=F)
# ./munge_sumstats.py --sumstats onco_result.txt --N 26110.39 --out onco --merge-alleles w_hm3.snplist --info-min 0.3
# ./ldsc.py --h2 onco.sumstats.gz --ref-ld-chr eur_w_ld_chr/ --w-ld-chr eur_w_ld_chr/ --out onco
# less onco.log

