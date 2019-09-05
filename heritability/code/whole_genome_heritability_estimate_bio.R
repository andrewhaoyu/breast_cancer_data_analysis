#Goal estimate heritability using standard log odds ratio

#prepare the data for ld score regression file
library(data.table)
library(dplyr)
setwd('/spin1/users/zhangh24/breast_cancer_data_analysis/')
standard_result <- as.data.frame(fread("./discovery_SNP/result/ResultsMeta_GWAS_iCOGs_Onco_filter_R2_MAF.txt"))
#all_result <- fread("./data/oncoarray_bcac_public_release_oct17.txt",header=T)
#create snp.id with official name
#SNPs with rs number will be listed as rs SNP
#SNPs without rs number will be listed as chr:pos
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
#generate meta-analysis odds ratio for icog and onco
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
write.table(onco_result,file="/spin1/users/zhangh24/ldsc/onco_result.txt",col.names = T,quote=F)
./munge_sumstats.py --sumstats onco_result.txt --out onco --merge-alleles w_hm3.snplist --info-min 0.3
./ldsc.py --h2 onco.sumstats.gz --ref-ld-chr eur_w_ld_chr/ --w-ld-chr eur_w_ld_chr/ --out onco
less onco.log

#write.table(meta_result,file="./heritability/result/meta_result.txt",col.names = T,quote=F)








bcac_result <- standard_result %>% mutate(
  z = BCAC_meta_beta/sqrt(BCAC_meta_var),
 # sample_size = 1/(BCAC_meta_var*2*EAFcontrols.Onco*(1-EAFcontrols.Onco)),
  P = 2*pnorm(-abs(z))) %>% 
  select(snp.id,chr.Onco,Position.Onco,Effect.Onco,Baseline.Onco,z,P,
         r2.Onco,
         EAFcontrols.Onco)
colnames(bcac_result) <- c("snpid",
                           "CHR",
                           "bp",
                           "A2",
                           "A1",
                           "Z",
                           "P",
                           "info",
                           "freq")
write.table(bcac_result,file="/spin1/users/zhangh24/ldsc/bcac_result.txt",col.names = T,quote=F)

conda env create --file environment.yml
source activate ldsc
./munge_sumstats.py --sumstats bcac_result.txt --out meta --merge-alleles w_hm3.snplist --info-min 0.3  --signed-sumstats Z,0 --frq freq  --N 53091.03
./ldsc.py --h2 meta.sumstats.gz --ref-ld-chr eur_w_ld_chr/ --w-ld-chr eur_w_ld_chr/ --out meta
less meta.log

cimba <- as.data.frame(fread("/spin1/users/zhangh24/breast_cancer_data_analysis/data/brca1_bcac_tn_meta.txt"))
idx <- which(cimba$CHR==6&cimba$position==33239869)
cimba[idx,]
#sample size calculated
bcac_result <- standard_result %>% mutate(
  z = BCAC_meta_beta/sqrt(BCAC_meta_var),
  sample_size = 1/(BCAC_meta_var*2*EAFcontrols.Onco*(1-EAFcontrols.Onco)),
  P = 2*pnorm(-abs(z)),
  or = exp(BCAC_meta_beta),
  sebcac = sqrt(BCAC_meta_var)) %>% 
  select(snp.id,chr.Onco,Position.Onco,Effect.Onco,Baseline.Onco,or,sebcac,P,
         r2.Onco,
         EAFcontrols.Onco,
         sample_size)
colnames(bcac_result) <- c("snpid",
                           "CHR",
                           "bp",
                           "A2",
                           "A1",
                           "or",
                           "se",
                           "P",
                           "info",
                           "freq",
                           "N")
write.table(bcac_result,file="/spin1/users/zhangh24/ldsc/bcac_result_new.txt",col.names = T,quote=F)

# conda env create --file environment.yml
# source activate ldsc
./munge_sumstats.py --sumstats bcac_result_new.txt --out meta_new --merge-alleles w_hm3.snplist --info-min 0.3 --frq freq
./ldsc.py --h2 meta_new.sumstats.gz --ref-ld-chr eur_w_ld_chr/ --w-ld-chr eur_w_ld_chr/ --out meta_new
less meta_new.log

#sample size plug in
bcac_result <- standard_result %>% mutate(
  z = BCAC_meta_beta/sqrt(BCAC_meta_var),
  sample_size = 1/(BCAC_meta_var*2*EAFcontrols.Onco*(1-EAFcontrols.Onco)),
  P = 2*pnorm(-abs(z)),
  or = exp(BCAC_meta_beta),
  sebcac = sqrt(BCAC_meta_var)) %>% 
  select(snp.id,chr.Onco,Position.Onco,Effect.Onco,Baseline.Onco,or,sebcac,P,
         r2.Onco,
         EAFcontrols.Onco)
colnames(bcac_result) <- c("snpid",
                           "CHR",
                           "bp",
                           "A2",
                           "A1",
                           "or",
                           "se",
                           "P",
                           "info",
                           "freq")
write.table(bcac_result,file="/spin1/users/zhangh24/ldsc/bcac_result_new2.txt",col.names = T,quote=F)

# conda env create --file environment.yml
# source activate ldsc
./munge_sumstats.py --sumstats bcac_result_new2.txt --out meta_new2 --merge-alleles w_hm3.snplist --info-min 0.3 --frq freq --N 53091.03
./ldsc.py --h2 meta_new2.sumstats.gz --ref-ld-chr eur_w_ld_chr/ --w-ld-chr eur_w_ld_chr/ --out meta_new2
less meta_new2.log













######meta result
meta_result <- standard_result %>% mutate(
  or = exp(Beta.meta),
  sample_size = 1/(var.meta*2*EAFcontrols.Onco*(1-EAFcontrols.Onco))
) %>% 
  select(snp.id,chr.Onco,Position.Onco,Effect.Meta,Baseline.Meta,or,sdE.meta,p.meta,
         r2.Onco,
         EAFcontrols.Onco,
         sample_size)
colnames(meta_result) <- c("snpid",
                           "CHR",
                           "bp",
                           "A1",
                           "A2",
                           "or",
                           "se",
                           "P-value",
                           "info",
                           "freq",
                           "N")
write.table(meta_result,file="/spin1/users/zhangh24/ldsc/meta_result.txt",col.names = T,quote=F)
./munge_sumstats.py --sumstats bcac_result.txt --out meta --merge-alleles w_hm3.snplist --info-min 0.3  --signed-sumstats Z,0 --frq freq 
./ldsc.py --h2 meta.sumstats.gz --ref-ld-chr eur_w_ld_chr/ --w-ld-chr eur_w_ld_chr/ --out meta
less meta.log









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



setwd("/spin1/users/zhangh24/breast_cancer_data_analysis/")
load(paste0("./whole_genome_age/ICOG/standard_analysis/result/meta_result_shared_1p.Rdata"))
all.snp <- as.data.frame(fread("./discovery_SNP/result/ResultsMeta_GWAS_iCOGs_Onco_filter_R2_MAF.txt"))
library(dplyr)
all.snp = all.snp %>% select(var_name,Effect.Onco,Baseline.Onco)
standard_result <- meta_result_shared_1p
temp.str <- strsplit(standard_result$rs_id,":")
n <- nrow(standard_result)
snpid <- rep("c",n)
for(i in 1:n){
  snpid[i] <- temp.str[[i]][1]
}
idx <- which(c(1:n)%in%grep("rs",snpid)!=T)
snpid[idx] <- paste0("chr",standard_result$CHR,":",
                     standard_result$position)[idx]
standard_result$snp.id <- snpid

standard_result = merge(standard_result,
                        all.snp,
                        by="var_name")


bcac_result <- standard_result %>% mutate(
  z = log.odds/sqrt(sigma),
  or = exp(log.odds),
  se = sqrt(sigma),
  sample_size = 1/(sigma*2*exp_freq_a1*(1-exp_freq_a1))) %>% 
  select(snp.id,CHR,position,Effect.Onco,Baseline.Onco,or,se,p.value,
         info,
         exp_freq_a1,sample_size)
colnames(bcac_result) <- c("snpid",
                           "CHR",
                           "bp",
                           "A2",
                           "A1",
                           "or",
                           "se",
                           "P",
                           "info",
                           "freq",
                           "N")
write.table(bcac_result,file="/spin1/users/zhangh24/ldsc/bcac_result_noc.txt",col.names = T,quote=F)

conda env create --file environment.yml
source activate ldsc
./munge_sumstats.py --sumstats bcac_result_noc.txt --out meta_noc --merge-alleles w_hm3.snplist --info-min 0.3 --frq freq
./ldsc.py --h2 meta_noc.sumstats.gz --ref-ld-chr eur_w_ld_chr/ --w-ld-chr eur_w_ld_chr/ --out meta_noc
less meta_noc.log