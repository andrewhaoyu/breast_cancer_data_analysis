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
#estimate is 0.4227 with info 0.3 cutoff
#estimate is 0.4227 with no info cutoff
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
write.table(bcac_result,file="/spin1/users/zhangh24/ldsc/bcac_result.txt",col.names = T,quote=F)

conda env create --file environment.yml
source activate ldsc
./munge_sumstats.py --sumstats bcac_result.txt --out bcac --merge-alleles w_hm3.snplist --signed-sumstats Z,0 --frq freq --info-min 0.3  
./ldsc.py --h2 bcac.sumstats.gz --ref-ld-chr eur_w_ld_chr/ --w-ld-chr eur_w_ld_chr/ --out bcac
less bcac.log
#estimate is 0.4777 with no info cutoff
#estimate is 0.4777 with cutoff 0.3
#estimate is 0.4584 with cutoff 0.9



#match the SNPs to hapmap3 SNPs in JHPCE
load("/spin1/users/zhangh24/ldsc/BCAC.meta.data.result.jhpce.Rdata")
#merge BCAC results as snp.infor
snp.infor <- cbind(BCAC.meta.result[[1]],
                   BCAC.meta.result[[2]],
                   BCAC.meta.result[[3]],
                   BCAC.meta.result[[4]])
snp.infor$chr.pos <- paste0(snp.infor$CHR.x.x,"_",snp.infor$BP.x)
library(dplyr)
library(data.table)
#load in all SNPs information
# all.snp <- fread("../breast_cancer_data_analysis/discovery_SNP/result/ResultsMeta_GWAS_iCOGs_Onco_filter_R2_MAF.txt",header=T)
all.snp = standard_result
all.snp = all.snp %>% 
  mutate(chr.pos=paste0(chr.Onco,"_",Position.Onco))
#merge the two datasets together
snp.infor.merge <- merge(snp.infor,all.snp,
                         by="chr.pos")
bcac_result_jh <- snp.infor.merge %>% mutate(
  z = BCAC_meta_beta/sqrt(BCAC_meta_var),
  sample_size = 1/(BCAC_meta_var*2*EAFcontrols.Onco*(1-EAFcontrols.Onco)),
  P = 2*pnorm(-abs(z))) %>% 
  select(snp.id,chr.Onco,Position.Onco,Effect.Onco,Baseline.Onco,z,P,
         r2.Onco,
         EAFcontrols.Onco,
         sample_size)
colnames(bcac_result_jh) <- c("snpid",
                           "CHR",
                           "bp",
                           "A2",
                           "A1",
                           "Z",
                           "P",
                           "info",
                           "freq",
                           "N")
write.table(bcac_result_jh,file="/spin1/users/zhangh24/ldsc/bcac_result_jh.txt",col.names = T,quote=F)

conda env create --file environment.yml
source activate ldsc
./munge_sumstats.py --sumstats bcac_result_jh.txt --out bcac_jh --merge-alleles w_hm3.snplist --signed-sumstats Z,0 --frq freq --info-min 0.3 
./ldsc.py --h2 bcac_jh.sumstats.gz --ref-ld-chr eur_w_ld_chr/ --w-ld-chr eur_w_ld_chr/ --out bcac_jh
less bcac_jh.log
#estimate is 0.4592 with info cutoff as 0.9
#estimate is 0.4789 with info cutoff as 0.3
#estimate is 0.4789 with no info cutoff

#load public avaiable data from publichsed nature 2017
snp.data.onco <- as.data.frame(fread("./data/oncoarray_bcac_public_release_oct17.txt"))
colnames(snp.data.onco)
#get the rs id
snp.id.infor <- standard_result %>% 
  select(var_name,snp.id,Effect.iCOGs,Baseline.iCOGs)
snp.data.onco <- merge(snp.data.onco,snp.id.infor,by ="var_name")
onco_result_pub <- snp.data.onco %>% mutate(
  bcac_onco2_se = as.numeric(bcac_onco2_se),
  bcac_onco2_beta = as.numeric(bcac_onco2_beta),
  bcac_onco2_eaf_controls = as.numeric(bcac_onco2_eaf_controls),
  z = bcac_onco2_beta/bcac_onco2_se,
  sample_size = 1/(bcac_onco2_se^2*2*bcac_onco2_eaf_controls*(1-bcac_onco2_eaf_controls))) %>% 
  select(snp.id,chr,position_b37,a0,a1,z,bcac_onco2_P1df_Wald,
         bcac_onco2_r2,
         bcac_onco2_eaf_controls,
         sample_size)


colnames(onco_result_pub) <- c("snpid",
                              "CHR",
                              "bp",
                              "A2",
                              "A1",
                              "Z",
                              "P",
                              "info",
                              "freq",
                              "N")
write.table(onco_result_pub,file="/spin1/users/zhangh24/ldsc/onco_result_pub.txt",col.names = T,quote=F)
./munge_sumstats.py --sumstats onco_result_pub.txt --out onco_pub --merge-alleles w_hm3.snplist --info-min 0.3
./ldsc.py --h2 onco_pub.sumstats.gz --ref-ld-chr eur_w_ld_chr/ --w-ld-chr eur_w_ld_chr/ --out onco_pub
less onco_pub.log
#estimate 0.4158 with info cutoff 0.9
#estimate 0.4341 with info cutoff as 0.3
#estimate 0.4341 with no info cutoff


#use or instead of Z statistics (conclusion, the same as Z)
#use put in sample size instead of calculated
onco_result_pub <- snp.data.onco %>% mutate(
  bcac_onco2_se = as.numeric(bcac_onco2_se),
  bcac_onco2_beta = as.numeric(bcac_onco2_beta),
  bcac_onco2_or = exp(bcac_onco2_beta),
  bcac_onco2_eaf_controls = as.numeric(bcac_onco2_eaf_controls),
  z = bcac_onco2_beta/bcac_onco2_se,
  sample_size = 1/(bcac_onco2_se^2*2*bcac_onco2_eaf_controls*(1-bcac_onco2_eaf_controls))) %>% 
  select(snp.id,chr,position_b37,a0,a1,bcac_onco2_beta,bcac_onco2_P1df_Wald,
         bcac_onco2_r2)


colnames(onco_result_pub) <- c("snpid",
                               "CHR",
                               "bp",
                               "A1",
                               "A2",
                               "beta",
                               "P",
                               "info")
write.table(onco_result_pub,file="/spin1/users/zhangh24/ldsc/onco_result_pub_or.txt",col.names = T,quote=F)
./munge_sumstats.py --sumstats onco_result_pub_or.txt --out onco_pub_or --merge-alleles w_hm3.snplist --info-min 0.3
./ldsc.py --h2 onco_pub_or.sumstats.gz --ref-ld-chr eur_w_ld_chr/ --w-ld-chr eur_w_ld_chr/ --out onco_pub_or
less onco_pub_or.log
#estimate 0.3853 with info cutoff 0.9
#estimate 0.4041 with info cutoff as 0.3
#estimate 0.4041 with no info cutoff













#load whole genome intrinsic subtypes data
load("./whole_genome_age/ICOG/Intrinsic_subtypes/result/meta_result_shared_1p_082119.Rdata")
#get the rs id
snp.id.infor <- standard_result %>% 
  select(var_name,snp.id,Effect.iCOGs,Baseline.iCOGs)
intrinsic_data <- merge(meta_result_shared_1p,snp.id.infor,by ="var_name")
colnames(intrinsic_data)[21:45] <- paste0("cov",c(1:25))
colnames(intrinsic_data)[16:20] <- paste0("beta",c(1:5))
intrinsic_data <- intrinsic_data %>% mutate(
  z = beta1/sqrt(cov1),
  p = 2*pnorm(-abs(z),lower.tail = T),
  sample_size = 1/(cov1*2*exp_freq_a1*(1-exp_freq_a1))) %>% 
  select(snp.id,CHR,position,Baseline.iCOGs,Effect.iCOGs,z,p,
         info,
         exp_freq_a1,
         sample_size)


colnames(intrinsic_data) <- c("snpid",
                               "CHR",
                               "bp",
                               "A2",
                               "A1",
                               "Z",
                               "P",
                               "info",
                               "freq",
                               "N")
idx <- which(meta_result_shared_1p$CHR==1&
               meta_result_shared_1p$position==100000827)
write.table(intrinsic_data,file="/spin1/users/zhangh24/ldsc/intrinsic_data.txt",col.names = T,quote=F)
./munge_sumstats.py --sumstats intrinsic_data.txt --out intrinsic_data --merge-alleles w_hm3.snplist --info-min 0.3
./ldsc.py --h2 intrinsic_data.sumstats.gz --ref-ld-chr eur_w_ld_chr/ --w-ld-chr eur_w_ld_chr/ --out intrinsic_data
less intrinsic_data.log
#estimate is 0.7316 with info cutoff as 0.9
#estimate is 0.6202 with info cutoff as 0.3
#estimate is 0.6202 with no info cutoff

#load hapmap3 intrinsic subtypes data used by guanghao
load("/spin1/users/zhangh24/ldsc/BCAC_CIMBABRCA1_082119.Rdata")

BCAC.meta.result.new <- cbind(BCAC.meta.result.new[[1]],
                              BCAC.meta.result.new[[2]],
                              BCAC.meta.result.new[[3]],
                              BCAC.meta.result.new[[4]])
colnames(BCAC.meta.result.new)[7] <- "luminal_a"
#get the rs id
BCAC.meta.result.new$chr.pos <- paste0(BCAC.meta.result.new$CHR,"_",BCAC.meta.result.new$Position)
standard_result = standard_result %>% 
  mutate(chr.pos=paste0(chr.Onco,"_",Position.Onco))
BCAC.meta.result.new <- merge(BCAC.meta.result.new,all.snp,
                         by="chr.pos")
intrinsic_data_lua <- BCAC.meta.result.new %>% mutate(
  z = luminal_a/sqrt(cov7),
  p = 2*pnorm(-abs(z),lower.tail = T),
  sample_size = 1/(cov7*2*freq.onco*(1-freq.onco))) %>% 
  select(SNP,CHR,Position,Reference_allele,Effect_allele,z,p,r2.Onco,
         freq.onco,
         sample_size)


colnames(intrinsic_data_lua) <- c("snpid",
                              "CHR",
                              "bp",
                              "A2",
                              "A1",
                              "Z",
                              "P",
                              "info",
                              "freq",
                              "N")
write.table(intrinsic_data_lua,file="/spin1/users/zhangh24/ldsc/intrinsic_data_lua.txt",col.names = T,quote=F)
./munge_sumstats.py --sumstats intrinsic_data_lua.txt --out intrinsic_data_lua --merge-alleles w_hm3.snplist --info-min 0.3
./ldsc.py --h2 intrinsic_data_lua.sumstats.gz --ref-ld-chr eur_w_ld_chr/ --w-ld-chr eur_w_ld_chr/ --out intrinsic_data_lua
less intrinsic_data_lua.log
#estimate 0.7463 with no info cutoff
#estimate 0.7463 with info cutoff as 0.3
#estimate 0.6973 without info cutoff as 0.9





























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
./munge_sumstats.py --sumstats bcac_result_noc.txt --out meta_noc --merge-alleles w_hm3.snplist --info-min 0.3 --frq freq --info-min 0.3
./ldsc.py --h2 meta_noc.sumstats.gz --ref-ld-chr eur_w_ld_chr/ --w-ld-chr eur_w_ld_chr/ --out meta_noc
less meta_noc.log