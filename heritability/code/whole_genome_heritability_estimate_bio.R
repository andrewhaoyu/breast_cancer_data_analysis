#Goal estimate heritability using overall estimates provided by Doug and heritability estimates using missing data algorithm

#prepare the data for ld score regression file
library(data.table)
library(dplyr)
setwd('/data/zhangh24/breast_cancer_data_analysis/')
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
  beta_meta <- (var_meta)*(beta1/var1+
                             beta2/var2)
  return(list(beta_meta,var_meta))
}
#generate meta-analysis odds ratio for icog and onco
standard_result = standard_result %>% 
  mutate(BCAC_meta_beta = Twometa(beta.iCOGs,SE.iCOGs^2,beta.Onco,SE.Onco^2)[[1]],
         BCAC_meta_var = Twometa(beta.iCOGs,SE.iCOGs^2,beta.Onco,SE.Onco^2)[[2]])



onco_result <- standard_result %>% mutate(
  or = exp(beta.Onco),
  sample_size = 1/(SE.Onco^2*2*EAFcontrols.Onco*(1-EAFcontrols.Onco)),
  MAF = ifelse(EAFcontrols.Onco>0.5,1-EAFcontrols.Onco,
               EAFcontrols.Onco)) %>%
  select(snp.id,chr.Onco,Position.Onco,Effect.Onco,Baseline.Onco,or,SE.Onco,P1df_risk_chi.Onco,
         r2.Onco,
         MAF,
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
                           "MAF",
                           "N")
write.table(onco_result,file="/data/zhangh24/ldsc/onco_result.txt",col.names = T,quote=F)
./munge_sumstats.py --sumstats onco_result.txt --out onco --merge-alleles w_hm3.snplist --info-min 0.3 --maf-min 0.01
./ldsc.py --h2 onco.sumstats.gz --ref-ld-chr eur_w_ld_chr/ --w-ld-chr eur_w_ld_chr/ --out onco
less onco.log
#estimate is 0.4227 with info 0.3 cutoff
#estimate is 0.4227 with no info cutoff
#estimate is 0.4196 with info 0.3 cutoff and maf 0.01
#write.table(meta_result,file="./heritability/result/meta_result.txt",col.names = T,quote=F)








bcac_result <- standard_result %>% mutate(
  z = BCAC_meta_beta/sqrt(BCAC_meta_var),
 sample_size = 1/(BCAC_meta_var*2*EAFcontrols.Onco*(1-EAFcontrols.Onco)),
  P = 2*pnorm(-abs(z)),
 MAF = ifelse(EAFcontrols.Onco>0.5,1-EAFcontrols.Onco,
              EAFcontrols.Onco)) %>% 
  select(snp.id,chr.Onco,Position.Onco,Effect.Onco,Baseline.Onco,z,P,
         r2.Onco,
         MAF,
         sample_size)
colnames(bcac_result) <- c("snpid",
                           "CHR",
                           "bp",
                           "A2",
                           "A1",
                           "Z",
                           "P",
                           "info",
                           "MAF",
                           "N")
write.table(bcac_result,file="/data/zhangh24/ldsc/bcac_result.txt",col.names = T,quote=F)

bcac_result = fread("/data/zhangh24/ldsc/bcac_result.txt")
write.table(bcac_result_temp,file="/data/zhangh24/ldsc/bcac_result.txt",col.names = T,quote=F)

conda env create --file environment.yml
source activate ldsc
./munge_sumstats.py --sumstats bcac_result.txt --out bcac --merge-alleles w_hm3.snplist --signed-sumstats Z,0 --info-min 0.3 --maf-min 0.01
./ldsc.py --h2 bcac.sumstats.gz --ref-ld-chr eur_w_ld_chr/ --w-ld-chr eur_w_ld_chr/ --out bcac
less bcac.log
#estimate is 0.4777 with no info cutoff with maf cutoff 0.01
#estimate is 0.4777 (se 0.0374) with cutoff 0.3 with maf cutoff 0.01
#estimate is 0.4584 with cutoff 0.9





bcac_result_all <- standard_result %>% mutate(
  z = Beta.meta/sqrt(var.meta),
  sample_size = 1/(BCAC_meta_var*2*EAFcontrols.Onco*(1-EAFcontrols.Onco)),
  P = 2*pnorm(-abs(z)),
  MAF = ifelse(EAFcontrols.Onco>0.5,1-EAFcontrols.Onco,
               EAFcontrols.Onco)) %>% 
  select(snp.id,chr.Onco,Position.Onco,Effect.Onco,Baseline.Onco,z,P,
         r2.Onco,
         MAF,
         sample_size)
colnames(bcac_result_all) <- c("snpid",
                           "CHR",
                           "bp",
                           "A2",
                           "A1",
                           "Z",
                           "P",
                           "info",
                           "MAF",
                           "N")
write.table(bcac_result_all,file="/data/zhangh24/ldsc/bcac_result_all.txt",col.names = T,quote=F)

./munge_sumstats.py --sumstats bcac_result_all.txt --out bcac_all --merge-alleles w_hm3.snplist --signed-sumstats Z,0 --info-min 0.3 --maf-min 0.01
./ldsc.py --h2 bcac_all.sumstats.gz --ref-ld-chr eur_w_ld_chr/ --w-ld-chr eur_w_ld_chr/ --out bcac_all
less bcac_all.log
#estimate is 0.5593 (0.0355) with info cutoff 0.3 with maf cutoff 0.01










#match the SNPs to hapmap3 SNPs in JHPCE
load("/data/zhangh24/ldsc/BCAC.meta.data.result.jhpce.Rdata")
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
  P = 2*pnorm(-abs(z)),
  MAF = ifelse(EAFcontrols.Onco>0.5,1-EAFcontrols.Onco,
               EAFcontrols.Onco)) %>% 
  select(snp.id,chr.Onco,Position.Onco,Effect.Onco,Baseline.Onco,z,P,
         r2.Onco,
         MAF,
         sample_size)
colnames(bcac_result_jh) <- c("snpid",
                           "CHR",
                           "bp",
                           "A2",
                           "A1",
                           "Z",
                           "P",
                           "info",
                           "MAF",
                           "N")
write.table(bcac_result_jh,file="/data/zhangh24/ldsc/bcac_result_jh.txt",col.names = T,quote=F)

conda env create --file environment.yml
source activate ldsc
./munge_sumstats.py --sumstats bcac_result_jh.txt --out bcac_jh --merge-alleles w_hm3.snplist --signed-sumstats Z,0 --info-min 0.3 --maf-min 0.01
./ldsc.py --h2 bcac_jh.sumstats.gz --ref-ld-chr eur_w_ld_chr/ --w-ld-chr eur_w_ld_chr/ --out bcac_jh
less bcac_jh.log
#estimate is 0.4592 with info cutoff as 0.9
#estimate is 0.4789 with info cutoff as 0.3 maf 0.01
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
  MAF = ifelse(bcac_onco2_eaf_controls>0.5,1-bcac_onco2_eaf_controls,bcac_onco2_eaf_controls),
  z = bcac_onco2_beta/bcac_onco2_se,
  sample_size = 1/(bcac_onco2_se^2*2*bcac_onco2_eaf_controls*(1-bcac_onco2_eaf_controls))) %>% 
  select(snp.id,chr,position_b37,a0,a1,z,bcac_onco2_P1df_Wald,
         bcac_onco2_r2,
         MAF,
         sample_size)


colnames(onco_result_pub) <- c("snpid",
                              "CHR",
                              "bp",
                              "A2",
                              "A1",
                              "Z",
                              "P",
                              "info",
                              "MAF",
                              "N")
write.table(onco_result_pub,file="/data/zhangh24/ldsc/onco_result_pub.txt",col.names = T,quote=F)
./munge_sumstats.py --sumstats onco_result_pub.txt --out onco_pub --merge-alleles w_hm3.snplist --info-min 0.3 --maf-min 0.01
./ldsc.py --h2 onco_pub.sumstats.gz --ref-ld-chr eur_w_ld_chr/ --w-ld-chr eur_w_ld_chr/ --out onco_pub
less onco_pub.log
#estimate 0.4158 with info cutoff 0.9
#estimate 0.4341 with info cutoff as 0.3
#estimate 0.4341 with no info cutoff
#estimate 0.4308 with info cutoff as 0.3 maf 0.01


#public avaiable data with only icogs and oncoarray meta-analysis
bcac_result_pub <- snp.data.onco %>% mutate(
  bcac_onco2_beta = as.numeric(bcac_onco2_beta),
  bcac_onco2_se = as.numeric(bcac_onco2_se),
  bcac_onco2_eaf_controls = as.numeric(bcac_onco2_eaf_controls),
  bcac_icogs2_se = as.numeric(bcac_icogs2_se),
  bcac_icogs2_beta = as.numeric(bcac_icogs2_beta),
  MAF = ifelse(bcac_onco2_eaf_controls>0.5,1-bcac_onco2_eaf_controls,bcac_onco2_eaf_controls),
  bcac_meta_beta = Twometa(bcac_icogs2_beta,bcac_icogs2_se^2,bcac_onco2_beta,bcac_onco2_se^2)[[1]],
  bcac_meta_var = Twometa(bcac_icogs2_beta,bcac_icogs2_se^2,bcac_onco2_beta,bcac_onco2_se^2)[[2]],
  z = bcac_meta_beta/sqrt(bcac_meta_var),
  p = 2*pnorm(-abs(z),lower.tail = T),
  sample_size = 1/(bcac_meta_var*2*bcac_onco2_eaf_controls*(1-bcac_onco2_eaf_controls))) %>% 
  select(snp.id,chr,position_b37,a0,a1,z,p,
         bcac_onco2_r2,
         MAF,
         sample_size)


colnames(bcac_result_pub) <- c("snpid",
                               "CHR",
                               "bp",
                               "A2",
                               "A1",
                               "Z",
                               "P",
                               "info",
                               "MAF",
                               "N")
write.table(bcac_result_pub,file="/data/zhangh24/ldsc/bcac_result_pub.txt",col.names = T,quote=F)
./munge_sumstats.py --sumstats bcac_result_pub.txt --out bcac_pub --merge-alleles w_hm3.snplist --info-min 0.3 --maf-min 0.01
./ldsc.py --h2 bcac_pub.sumstats.gz --ref-ld-chr eur_w_ld_chr/ --w-ld-chr eur_w_ld_chr/ --out bcac_pub
less bcac_pub.log
#estimate 0.496 (0.042) with info cutoff as 0.3 maf 0.01
#estimate 0.4835 (0.0414)with info cutoff as 0.9 maf 0.01
#estimate 0.4958 (0.0419) with info cutoff as 0.3 maf 0.001


bcac_result_pub_nplug <- snp.data.onco %>% mutate(
  bcac_onco2_beta = as.numeric(bcac_onco2_beta),
  bcac_onco2_se = as.numeric(bcac_onco2_se),
  bcac_onco2_eaf_controls = as.numeric(bcac_onco2_eaf_controls),
  bcac_icogs2_se = as.numeric(bcac_icogs2_se),
  bcac_icogs2_beta = as.numeric(bcac_icogs2_beta),
  MAF = ifelse(bcac_onco2_eaf_controls>0.5,1-bcac_onco2_eaf_controls,bcac_onco2_eaf_controls),
  bcac_meta_beta = Twometa(bcac_icogs2_beta,bcac_icogs2_se^2,bcac_onco2_beta,bcac_onco2_se^2)[[1]],
  bcac_meta_var = Twometa(bcac_icogs2_beta,bcac_icogs2_se^2,bcac_onco2_beta,bcac_onco2_se^2)[[2]],
  z = bcac_meta_beta/sqrt(bcac_meta_var),
  p = 2*pnorm(-abs(z),lower.tail = T),
  sample_size = 1/(bcac_meta_var*2*bcac_onco2_eaf_controls*(1-bcac_onco2_eaf_controls))) %>% 
  select(snp.id,chr,position_b37,a0,a1,z,p,
         bcac_onco2_r2,
         MAF)


colnames(bcac_result_pub_nplug) <- c("snpid",
                               "CHR",
                               "bp",
                               "A2",
                               "A1",
                               "Z",
                               "P",
                               "info",
                               "MAF")
write.table(bcac_result_pub_nplug,file="/data/zhangh24/ldsc/bcac_result_pub_nplug.txt",col.names = T,quote=F)

./munge_sumstats.py --sumstats bcac_result_pub_nplug.txt --out bcac_result_pub_nplug --merge-alleles w_hm3.snplist --info-min 0.3 --maf-min 0.01 --N  48487.39
./ldsc.py --h2 bcac_result_pub_nplug.sumstats.gz --ref-ld-chr eur_w_ld_chr/ --w-ld-chr eur_w_ld_chr/ --out bcac_result_pub_nplug
less bcac_result_pub_nplug.log
#icogs sample size 46785 cases, 42892 controls, effective sample 22377
#onco sample size 61282 cases, 45494 controls, effective sample 26110.39
#gwas sample size 14910 cases, 17588 controls, effective sample 8069.33
#estimate 0.4367 (0.036) with info cutoff as 0.3 maf 0.01








#public avaiable data with  icogs, oncoarray and gwas meta-analysis
bcac_result_pub_all <- snp.data.onco %>% mutate(
  bcac_onco2_eaf_controls = as.numeric(bcac_onco2_eaf_controls),
  MAF = ifelse(bcac_onco2_eaf_controls>0.5,1-bcac_onco2_eaf_controls,bcac_onco2_eaf_controls),
  bcac_onco_icogs_gwas_beta = as.numeric(bcac_onco_icogs_gwas_beta),
  bcac_onco_icogs_gwas_se = as.numeric(bcac_onco_icogs_gwas_se),
  z = bcac_onco_icogs_gwas_beta/bcac_onco_icogs_gwas_se,
  p = bcac_onco_icogs_gwas_P1df,
  sample_size = 1/(bcac_onco_icogs_gwas_se^2*2*bcac_onco2_eaf_controls*(1-bcac_onco2_eaf_controls))) %>% 
  select(snp.id,chr,position_b37,a0,a1,z,p,
         bcac_onco2_r2,
         MAF,
         sample_size)


colnames(bcac_result_pub_all) <- c("snpid",
                               "CHR",
                               "bp",
                               "A2",
                               "A1",
                               "Z",
                               "P",
                               "info",
                               "MAF",
                               "N")
write.table(bcac_result_pub_all,file="/data/zhangh24/ldsc/bcac_result_pub_all.txt",col.names = T,quote=F)
./munge_sumstats.py --sumstats bcac_result_pub_all.txt --out bcac_pub_all --merge-alleles w_hm3.snplist --info-min 0.3 --maf-min 0.01
./ldsc.py --h2 bcac_pub_all.sumstats.gz --ref-ld-chr eur_w_ld_chr/ --w-ld-chr eur_w_ld_chr/ --out bcac_pub_all
less bcac_pub_all.log
#estimate  0.5372 (0.0408) with info cutoff as 0.3 maf 0.01

bcac_result_pub_all_nplug <- snp.data.onco %>% mutate(
  bcac_onco2_eaf_controls = as.numeric(bcac_onco2_eaf_controls),
  MAF = ifelse(bcac_onco2_eaf_controls>0.5,1-bcac_onco2_eaf_controls,bcac_onco2_eaf_controls),
  bcac_onco_icogs_gwas_beta = as.numeric(bcac_onco_icogs_gwas_beta),
  bcac_onco_icogs_gwas_se = as.numeric(bcac_onco_icogs_gwas_se),
  z = bcac_onco_icogs_gwas_beta/bcac_onco_icogs_gwas_se,
  p = bcac_onco_icogs_gwas_P1df,
  sample_size = 1/(bcac_onco_icogs_gwas_se^2*2*bcac_onco2_eaf_controls*(1-bcac_onco2_eaf_controls))) %>% 
  select(snp.id,chr,position_b37,a0,a1,z,p,
         bcac_onco2_r2,
         MAF)


colnames(bcac_result_pub_all_nplug) <- c("snpid",
                                   "CHR",
                                   "bp",
                                   "A2",
                                   "A1",
                                   "Z",
                                   "P",
                                   "info",
                                   "MAF")
write.table(bcac_result_pub_all_nplug,file="/data/zhangh24/ldsc/bcac_result_pub_all_nplug.txt",col.names = T,quote=F)
./munge_sumstats.py --sumstats bcac_result_pub_all_nplug.txt --out bcac_pub_all_nplug --merge-alleles w_hm3.snplist --info-min 0.3 --maf-min 0.01 --N 56556.72
./ldsc.py --h2 bcac_pub_all_nplug.sumstats.gz --ref-ld-chr eur_w_ld_chr/ --w-ld-chr eur_w_ld_chr/ --out bcac_pub_all_nplug
less bcac_pub_all.log
#estimate  0.4657 (0.0359) with info cutoff as 0.3 maf 0.01




#load whole genome intrinsic subtypes data
load("./whole_genome_age/ICOG/Intrinsic_subtypes/result/meta_result_shared_1p_082119.Rdata")
#get the rs id
snp.id.infor <- standard_result %>% 
  select(var_name,snp.id,Effect.iCOGs,Baseline.iCOGs)
intrinsic_data <- merge(meta_result_shared_1p,snp.id.infor,by ="var_name")
colnames(intrinsic_data)[21:45] <- paste0("cov",c(1:25))
colnames(intrinsic_data)[16:20] <- paste0("beta",c(1:5))
lua <- intrinsic_data %>% mutate(
  z = beta1/sqrt(cov1),
  sample_size = 1/(cov1*2*exp_freq_a1*(1-exp_freq_a1)),
  P = 2*pnorm(-abs(z),lower.tail =T),
  MAF = ifelse(exp_freq_a1>0.5,1-exp_freq_a1,
               exp_freq_a1)) %>% 
  select(snp.id,CHR,position,Baseline.iCOGs,Effect.iCOGs,z,P,
         info,
         MAF,
         sample_size)
colnames(lua) <- c("snpid",
                   "CHR",
                   "bp",
                   "A2",
                   "A1",
                   "Z",
                   "P",
                   "info",
                   "MAF",
                   "N")
write.table(lua,file="/data/zhangh24/ldsc/lua.txt",col.names = T,quote=F)
lub <- intrinsic_data %>% mutate(
  z = beta2/sqrt(cov7),
  sample_size = 1/(cov7*2*exp_freq_a1*(1-exp_freq_a1)),
  P = 2*pnorm(-abs(z),lower.tail =T),
  MAF = ifelse(exp_freq_a1>0.5,1-exp_freq_a1,
               exp_freq_a1)) %>% 
  select(snp.id,CHR,position,Baseline.iCOGs,Effect.iCOGs,z,P,
         info,
         MAF,
         sample_size)
colnames(lub) <- c("snpid",
                   "CHR",
                   "bp",
                   "A2",
                   "A1",
                   "Z",
                   "P",
                   "info",
                   "MAF",
                   "N")
write.table(lub,file="/data/zhangh24/ldsc/lub.txt",col.names = T,quote=F)
lubher2 <- intrinsic_data %>% mutate(
  z = beta3/sqrt(cov13),
  sample_size = 1/(cov13*2*exp_freq_a1*(1-exp_freq_a1)),
  P = 2*pnorm(-abs(z),lower.tail =T),
  MAF = ifelse(exp_freq_a1>0.5,1-exp_freq_a1,
               exp_freq_a1)) %>% 
  select(snp.id,CHR,position,Baseline.iCOGs,Effect.iCOGs,z,P,
         info,
         MAF,
         sample_size)
colnames(lubher2) <- c("snpid",
                       "CHR",
                       "bp",
                       "A2",
                       "A1",
                       "Z",
                       "P",
                       "info",
                       "MAF",
                       "N")
write.table(lubher2,file="/data/zhangh24/ldsc/lubher2.txt",col.names = T,quote=F)
her2 <- intrinsic_data %>% mutate(
  z = beta4/sqrt(cov19),
  sample_size = 1/(cov19*2*exp_freq_a1*(1-exp_freq_a1)),
  P = 2*pnorm(-abs(z),lower.tail =T),
  MAF = ifelse(exp_freq_a1>0.5,1-exp_freq_a1,
               exp_freq_a1)) %>% 
  select(snp.id,CHR,position,Baseline.iCOGs,Effect.iCOGs,z,P,
         info,
         MAF,
         sample_size)
colnames(her2) <- c("snpid",
                    "CHR",
                    "bp",
                    "A2",
                    "A1",
                    "Z",
                    "P",
                    "info",
                    "MAF",
                    "N")
write.table(her2,file="/data/zhangh24/ldsc/her2.txt",col.names = T,quote=F)
tn <- intrinsic_data %>% mutate(
  z = beta5/sqrt(cov25),
  sample_size = 1/(cov25*2*exp_freq_a1*(1-exp_freq_a1)),
  P = 2*pnorm(-abs(z),lower.tail =T),
  MAF = ifelse(exp_freq_a1>0.5,1-exp_freq_a1,
               exp_freq_a1)) %>% 
  select(snp.id,CHR,position,Baseline.iCOGs,Effect.iCOGs,z,P,
         info,
         MAF,
         sample_size)
colnames(tn) <- c("snpid",
                  "CHR",
                  "bp",
                  "A2",
                  "A1",
                  "Z",
                  "P",
                  "info",
                  "MAF",
                  "N")
write.table(tn,file="/data/zhangh24/ldsc/tn.txt",col.names = T,quote=F)



# CIMBA<- as.data.frame(fread("./data/brca1_bc.txt",header=T))
# colnames(CIMBA)[5] = "beta_cimba"
# 
# temp.str <- strsplit(CIMBA$MarkerName,"_")
# n <- nrow(CIMBA)
# eff_allele = rep("c",n)
# ref_allele = rep("c",n)
# for(i in 1:n){
#   ref_allele[i] <- temp.str[[i]][3]
#   eff_allele[i] <- temp.str[[i]][4]
#   
# }
# idx <- which(eff_allele!=toupper(CIMBA$Allele2))
# CIMBA$beta_cimba[idx] = -CIMBA$beta_cimba[idx]
# CIMBA$Freq1[idx] = 1-CIMBA$Freq1[idx]
# CIMBA$ref_allele = ref_allele
# CIMBA$eff_allele = eff_allele
# CIMBA = CIMBA %>% select(var_name= MarkerName,ref_allele,eff_allele,freq_ref_allele = Freq1,beta_cimba,StdErr)
#write.table(CIMBA,file = "./data/brca1_bc_alligned_with_BCAC.txt",col.names = T,row.names = F,quote=F)
#load CIMBA data
CIMBA <- as.data.frame(fread("./data/brca1_bc_alligned_with_BCAC.txt"))
snp.id.infor <- standard_result %>% 
  select(var_name,snp.id,Effect.iCOGs,Baseline.iCOGs,r2.Onco,chr.Onco,Position.Onco)

cimba = merge(snp.id.infor,CIMBA,by="var_name")

cimba <- cimba %>% mutate(
  z = beta_cimba/StdErr,
  sample_size = 1/(StdErr^2*2*freq_ref_allele*(1-freq_ref_allele)),
  P = 2*pnorm(-abs(z),lower.tail =T),
  MAF = ifelse(freq_ref_allele>0.5,1-freq_ref_allele,
               freq_ref_allele)) %>% 
  select(snp.id,chr.Onco,Position.Onco,ref_allele,eff_allele,z,P,
         r2.Onco,
         MAF,
         sample_size)
colnames(cimba) <- c("snpid",
                  "CHR",
                  "bp",
                  "A2",
                  "A1",
                  "Z",
                  "P",
                  "info",
                  "MAF",
                  "N")
write.table(cimba,file="/data/zhangh24/ldsc/cimba.txt",col.names = T,quote=F)
./munge_sumstats.py --sumstats lua.txt --out lua --merge-alleles w_hm3.snplist --info-min 0.3 --maf-min 0.01
./munge_sumstats.py --sumstats lub.txt --out lub --merge-alleles w_hm3.snplist --info-min 0.3 --maf-min 0.01
./munge_sumstats.py --sumstats lubher2.txt --out lubher2 --merge-alleles w_hm3.snplist --info-min 0.3 --maf-min 0.01
./munge_sumstats.py --sumstats her2.txt --out her2 --merge-alleles w_hm3.snplist --info-min 0.3 --maf-min 0.01
./munge_sumstats.py --sumstats tn.txt --out tn --merge-alleles w_hm3.snplist --info-min 0.3 --maf-min 0.01
./munge_sumstats.py --sumstats cimba.txt --out cimba --merge-alleles w_hm3.snplist --info-min 0.3 --maf-min 0.01










./ldsc.py --h2 lua.sumstats.gz --ref-ld-chr eur_w_ld_chr/ --w-ld-chr eur_w_ld_chr/ --out lua
#estimate lua heritability 0.6201 (0.0564)
./ldsc.py --h2 lubher2.sumstats.gz --ref-ld-chr eur_w_ld_chr/ --w-ld-chr eur_w_ld_chr/ --out lubher2
#estimate lubher2 heritability 0.5971 (0.0768)
./ldsc.py --h2 lub.sumstats.gz --ref-ld-chr eur_w_ld_chr/ --w-ld-chr eur_w_ld_chr/ --out lub
#estimate lub heritability 0.7398 (0.0927)
./ldsc.py --h2 her2.sumstats.gz --ref-ld-chr eur_w_ld_chr/ --w-ld-chr eur_w_ld_chr/ --out her2
#estimate her2 heritability 0.6885 (0.154)
./ldsc.py --h2 tn.sumstats.gz --ref-ld-chr eur_w_ld_chr/ --w-ld-chr eur_w_ld_chr/ --out tn
#estimate tn heritability  0.4917 (0.0719)
./ldsc.py --h2 cimba.sumstats.gz --ref-ld-chr eur_w_ld_chr/ --w-ld-chr eur_w_ld_chr/ --out cimba
#estimate tn heritability 0.3094 (0.0813)


./ldsc.py --rg lua.sumstats.gz,lubher2.sumstats.gz,lub.sumstats.gz,her2.sumstats.gz,tn.sumstats.gz,cimba.sumstats.gz --ref-ld-chr eur_w_ld_chr/ --w-ld-chr eur_w_ld_chr/ --out gr_cov_lua
#genetic correlation
#lua lubher2 0.7951(0.0461)
#lua lub 0.7389 (0.0490)
#lua her2 0.5677 (0.0700)
#lua tn 0.4557(0.0521)
#lua cimba 0.3865(0.0853)
./ldsc.py --rg lubher2.sumstats.gz,lub.sumstats.gz,her2.sumstats.gz,tn.sumstats.gz,cimba.sumstats.gz  --ref-ld-chr eur_w_ld_chr/ --w-ld-chr eur_w_ld_chr/ --out gr_cov_lubher2
#genetic correlation
#lubher2 lub 0.6932(0.0722)
#lubher2 her2 0.5935 (0.1073)
#lubher2 tn 0.4041 (0.0750)
#lubher2 cimba 0.3109(0.1236)

./ldsc.py --rg lub.sumstats.gz,her2.sumstats.gz --ref-ld-chr eur_w_ld_chr/ --w-ld-chr eur_w_ld_chr/ --out lub_her2:q
:q



./ldsc.py --rg lub.sumstats.gz,her2.sumstats.gz,tn.sumstats.gz,cimba.sumstats.gz  --ref-ld-chr eur_w_ld_chr/ --w-ld-chr eur_w_ld_chr/ --out gr_cov_lub
#genetic correlation
#lub her2 0.3499 (0.1043)
#lub tn 0.6027 (0.0840)
#lub cimba 0.3814(0.1567)

./ldsc.py --rg her2.sumstats.gz,tn.sumstats.gz,cimba.sumstats.gz  --ref-ld-chr eur_w_ld_chr/ --w-ld-chr eur_w_ld_chr/ --out gr_cov_her2
#genetic correlation
#her2 tn 0.5572 (0.1256)
#her2 cimba 0.8049(0.2400)

./ldsc.py --rg tn.sumstats.gz,cimba.sumstats.gz  --ref-ld-chr eur_w_ld_chr/ --w-ld-chr eur_w_ld_chr/ --out gr_cov_tn
#genetic correlation
#tn cimba 0.8401(0.154)









./ldsc.py --h2 lua_c.sumstats.gz --ref-ld-chr eur_w_ld_chr/ --w-ld-chr eur_w_ld_chr/ --out lua_c
less lua_c.log
#estimate is 0.592 (0.063)

write.table(intrinsic_data,file="/data/zhangh24/ldsc/lua.txt",col.names = T,quote=F)
./munge_sumstats.py --sumstats lua.txt --out lua --merge-alleles w_hm3.snplist --info-min 0.3
./ldsc.py --h2 lua.sumstats.gz --ref-ld-chr eur_w_ld_chr/ --w-ld-chr eur_w_ld_chr/ --out lua
less lua.log
#

#load hapmap3 intrinsic subtypes data used by guanghao
load("/data/zhangh24/ldsc/BCAC_CIMBABRCA1_082119.Rdata")

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
write.table(intrinsic_data_lua,file="/data/zhangh24/ldsc/intrinsic_data_lua.txt",col.names = T,quote=F)
./munge_sumstats.py --sumstats intrinsic_data_lua.txt --out intrinsic_data_lua --merge-alleles w_hm3.snplist --info-min 0.3
./ldsc.py --h2 intrinsic_data_lua.sumstats.gz --ref-ld-chr eur_w_ld_chr/ --w-ld-chr eur_w_ld_chr/ --out intrinsic_data_lua
less intrinsic_data_lua.log
#estimate 0.7463 with no info cutoff
#estimate 0.7463 with info cutoff as 0.3
#estimate 0.6973 without info cutoff as 0.9

