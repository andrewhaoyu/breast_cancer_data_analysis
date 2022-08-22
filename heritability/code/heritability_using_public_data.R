#Goal estimate heritability using overall estimates provided by Doug and heritability estimates using missing data algorithm
#two data meta-analysis
Twometa <- function(beta1,var1,beta2,var2){
  var_meta <- 1/(1/var1+1/var2)
  beta_meta <- (var_meta)*(beta1/var1+
                             beta2/var2)
  return(list(beta_meta,var_meta))
}
#load the current gwas summary stat
standard_result <- as.data.frame(fread("./discovery_SNP/result/ResultsMeta_GWAS_iCOGs_Onco_filter_R2_MAF.txt"))
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
snp.id.infor <- standard_result %>% 
  select(var_name,snp.id,Effect.iCOGs,Baseline.iCOGs)

#load public avaiable data from publichsed nature 2017
library(data.table)
library(dplyr)
snp.data.onco <- as.data.frame(fread("./data/oncoarray_bcac_public_release_oct17.txt"))
snp.data.onco <- merge(snp.data.onco,snp.id.infor,by ="var_name")
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
