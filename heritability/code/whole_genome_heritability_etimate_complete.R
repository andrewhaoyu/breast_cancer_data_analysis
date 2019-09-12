#Goal estimate heritability using overall estimates not adjusted by country with same samples as subtypes analysis and subtypes analysis results from complete data

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
  beta_meta <- (var_meta)*(beta1/var1+
                             beta2/var2)
  return(list(beta_meta,var_meta))
}

#heritability estimate using icog+onco without adjust for country, all used the data with only subtypes
load(paste0("./whole_genome_age/ICOG/standard_analysis/result/meta_result_shared_1p_s.Rdata"))

colnames(meta_result_shared_1p)[15:32] = c(paste0("beta_",c("overall","Luminial_A","Luminal_B",
                                                            "Luminal_B_HER2Neg",
                                                            "HER2Enriched",
                                                            "TripleNegative")),
                                           paste0("var_",c("overall","Luminial_A","Luminal_B",
                                                           "Luminal_B_HER2Neg",
                                                           "HER2Enriched",
                                                           "TripleNegative")),
                                           paste0("p_",c("overall","Luminial_A","Luminal_B",
                                                         "Luminal_B_HER2Neg",
                                                         "HER2Enriched",
                                                         "TripleNegative")))

standard_result <- merge(standard_result,meta_result_shared_1p,by="var_name")


names <- c("overall","Luminial_A","Luminal_B",
           "Luminal_B_HER2Neg",
           "HER2Enriched",
           "TripleNegative")
#overall analysis without country
bcac_result_noc <- standard_result %>% mutate(
  z = beta_overall/sqrt(var_overall),
  sample_size = 1/(var_overall*2*EAFcontrols.Onco*(1-EAFcontrols.Onco)),
  P = 2*pnorm(-abs(z)),
  MAF = ifelse(EAFcontrols.Onco>0.5,1-EAFcontrols.Onco,
               EAFcontrols.Onco)) %>% 
  select(snp.id,chr.Onco,Position.Onco,Effect.Onco,Baseline.Onco,z,P,
         r2.Onco,
         MAF,
         sample_size)
colnames(bcac_result_noc) <- c("snpid",
                               "CHR",
                               "bp",
                               "A2",
                               "A1",
                               "Z",
                               "P",
                               "info",
                               "MAF",
                               "N")
write.table(bcac_result_noc,file="/spin1/users/zhangh24/ldsc/bcac_result_nocountry.txt",col.names = T,quote=F)
./munge_sumstats.py --sumstats bcac_result_nocountry.txt --out bcac_noc --merge-alleles w_hm3.snplist --info-min 0.3 --maf-min 0.01
./ldsc.py --h2 bcac_noc.sumstats.gz --ref-ld-chr eur_w_ld_chr/ --w-ld-chr eur_w_ld_chr/ --out bcac_noc
less bcac_noc.log
#estimate is 0.5149 (0.039)

lua_c <- standard_result %>% mutate(
  z = beta_Luminial_A/sqrt(var_Luminial_A),
  sample_size = 1/(var_Luminial_A*2*EAFcontrols.Onco*(1-EAFcontrols.Onco)),
  P = p_Luminial_A,
  MAF = ifelse(EAFcontrols.Onco>0.5,1-EAFcontrols.Onco,
               EAFcontrols.Onco)) %>% 
  select(snp.id,chr.Onco,Position.Onco,Effect.Onco,Baseline.Onco,z,P,
         r2.Onco,
         MAF,
         sample_size)
colnames(lua_c) <- c("snpid",
                     "CHR",
                     "bp",
                     "A2",
                     "A1",
                     "Z",
                     "P",
                     "info",
                     "MAF",
                     "N")
write.table(lua_c,file="/spin1/users/zhangh24/ldsc/lua_c.txt",col.names = T,quote=F)
./munge_sumstats.py --sumstats lua_c.txt --out lua_c --merge-alleles w_hm3.snplist --info-min 0.3 --maf-min 0.01
./ldsc.py --h2 lua_c.sumstats.gz --ref-ld-chr eur_w_ld_chr/ --w-ld-chr eur_w_ld_chr/ --out lua_c
less lua_c.log
#estimate is 0.592 (0.063)


lub_c <- standard_result %>% mutate(
  z = beta_Luminal_B/sqrt(var_Luminal_B),
  sample_size = 1/(var_Luminal_B*2*EAFcontrols.Onco*(1-EAFcontrols.Onco)),
  P = p_Luminal_B,
  MAF = ifelse(EAFcontrols.Onco>0.5,1-EAFcontrols.Onco,
               EAFcontrols.Onco)) %>% 
  select(snp.id,chr.Onco,Position.Onco,Effect.Onco,Baseline.Onco,z,P,
         r2.Onco,
         MAF,
         sample_size)
colnames(lub_c) <- c("snpid",
                     "CHR",
                     "bp",
                     "A2",
                     "A1",
                     "Z",
                     "P",
                     "info",
                     "MAF",
                     "N")
write.table(lub_c,file="/spin1/users/zhangh24/ldsc/lub_c.txt",col.names = T,quote=F)
./munge_sumstats.py --sumstats lub_c.txt --out lub_c --merge-alleles w_hm3.snplist --info-min 0.3 --maf-min 0.01
./ldsc.py --h2 lub_c.sumstats.gz --ref-ld-chr eur_w_ld_chr/ --w-ld-chr eur_w_ld_chr/ --out lub_c
less lub_c.log
#estimate is 0.751 (0.093)

lubher2_c <- standard_result %>% mutate(
  z = beta_Luminal_B_HER2Neg/sqrt(var_Luminal_B_HER2Neg),
  sample_size = 1/(var_Luminal_B_HER2Neg*2*EAFcontrols.Onco*(1-EAFcontrols.Onco)),
  P = p_Luminal_B_HER2Neg,
  MAF = ifelse(EAFcontrols.Onco>0.5,1-EAFcontrols.Onco,
               EAFcontrols.Onco)) %>% 
  select(snp.id,chr.Onco,Position.Onco,Effect.Onco,Baseline.Onco,z,P,
         r2.Onco,
         MAF,
         sample_size)
colnames(lubher2_c) <- c("snpid",
                         "CHR",
                         "bp",
                         "A2",
                         "A1",
                         "Z",
                         "P",
                         "info",
                         "MAF",
                         "N")
write.table(lubher2_c,file="/spin1/users/zhangh24/ldsc/lubher2_c.txt",col.names = T,quote=F)
./munge_sumstats.py --sumstats lubher2_c.txt --out lubher2_c --merge-alleles w_hm3.snplist --info-min 0.3 --maf-min 0.01
./ldsc.py --h2 lubher2_c.sumstats.gz --ref-ld-chr eur_w_ld_chr/ --w-ld-chr eur_w_ld_chr/ --out lubher2_c
less lubher2_c.log
#estimate is 0.563 (0.063)


her2_c <- standard_result %>% mutate(
  z = beta_HER2Enriched/sqrt(var_HER2Enriched),
  sample_size = 1/(var_HER2Enriched*2*EAFcontrols.Onco*(1-EAFcontrols.Onco)),
  P = p_HER2Enriched,
  MAF = ifelse(EAFcontrols.Onco>0.5,1-EAFcontrols.Onco,
               EAFcontrols.Onco)) %>% 
  select(snp.id,chr.Onco,Position.Onco,Effect.Onco,Baseline.Onco,z,P,
         r2.Onco,
         MAF,
         sample_size)
colnames(her2_c) <- c("snpid",
                      "CHR",
                      "bp",
                      "A2",
                      "A1",
                      "Z",
                      "P",
                      "info",
                      "MAF",
                      "N")
write.table(her2_c,file="/spin1/users/zhangh24/ldsc/her2_c.txt",col.names = T,quote=F)
./munge_sumstats.py --sumstats her2_c.txt --out her2_c --merge-alleles w_hm3.snplist --info-min 0.3 --maf-min 0.01
./ldsc.py --h2 her2_c.sumstats.gz --ref-ld-chr eur_w_ld_chr/ --w-ld-chr eur_w_ld_chr/ --out her2_c
less her2_c.log
#estimate is 0.612 (0.063)

TN_c <- standard_result %>% mutate(
  z = beta_TripleNegative/sqrt(var_TripleNegative),
  sample_size = 1/(var_TripleNegative*2*EAFcontrols.Onco*(1-EAFcontrols.Onco)),
  P = p_TripleNegative,
  MAF = ifelse(EAFcontrols.Onco>0.5,1-EAFcontrols.Onco,
               EAFcontrols.Onco)) %>% 
  select(snp.id,chr.Onco,Position.Onco,Effect.Onco,Baseline.Onco,z,P,
         r2.Onco,
         MAF,
         sample_size)
colnames(TN_c) <- c("snpid",
                    "CHR",
                    "bp",
                    "A2",
                    "A1",
                    "Z",
                    "P",
                    "info",
                    "MAF",
                    "N")
write.table(TN_c,file="/spin1/users/zhangh24/ldsc/TN_c.txt",col.names = T,quote=F)
./munge_sumstats.py --sumstats TN_c.txt --out TN_c --merge-alleles w_hm3.snplist --info-min 0.3 --maf-min 0.01
./ldsc.py --h2 TN_c.sumstats.gz --ref-ld-chr eur_w_ld_chr/ --w-ld-chr eur_w_ld_chr/ --out TN_c
less TN_c.log
#estimate is 0.52 (0.085)