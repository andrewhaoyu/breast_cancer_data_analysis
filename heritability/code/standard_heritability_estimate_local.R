#Goal: estimate heritability using summary level statistics from standard analysis

setwd("/dcl01/chatterj/data/hzhang1/breast_intrinsic/")
#load BCAC data
load("./two_stage_model_results/BCAC.meta.result.Rdata")
#merge BCAC results as snp.infor
snp.infor <- cbind(BCAC.meta.result[[1]],
                   BCAC.meta.result[[2]],
                   BCAC.meta.result[[3]],
                   BCAC.meta.result[[4]])
snp.infor$chr.pos <- paste0(snp.infor$CHR.x.x,"_",snp.infor$BP.x)
library(dplyr)
library(data.table)
#load in all SNPs information
all.snp <- fread("../breast_cancer_data_analysis/discovery_SNP/result/ResultsMeta_GWAS_iCOGs_Onco_filter_R2_MAF.txt",header=T)

all.snp = all.snp %>% 
  mutate(chr.pos=paste0(chr.Onco,"_",Position.Onco))
#merge the two datasets together
snp.infor.merge <- left_join(snp.infor,all.snp,
                         by="chr.pos")
idx <- which(snp.infor.merge$Effect.Meta!=
               snp.infor.merge$alleles4)
#remove the SNPs that uses different alleles in hapmap3 and BCAC
snp.infor.merge <- snp.infor.merge[-idx,]
idx <- which(snp.infor.merge$Baseline.Meta!=
                     snp.infor.merge$alleles3)

#onco array data
Twometa <- function(beta1,var1,beta2,var2){
  var_meta <- 1/(1/var1+1/var2)
  beta_meta <- (var_meta)*(1/var1*beta1+
                               1/var2*beta2)
  return(list(beta_meta,var_meta))
}


colnames(snp.infor.merge)[c(7,17)] <- c("Luminal_A",
                                        "var_Luminal_A")

snp.infor.merge = snp.infor.merge %>% 
  mutate(BCAC_meta_beta = Twometa(beta.iCOGs,SE.iCOGs^2,beta.Onco,SE.Onco^2)[[1]],
         BCAC_meta_var = Twometa(beta.iCOGs,SE.iCOGs^2,beta.Onco,SE.Onco^2)[[2]],
         Z = BCAC_meta_beta/sqrt(BCAC_meta_var),
         P = 2*pnorm(-abs(Z)),
         sample_size = 1/(BCAC_meta_var*2*EAFcontrols.Onco*(1-EAFcontrols.Onco)),
         or_BCAC_meta = exp(BCAC_meta_beta),
         se_BCAC_meta = sqrt(BCAC_meta_var),
         Z_new = Luminal_A/sqrt(var_Luminal_A),
         P_new = 2*pnorm(-abs(Z_new)),
         sample_size_new = 1/(var_Luminal_A*2*EAFcontrols.Onco*(1-EAFcontrols.Onco))
         )









snp.infor.merge.new <- snp.infor.merge %>% 
  select(SNP.x,CHR.x.x,BP.x,Baseline.Meta,Effect.Meta,
         Z_new,
         P_new,
         r2.Onco,
         EAFcontrols.Onco,
         sample_size_new
         )
colnames(snp.infor.merge.new) <- c("snpid",
                           "CHR",
                           "bp",
                           "A1",
                           "A2",
                           "Z",
                           "P",
                           "info",
                           "freq",
                           "N")
idx <- which(is.na(snp.infor.merge.new$A1))
snp.infor.merge.new <- snp.infor.merge.new[-idx,]

idx <- which(snp.infor.merge.new$Z^2>=80)
snp.infor.merge.new <- snp.infor.merge.new[-idx,]
write.table(snp.infor.merge.new,file="/dcl01/chatterj/data/hzhang1/ldsc/luminal_A_hapmap3.txt",col.names = T,quote=F)
source activate ldsc
./munge_sumstats.py --sumstats luminal_A_hapmap3.txt --out luminal_A_hapmap3 --merge-alleles w_hm3.snplist --info-min 0.3  --signed-sumstats Z,0 --frq freq
./ldsc.py --h2 luminal_A_hapmap3.sumstats.gz --ref-ld-chr eur_w_ld_chr/ --w-ld-chr eur_w_ld_chr/ --out luminal_A_hapmap3
less luminal_A_hapmap3.log


snp.infor.merge.new <- snp.infor.merge %>% 
  select(SNP.x,CHR.x.x,BP.x,Baseline.Meta,Effect.Meta,Z,P,
         r2.Onco,
         EAFcontrols.Onco,
         sample_size
  )

min(snp.infor.merge.new$freq)
colnames(snp.infor.merge.new) <- c("snpid",
                                   "CHR",
                                   "bp",
                                   "A1",
                                   "A2",
                                   "Z",
                                   "P",
                                   "info",
                                   "freq",
                                   "N")
idx <- which(is.na(snp.infor.merge.new$A1))
snp.infor.merge.new <- snp.infor.merge.new[-idx,]


write.table(snp.infor.merge.new,file="/dcl01/chatterj/data/hzhang1/ldsc/standard_result_hapmap3.txt",col.names = T,quote=F)
source activate ldsc
./munge_sumstats.py --sumstats standard_result_hapmap3.txt --out standard_rseult_hapmap3 --merge-alleles w_hm3.snplist --info-min 0.3  --signed-sumstats Z,0 --frq freq 
./ldsc.py --h2 standard_rseult_hapmap3.sumstats.gz --ref-ld-chr eur_w_ld_chr/ --w-ld-chr eur_w_ld_chr/ --out standard_rseult_hapmap3
less standard_rseult_hapmap3.log