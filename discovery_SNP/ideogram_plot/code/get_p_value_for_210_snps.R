#get p-value for 210 known snps
setwd("/data/zhangh24/breast_cancer_data_analysis/")
snp <- read.csv("./data/210_known_discovery_snp_paper_order.csv")
library(dplyr)
library(data.table)
snp = snp %>% 
  mutate(chr.pos = paste0(CHR,":",position))
data <- as.data.frame(fread("./discovery_SNP/result/ResultsMeta_GWAS_iCOGs_Onco_filter_R2_MAF.txt"))
data = data %>% 
  mutate(chr.pos = paste0(chr.Onco,":",Position.Onco))

snp_all = left_join(snp,data,by="chr.pos")
#remove two multi-allelic snps that are not reported before
idx <- c(62,119)
snp_all = snp_all[-idx,]
snp_overall = snp_all %>% 
  select(Best.published.SNP,CHR,position,p.meta) %>% 
  rename(SNP=Best.published.SNP,P=p.meta)


meta_result_shared_1p <- as.data.frame(fread("./whole_genome_age/ICOG/ERPRHER2GRADE_fixed_baseline/result/meta_result_shared_1p_fixed.txt",header=T))
ftop <- meta_result_shared_1p
ftop.p <- ftop$p.value

meta_result_shared_1p <- as.data.frame(fread("./whole_genome_age/ICOG/ERPRHER2GRADE_casecase/result/meta_result_shared_1p_mixed.txt",header=T))
mtop <- meta_result_shared_1p
mtop.p <- mtop$p.value
subtypes.p <- apply(cbind(ftop.p,mtop.p),1,min)

ftop$subtypes.p <- subtypes.p
ftop = ftop %>% 
  mutate(chr.pos = paste0(CHR,":",position))
snp_all = left_join(snp,ftop,by="chr.pos")
#remove two multi-allelic snps that are not reported before
idx <- c(62,119)
snp_all = snp_all[-idx,]
snp_subtypes = snp_all %>% 
  select(Best.published.SNP,CHR.x,position.x,subtypes.p) %>% 
  rename(SNP = Best.published.SNP, CHR =CHR.x,
         posiiton = position.x,P=subtypes.p )

cimba_result_all <- as.data.frame(fread("./data/CIMBA_BCAC_meta_analysis_083019.txt",header = T))
cimba_result_all = cimba_result_all %>% 
  mutate(chr.pos = paste0(CHR,":",position))
snp_all = left_join(snp,cimba_result_all,by="chr.pos")
#remove three multi-allelic snps that are not reported before

idx <- c(62,119,78)
snp_all = snp_all[-idx,]
colnames(snp_all)[14] <- "P"
snp_cimba = snp_all %>% 
  select(Best.published.SNP,CHR.x,position.x,P) %>%
  rename(SNP = Best.published.SNP, CHR =CHR.x,
         posiiton = position.x)

colnames(snp_cimba)[4] <- "TN-BRCA1 meta P"
colnames(snp_subtypes)[4] <- "Subtypes analysis P"
colnames(snp_overall)[4] <- "Overall analysis P"

snp_result <- cbind(snp_overall,snp_subtypes[,4],
                    snp_cimba[,4])
write.csv(snp_result,file = "./discovery_SNP/ideogram_plot/result/210_SNPs_result.csv")
