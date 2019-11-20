#prepare the data for creating ideogram using Phenogram
setwd("/data/zhangh24/breast_cancer_data_analysis/")
snp <- read.csv("./data/210_known_discovery_snp_paper_order.csv")
cate <- rep(0,nrow(snp))
cate[1:178] <- c("Known variants")
cate[179:200] <- c("Variants detected in overall analysis")
cate[201:208] <- c("Variants detected in subtypes analysis")
cate[209:210] <- c("Variants detected in subtypes analysis")
snp <- data.frame(snp,cate,stringsAsFactors = F)
colnames(snp) <- c("snp","chr","pos","phenotype")
write.table(snp,file = paste0("./discovery_SNP/ideogram_plot/result/ideogram_data.txt"),col.names = T,quote=F,row.names = F,sep = "\t")





#label all the snps that are significant in different analysis
#find snps that are significant in overall analysis
library(data.table)
setwd("/data/zhangh24/breast_cancer_data_analysis/")
snp <- read.csv("./data/210_known_discovery_snp_paper_order.csv")
snp_known <- snp[1:178,]
snp_known$cate = rep("Known variants",178)
library(dplyr)
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
  filter(p.meta<=5E-08) %>% 
  select(Best.published.SNP,CHR,position)
snp_overall$cate = rep("Sigficant in overall analysis")
#find snps that are significant in subtypes analysis

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
  filter(subtypes.p<=5E-08) %>% 
  select(Best.published.SNP,CHR.x,position.x)
snp_subtypes$cate = rep("Sigficant in subtypes analysis")

#find snps that are significant in BCAC-TN and CIMBA BRCA1  meta analysis

cimba_result_all <- as.data.frame(fread("./data/CIMBA_BCAC_meta_analysis_083019.txt",header = T))
cimba_result_all = cimba_result_all %>% 
  mutate(chr.pos = paste0(CHR,":",position))
snp_all = left_join(snp,cimba_result_all,by="chr.pos")
#remove three multi-allelic snps that are not reported before

idx <- c(62,119,78)
snp_all = snp_all[-idx,]
colnames(snp_all)[14] <- "P"
snp_cimba = snp_all %>% 
  filter(P<=5E-08) %>% 
  select(Best.published.SNP,CHR.x,position.x)
snp_cimba$cate = rep("Sigficant in BCAC TN and CIMBA BRCA1 meta-analysis analysis")
colnames(snp_cimba) <- c("SNP","CHR",
                         "position",
                         "phenotype")
colnames(snp_subtypes) <- c("SNP","CHR",
                         "position",
                         "phenotype")
colnames(snp_overall) <- c("SNP","CHR",
                            "position",
                           "phenotype")
colnames(snp_known) <- c("SNP","CHR",
                         "position",
                         "phenotype")


snp_all <- rbind(snp_known,
                 snp_overall,
                 snp_subtypes,
                 snp_cimba)
write.table(snp_all,file = paste0("./discovery_SNP/ideogram_plot/result/ideogram_data_all_analysis.txt"),col.names = T,quote=F,row.names = F,sep = "\t")


