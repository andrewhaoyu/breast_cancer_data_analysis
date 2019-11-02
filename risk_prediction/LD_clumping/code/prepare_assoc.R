#-------------------------------------------------------------------
# Update Date: 11/20/2018
# Create Date: 11/20/2018
# Goal: prepare assoc dataset for plink
# Author: Haoyu Zhang
#-------------------------------------------------------------------
setwd('/data/zhangh24/breast_cancer_data_analysis/risk_prediction/')
load("./intrinsic_subtypes_whole_genome/ICOG/result/meta_result_shared_1p.Rdata")

#load("./EB_whole_genome/result/whole_gonome.rdata")


library(dplyr)
n <- nrow(whole_genome)
assoc <- whole_genome %>% mutate(p.min = pmin(p.value,FTOP_result),
                        TEST = rep("ADD",n),
                        NMISS = rep(0,n),
                        OR = exp(score),
                        STAT = rnorm(n)) %>%
                        select(CHR,SNP.ONCO,position,
                               effect_allele,
                               TEST,NMISS,
                               OR,STAT,p.min)

colnames(assoc) <- c("CHR",
                     "SNP",
                     "BP",
                     "A1",
                     "TEST",
                     "NMISS",
                     "OR",
                     "STAT",
                     "P")
write.table(assoc,file="/data/zhangh24/BCAC/impute_plink_onco/LD_assoc",col.names = T, row.names=F, quote = F)
