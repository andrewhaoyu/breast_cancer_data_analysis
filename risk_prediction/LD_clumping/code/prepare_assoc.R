#-------------------------------------------------------------------
# Update Date: 11/20/2018
# Create Date: 11/20/2018
# Goal: prepare assoc dataset for plink to perform LD clumping
# Author: Haoyu Zhang
#-------------------------------------------------------------------
setwd('/data/zhangh24/breast_cancer_data_analysis/risk_prediction/')
load("./intrinsic_subtypes_whole_genome/ICOG/result/whole_gonome.rdata")

#load("./EB_whole_genome/result/whole_gonome.rdata")


library(dplyr)
n <- nrow(whole_genome)
assoc <- whole_genome %>% mutate(p.min = pmin(stan_p,FTOP_result),
                        TEST = rep("ADD",n),
                        NMISS = rep(0,n),
                        OR = exp(stan_logodds),
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

#library(data.table)
#data <- fread("/data/zhangh24/BCAC/impute_plink_onco/LD_assoc")
#idx <- which(data$SNP=="rs56305903:52041531:T:C")
#data[idx,]
#idx <- which(whole_genome$SNP.ONCO=="rs56305903:52041531:T:C")
#whole_genome[idx,]