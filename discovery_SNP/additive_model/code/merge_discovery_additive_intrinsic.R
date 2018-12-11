#-------------------------------------------------------------------
# Update Date: 12/10/2018
# Goal: merge discovery additive and intrinsic results together
# Author: Haoyu Zhang
#-------------------------------------------------------------------
#---------------------------------------#---------------------------------------
additive <- read.csv("./discovery_SNP/additive_model/result/additive_model_final.csv")
intrinsic <- read.csv("./discovery_SNP/additive_model/result/intrinsic_subtype_final.csv")
colnames(additive)[1] <- "rs_id"
library(dplyr)
combind <- left_join(additive,intrinsic)
write.csv(combind,file= "./discovery_SNP/additive_model/result/discovery_snp_summary_121118.csv")
#select global heterogeneity significant with FDR <= 0.05 SNPs 
combind.2 <- combind[1:15,]
combind.2.p <- combind.2[,c(21,23,25,27,29)]
