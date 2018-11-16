#-------------------------------------------------------------------
# Update Date: 11/15/2018
# Create Date: 11/15/2018
# Goal: prepare the SNP id that is going to take out for test data
# Author: Haoyu Zhang
#-------------------------------------------------------------------
setwd('/spin1/users/zhangh24/breast_cancer_data_analysis/risk_prediction')
load("./EB_whole_genome/result/whole_gonome.rdata")
snp.onco.extract.id <- whole_genome[,14,drop=F]
write.table(snp.onco.extract.id,file = "/spin1/users/zhangh24/BCAC/impute_onco/onco_1p_shared_id.txt",quote = F,row.names=F)
