setwd('/spin1/users/zhangh24/breast_cancer_data_analysis/data')
library(data.table)
###########load CIMBA data
CIMBA <- fread('./data/brca1_bc.txt',header=T)
###########load BCAC intrinsic subtype data
load(paste0("/spin1/users/zhangh24/breast_cancer_data_analysis/whole_genome_age/ICOG/Intrinsic_subtypes/result/meta_result_shared_1p.Rdata"))
################change CIMBA MarkerName into BCAC rs_id style
#CIMBA_rs_id <- gsub("_",":",CIMBA$MarkerName)
#head(CIMBA_rs_id)
#CIMBA$rs_id <- CIMBA_rs_id
head(CIMBA)
colnames(meta_result_shared_1p)[21:45] <- paste0("var",c(1:25))
##############combine CIMBA and BCAC
combine.data <- merge(CIMBA,meta_result_shared_1p,by.x = "MarkerName",
                      by.y = "var_name")
##############save BCAC and CIMBA combine data for meta-analysis
save(combine.data,file = paste0(""))