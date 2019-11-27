#goal: align CIMBA data with the same order as BCAC
library(data.table)
setwd('/data/zhangh24/breast_cancer_data_analysis/')
CIMBA<- as.data.frame(fread("./data/brca1_bc.txt",header=T))
colnames(CIMBA)[5] = "beta_cimba"
#BCAC use the same order as the MarkerName for reference allele and effect allele
temp.str <- strsplit(CIMBA$MarkerName,"_")
n <- nrow(CIMBA)
eff_allele = rep("c",n)
ref_allele = rep("c",n)
for(i in 1:n){
  ref_allele[i] <- temp.str[[i]][3]
  eff_allele[i] <- temp.str[[i]][4]

}
idx <- which(eff_allele!=toupper(CIMBA$Allele2))
CIMBA$beta_cimba[idx] = -CIMBA$beta_cimba[idx]
CIMBA$Freq1[idx] = 1-CIMBA$Freq1[idx]
CIMBA$ref_allele = ref_allele
CIMBA$eff_allele = eff_allele
library(dplyr)
CIMBA = CIMBA %>% select(var_name= MarkerName,ref_allele,eff_allele,freq_ref_allele = Freq1,beta_cimba,StdErr)
write.table(CIMBA,file = "./data/brca1_bc_alligned_with_BCAC.txt",col.names = T,row.names = F,quote=F)
