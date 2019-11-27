args = commandArgs(trailingOnly = T)
i1 = as.numeric(args[[1]])
setwd('/data/zhangh24/breast_cancer_data_analysis/')
#load CIMBA data
CIMBA <- as.data.frame(fread("./data/brca1_bc_alligned_with_BCAC.txt"))
load(paste0("./whole_genome_age/ICOG/Intrinsic_subtypes/result/meta_result_shared_1p_082119.Rdata"))
colnames(meta_result_shared_1p)[16:20] <- 
  c("Luminial_A","Luminal_B",
    "Luminal_B_HER2Neg",
    "HER2_Enriched",
    "TN")

intrinsic_result = meta_result_shared_1p %>% 
  select(var_name,TN,"25") %>% 
  rename(TN_var="25")
library(dplyr)

BCAC_CIMBA <- inner_join(intrinsic_result,CIMBA,
                         by="var_name")
BCAC_CIMBA = BCAC_CIMBA %>% 
  mutate(beta_cimba_var = StdErr^2)
library(bc2, lib.loc ="/home/zhangh24/R/x86_64-pc-linux-gnu-library/3.6/")





n <- nrow(BCAC_CIMBA)
start.end <- startend(n,30,i1)
start <- start.end[1]
end <- start.end[2]
total <- end-start+1
############## Meta analysis between BCAC and CIMBA
############## CIMBA only get the triple negative results
############## Triple negative is the 5th subtypes in BCAC
meta.log.odds <- rep(0,total)
meta.sigma <- rep(0,total)
meta.p <- rep(0,total)
temp = 1
for(i in start:end){
  logodds1 = BCAC_CIMBA$TN[i]
  sigma1 = BCAC_CIMBA$TN_var[i]
  logodds2 = BCAC_CIMBA$beta_cimba[i]
  sigma2 = BCAC_CIMBA$beta_cimba_var[i]
  resul.temp <- LogoddsMetaAnalysis(logodds1,sigma1,logodds2,sigma2)
  meta.log.odds[temp] <- resul.temp[[1]]
  meta.sigma[temp] <- resul.temp[[2]]
  meta.p[temp] <- 2*pnorm(-abs(resul.temp[[1]]/sqrt(resul.temp[[2]])))
  temp = temp+1
}
result.sub <- list(meta.log.odds,
                   meta.sigma,
                   meta.p)
save(result.sub,file = paste0("./discovery_SNP/CIMBA_BCAC_meta_analysis/result/meta.result.sub",i1,".Rdata"))
