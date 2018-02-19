

library(data.table)
library(snow)
load("/spin1/users/zhangh20/match.Rdata")
load("/spin1/users/zhangh20/breast_cancer/standard_gwas/case_control/icog/icog_result_odds_sd.Rdata")
load("/spin1/users/zhangh20/breast_cancer/standard_gwas/case_control/onco/onco_result_odds_sd.Rdata")


idx = which(is.na(data$SNP.ICOGS)|is.na(data$SNP.ONCO)|is.na(data$var_name))

data_c = data[-idx,]


shared_rs_id = intersect(data_c$SNP.ICOGS,icog_result$rs_id)
shared_rs_id2=intersect(data_c$SNP.ONCO,onco_result$rs_id)
idx.icog_shared = which((icog_result$rs_id%in%shared_rs_id)==T)
icog_result_shared = icog_result[idx.icog_shared,]
idx.icog_match = match(shared_rs_id,icog_result_shared$rs_id)
icog_result_shared = icog_result_shared[idx.icog_match,]

idx.onco_shared = which((onco_result$rs_id%in%shared_rs_id2)==T)
onco_result_shared = onco_result[idx.onco_shared,]
idx.onco_match = match(shared_rs_id,onco_result_shared$rs_id)
onco_result_shared = onco_result_shared[idx.onco_match,]


n = nrow(icog_result_shared)
idx.data_c_shared = which((data_c$SNP.ICOGS%in%shared_rs_id)==T)
data_c_shared = data_c[idx.data_c_shared,]
all.equal(data_c_shared$SNP.ICOGS,shared_rs_id)
CHR = as.numeric(gsub("_.*","",data_c_shared$var_name))



save(CHR,file=paste0("/data/zhangh20/breast_cancer/standard_gwas/case_control/meta/CHR.Rdata"))
