

library(data.table)
library(snow)
load("/data/zhangh20/match.Rdata")
load("/data/zhangh20/breast_cancer/standard_gwas/case_control/icog/icog_result_odds_sd.Rdata")
load("/data/zhangh20/breast_cancer/standard_gwas/case_control/onco/onco_result_odds_sd.Rdata")


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



rs_id = rep("c",n)
position = rep(0,n)

CHR = rep(0,n)

logodds_meta = rep(0,n)
sd_meta= rep(0,n)
case_control_meta_result_final = data.frame(rs_id,CHR,position,logodds_meta,sd_meta,stringsAsFactors=F)

total.num = 0
for(i in 1:1000){
  print(i)
  load(paste0("/data/zhangh20/breast_cancer/standard_gwas/case_control/meta/case_control_meta_result",i,".Rdata"))
  temp = nrow(case_control_meta_result)
  case_control_meta_result_final[(1+total.num):(temp+total.num),]=case_control_meta_result
  total.num =total.num+temp
}
load("/data/zhangh20/breast_cancer/standard_gwas/case_control/meta/CHR.Rdata")

case_control_meta_result_final$CHR =CHR

save(case_control_meta_result_final,file=paste0("/data/zhangh20/breast_cancer/standard_gwas/case_control/meta/case_control_meta_result_final.Rdata"))
