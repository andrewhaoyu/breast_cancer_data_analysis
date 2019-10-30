library(data.table)
library(snow)
load("/data/zhangh20/match.Rdata")
load("/data/zhangh20/breast_cancer/standard_gwas/ER+_control/icog/icog_result_odds_sd.Rdata")
load("/data/zhangh20/breast_cancer/standard_gwas/ER+_control/onco/onco_result_odds_sd.Rdata")

p_value_function <- function(z){
  result <- NULL
  for(i in 1:length(z)){

      result <- c(result,2*(pnorm(-abs(z[i]))))

  }
  return(result)
}

meta_analysis_function = function(logodds1,sd1,logodds2,sd2){
  sigma1 = sd1^2
  sigma2 = sd2^2
  var_meta = (1/sigma1+1/sigma2)^-1
  logodds_meta = var_meta*(logodds1/sigma1+logodds2/sigma2)
  return(c(logodds_meta,sqrt(var_meta)))
}

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
ERPos_control_meta_result_final = data.frame(rs_id,CHR,position,logodds_meta,sd_meta,stringsAsFactors=F)

total.num = 0
for(i in 1:1000){
  print(i)
  load(paste0("/data/zhangh20/breast_cancer/standard_gwas/ER+_control/meta/ER+_control_meta_result",i,".Rdata"))
  temp = nrow(ERPos_control_meta_result)
  ERPos_control_meta_result_final[(1+total.num):(temp+total.num),]=ERPos_control_meta_result
  total.num = temp+total.num
}

load("/data/zhangh20/breast_cancer/standard_gwas/case_control/meta/CHR.Rdata")

ERPos_control_meta_result_final$CHR =CHR

save(ERPos_control_meta_result_final,file=paste0("/data/zhangh20/breast_cancer/standard_gwas/ER+_control/meta/ER+_control_meta_result_final.Rdata"))
