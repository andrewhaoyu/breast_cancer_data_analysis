args = commandArgs(trailingOnly=T)
i1 = as.numeric(args[[1]])

library(data.table)
library(snow)
load("/spin1/users/zhangh20/match.Rdata")
load("/spin1/users/zhangh20/breast_cancer/standard_gwas/ER-_control/icog/icog_result_odds_sd.Rdata")
load("/spin1/users/zhangh20/breast_cancer/standard_gwas/ER-_control/onco/onco_result_odds_sd.Rdata")

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

groups_id = as.numeric(cut(c(1:n),1000))

idx.group = which(groups_id == i1)
icog_result_shared = icog_result_shared[idx.group,]
onco_result_shared = onco_result_shared[idx.group,]

rs_id = icog_result_shared$rs_id
position = icog_result_shared$position

CHR = as.numeric(gsub(":.*","",rs_id))


logodds_sd_combind = cbind(icog_result_shared$logodds,icog_result_shared$sd,
  onco_result_shared$logodds,onco_result_shared$sd)

num_of_thread = 2

clus = makeCluster(num_of_thread)

clusterExport(clus,"meta_analysis_function")

aa = parRapply(clus,logodds_sd_combind,function(x)
{meta_analysis_function(x[1],x[2],x[3],x[4])})

n = nrow(icog_result_shared)
logodds_meta = aa[seq(1,2*n,2)]
sd_meta = aa[seq(2,2*n,2)]


ERNeg_control_meta_result = data.frame(rs_id,CHR,position,logodds_meta,sd_meta,stringsAsFactors=F)

save(ERNeg_control_meta_result,file=paste0("/data/zhangh20/breast_cancer/standard_gwas/ER-_control/meta/ER-_control_meta_result",i1,".Rdata"))
