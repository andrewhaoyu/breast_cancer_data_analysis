load("/spin1/users/zhangh24/breast_cancer_data_analysis/whole_genome/ICOG/ERPRHER2_fixed/result/Icog_result.Rdata")
load("/spin1/users/zhangh24/breast_cancer_data_analysis/whole_genome/ONCO/ERPRHER2_fixed/result/onco_result.Rdata")













load("/spin1/users/zhangh24/match.Rdata")





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
idx.onco_match = match(shared_rs_id2,onco_result_shared$rs_id)
onco_result_shared = onco_result_shared[idx.onco_match,]
save(icog_result_shared,file="/spin1/users/zhangh24/breast_cancer_data_analysis/whole_genome/ICOG/ERPRHER2_fixed/result/icog_result_shared.Rdata")
save(onco_result_shared,file="/spin1/users/zhangh24/breast_cancer_data_analysis/whole_genome/ONCO/ERPRHER2_fixed/result/onco_result_shared.Rdata")





