load("/spin1/users/zhangh24/breast_cancer_data_analysis/risk_prediction/FTOP_whole_genome/ICOG/result/Icog_result.Rdata")
load("/spin1/users/zhangh24/breast_cancer_data_analysis/risk_prediction/FTOP_whole_genome/ONCO/result/onco_result.Rdata")














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

####take out data_c
idx.shared_data_c <- which((data_c$SNP.ICOGS%in%shared_rs_id)==T)

data_c_shared <- data_c[idx.shared_data_c,]
idx.icog_match_data_c <- match(shared_rs_id,data_c_shared$SNP.ICOGS)
data_c_shared <- data_c_shared[idx.icog_match_data_c,]


#icog_result_shared <- icog_result_shared[,-ncol(icog_result_shared)]
icog_result_shared <- cbind(icog_result_shared,data_c_shared)
all.equal(icog_result_shared$rs_id,icog_result_shared$SNP.ICOGS)


#onco_result_shared <- onco_result_shared[,-ncol(onco_result_shared)]
onco_result_shared <- cbind(onco_result_shared,data_c_shared)
all.equal(onco_result_shared$rs_id,onco_result_shared$SNP.ONCO)
setwd('/spin1/users/zhangh24/breast_cancer_data_analysis/risk_prediction/FTOP_whole_genome/')
# save(icog_result_shared,file="./ICOG/result/icog_result_shared.Rdata")
# save(onco_result_shared,file="./ONCO/result/onco_result_shared.Rdata")

#load("/spin1/users/zhangh24/breast_cancer_data_analysis/whole_genome/ICOG/ERPRHER2_fixed/result/icog_result_shared.Rdata")
#load("/spin1/users/zhangh24/breast_cancer_data_analysis/whole_genome/ONCO/ERPRHER2_fixed/result/onco_result_shared.Rdata")



idx.filter <- which(icog_result_shared$exp_freq_a1>=0.01&
                      onco_result_shared$exp_freq_a1>=0.01&
                      icog_result_shared$exp_freq_a1<=0.99&
                      onco_result_shared$exp_freq_a1<=0.99)
icog_result_shared_1p <- icog_result_shared[idx.filter,]
onco_result_shared_1p <- onco_result_shared[idx.filter,]





save(icog_result_shared_1p,file="./ICOG/result/icog_result_shared_1p.Rdata")
save(onco_result_shared_1p,file="./ONCO/result/onco_result_shared_1p.Rdata")


























