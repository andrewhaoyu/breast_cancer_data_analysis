# load("/data/zhangh24/breast_cancer_data_analysis/whole_genome_age/ICOG/Intrinsic_subtypes/result/Icog_result_intrinsic_subtype.Rdata")
# load("/data/zhangh24/breast_cancer_data_analysis/whole_genome_age/ONCO/intrinsic_subtypes/result/onco_result_intrinsic_subtype.Rdata")
load("/data/zhangh24/breast_cancer_data_analysis/whole_genome_age/ICOG/Intrinsic_subtypes/result/Icog_result_intrinsic_subtype_sex.Rdata")
load("/data/zhangh24/breast_cancer_data_analysis/whole_genome_age/ONCO/intrinsic_subtypes/result/onco_result_intrinsic_subtype_sex.Rdata")
icog_result <- icog_result_casecase
onco_result <- onco_result_casecase
rm(icog_result_casecase)
rm(onco_result_casecase)
gc()












load("/data/zhangh24/match.Rdata")





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

# save(icog_result_shared,file="/data/zhangh24/breast_cancer_data_analysis/whole_genome_age/ICOG/Intrinsic_subtypes/result/icog_result_shared_sex.Rdata")
# save(onco_result_shared,file="/data/zhangh24/breast_cancer_data_analysis/whole_genome_age/ONCO/intrinsic_subtypes/result/onco_result_shared_sex.Rdata")

#load("/data/zhangh24/breast_cancer_data_analysis/whole_genome/ICOG/ERPRHER2_fixed/result/icog_result_shared.Rdata")
#load("/data/zhangh24/breast_cancer_data_analysis/whole_genome/ONCO/ERPRHER2_fixed/result/onco_result_shared.Rdata")



idx.filter <- which(icog_result_shared$exp_freq_a1>=0.01&
                      onco_result_shared$exp_freq_a1>=0.01&
                      icog_result_shared$exp_freq_a1<=0.99&
                      onco_result_shared$exp_freq_a1<=0.99)
icog_result_shared_1p <- icog_result_shared[idx.filter,]
onco_result_shared_1p <- onco_result_shared[idx.filter,]





save(icog_result_shared_1p,file="/data/zhangh24/breast_cancer_data_analysis/whole_genome_age/ICOG/Intrinsic_subtypes/result/icog_result_shared_1p_sex.Rdata")
save(onco_result_shared_1p,file="/data/zhangh24/breast_cancer_data_analysis/whole_genome_age/ONCO/intrinsic_subtypes/result/onco_result_shared_1p_sex.Rdata")









# idx.icog.only = which((!is.na(data$SNP.ICOGS))&(!is.na(data$var_name))&is.na(data$SNP.ONCO))
# 
# data_icog_only <- data[idx.icog.only,]
# 
# shared_rs_id_icog_only = intersect(data_icog_only$SNP.ICOGS,icog_result$rs_id)
# 
# 
# 
# 
# idx.icog.only_shared = which((icog_result$rs_id%in%shared_rs_id_icog_only)==T)
# icog_result_only_shared = icog_result[idx.icog.only_shared,]
# idx.icog.only_match_shared = match(shared_rs_id_icog_only,icog_result_only_shared$rs_id)
# icog_result_only_shared = icog_result_only_shared[idx.icog.only_match_shared,]
# 
# 
# ####take out data_c
# idx.icog.only.shared_data <- which((data_icog_only$SNP.ICOGS%in%shared_rs_id_icog_only)==T)
# 
# data_icog_only_shared <- data_icog_only[idx.icog.only.shared_data,]
# idx.icog_only_match_data <- match(shared_rs_id_icog_only,data_icog_only_shared$SNP.ICOGS)
# data_icog_only_shared <- data_icog_only_shared[idx.icog_only_match_data,]
# 
# 
# #icog_result_shared <- icog_result_shared[,-ncol(icog_result_shared)]
# icog_result_only_shared <- cbind(icog_result_only_shared,data_icog_only_shared )
# all.equal(icog_result_only_shared$rs_id,icog_result_only_shared$SNP.ICOGS)
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# idx.onco.only = which((!is.na(data$SNP.ONCO))&(!is.na(data$var_name))&is.na(data$SNP.ICOGS))
# 
# data_onco_only <- data[idx.onco.only,]
# 
# shared_rs_id_onco_only = intersect(data_onco_only$SNP.ONCO,onco_result$rs_id)
# 
# 
# 
# 
# idx.onco.only_shared = which((onco_result$rs_id%in%shared_rs_id_onco_only)==T)
# onco_result_only_shared = onco_result[idx.onco.only_shared,]
# idx.onco.only_match_shared = match(shared_rs_id_onco_only,onco_result_only_shared$rs_id)
# onco_result_only_shared = onco_result_only_shared[idx.onco.only_match_shared,]
# 
# ####take out data_c
# idx.onco.only.shared_data <- which((data_onco_only$SNP.ONCO%in%shared_rs_id_onco_only)==T)
# 
# data_onco_only_shared <- data_onco_only[idx.onco.only.shared_data,]
# idx.onco_only_match_data <- match(shared_rs_id_onco_only,data_onco_only_shared$SNP.ONCO)
# data_onco_only_shared <- data_onco_only_shared[idx.onco_only_match_data,]
# 
# 
# #icog_result_shared <- icog_result_shared[,-ncol(icog_result_shared)]
# onco_result_only_shared <- cbind(onco_result_only_shared,data_onco_only_shared )
# all.equal(onco_result_only_shared$rs_id,onco_result_only_shared$SNP.ONCO)
# 
# 
# 
# 
# 
# 
# idx.filter.icog.only <- which(icog_result_only_shared$exp_freq_a1>=0.01&
#                                 icog_result_only_shared$exp_freq_a1<=0.99)
# icog_result_only_shared_1p <- icog_result_only_shared[idx.filter.icog.only,]
# 
# idx.filter.onco.only <- which(onco_result_only_shared$exp_freq_a1>=0.01&
#                                 onco_result_only_shared$exp_freq_a1<=0.99)
# 
# onco_result_only_shared_1p <- onco_result_only_shared[idx.filter.onco.only,]
# 
# 
# 
# 
# 
# 
# 
# save(icog_result_only_shared_1p,file="/data/zhangh24/breast_cancer_data_analysis/whole_genome_age/ICOG/Intrinsic_subtypes/result/icog_result_only_shared_1p_082119.Rdata")
# save(onco_result_only_shared_1p,file="/data/zhangh24/breast_cancer_data_analysis/whole_genome_age/ONCO/intrinsic_subtypes/result/onco_result_only_shared_1p_082119.Rdata")
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 













