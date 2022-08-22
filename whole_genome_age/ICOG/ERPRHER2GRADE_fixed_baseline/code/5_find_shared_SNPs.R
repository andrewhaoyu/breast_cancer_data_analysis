#Goal: prepare the icog and onco array data for meta-analaysis
library(data.table)
library(dplyr)
#load icogs data
load("/data/zhangh24/breast_cancer_data_analysis/whole_genome_age/ICOG/ERPRHER2GRADE_fixed_baseline/result/Icog_result.Rdata")
#load onco array data
load("/data/zhangh24/breast_cancer_data_analysis/whole_genome_age/ONCO/ERPRHER2GRADE_fixed_baseline/result/onco_result.Rdata")
#load the index file to match icogs and oncoarray data
data = as.data.frame(fread("/data/zhangh24/breast_cancer_data_analysis/whole_genome_age/match"))
#first restrict to SNP that shared by icogs and oncoarray
idx = which(is.na(data$SNP.ICOGS)|is.na(data$SNP.ONCO)|is.na(data$var_name))
data_c = data[-idx,]

#combine icog result with SNP index in data_c
icog_result_shared = inner_join(data_c, icog_result, by = c("SNP.ICOGS" = "rs_id"))

#combine onco result with SNP index in data_c
onco_result_shared = inner_join(data_c, onco_result, by = c("SNP.ONCO" = "rs_id"))

#test whether icog and onco ID is the same
all.equal(icog_result_shared$var_name,onco_result_shared$var_name)
# 
# save(icog_result_shared,file="/data/zhangh24/breast_cancer_data_analysis/whole_genome_age/ICOG/ERPRHER2GRADE_fixed_baseline/result/icog_result_shared.Rdata")
# save(onco_result_shared,file="/data/zhangh24/breast_cancer_data_analysis/whole_genome_age/ONCO/ERPRHER2GRADE_fixed_baseline/result/onco_result_shared.Rdata")

#restrict to SNPs with allele frequency more than 0.01
idx.filter <- which(icog_result_shared$exp_freq_a1>=0.01&
                      onco_result_shared$exp_freq_a1>=0.01&
                      icog_result_shared$exp_freq_a1<=0.99&
                      onco_result_shared$exp_freq_a1<=0.99)
icog_result_shared_1p <- icog_result_shared[idx.filter,]
onco_result_shared_1p <- onco_result_shared[idx.filter,]


save(icog_result_shared_1p,file="/data/zhangh24/breast_cancer_data_analysis/whole_genome_age/ICOG/ERPRHER2GRADE_fixed_baseline/result/icog_result_shared_1p.Rdata")
save(onco_result_shared_1p,file="/data/zhangh24/breast_cancer_data_analysis/whole_genome_age/ONCO/ERPRHER2GRADE_fixed_baseline/result/onco_result_shared_1p.Rdata")

#find icog only SNPs
idx.icog.only = which((!is.na(data$SNP.ICOGS))&(!is.na(data$var_name))&is.na(data$SNP.ONCO))

data_icog_only <- data[idx.icog.only,]

icog_result_only_shared = inner_join(data_icog_only, icog_result, by = c("SNP.ICOGS" = "rs_id"))

#restrict to SNPs with allele frequency more than 0.01
idx.filter.icog.only <- which(icog_result_only_shared$exp_freq_a1>=0.01&
                                icog_result_only_shared$exp_freq_a1<=0.99)
icog_result_only_shared_1p <- icog_result_only_shared[idx.filter.icog.only,]

save(icog_result_only_shared_1p,file="/data/zhangh24/breast_cancer_data_analysis/whole_genome_age/ICOG/ERPRHER2GRADE_fixed_baseline/result/icog_result_only_shared_1p.Rdata")



#find onco only SNPs
idx.onco.only = which((!is.na(data$SNP.ONCO))&(!is.na(data$var_name))&is.na(data$SNP.ICOGS))

data_onco_only <- data[idx.onco.only,]

onco_result_only_shared = inner_join(data_onco_only, onco_result, by = c("SNP.ONCO" = "rs_id"))


#restrict to SNPs with allele frequency more than 0.01
idx.filter.onco.only <- which(onco_result_only_shared$exp_freq_a1>=0.01&
                      onco_result_only_shared$exp_freq_a1<=0.99)

onco_result_only_shared_1p <- onco_result_only_shared[idx.filter.onco.only,]

save(onco_result_only_shared_1p,file="/data/zhangh24/breast_cancer_data_analysis/whole_genome_age/ONCO/ERPRHER2GRADE_fixed_baseline/result/onco_result_only_shared_1p.Rdata")




#combine the result for icog and oncoarray
load("/data/zhangh24/breast_cancer_data_analysis/whole_genome_age/ICOG/ERPRHER2GRADE_fixed_baseline/result/icog_result_shared_1p.Rdata")
load("/data/zhangh24/breast_cancer_data_analysis/whole_genome_age/ONCO/ERPRHER2GRADE_fixed_baseline/result/onco_result_shared_1p.Rdata")
load("/data/zhangh24/breast_cancer_data_analysis/whole_genome_age/ICOG/ERPRHER2GRADE_fixed_baseline/result/icog_result_only_shared_1p.Rdata")
load("/data/zhangh24/breast_cancer_data_analysis/whole_genome_age/ONCO/ERPRHER2GRADE_fixed_baseline/result/onco_result_only_shared_1p.Rdata")

#specify the second.num given the model
#additive model second. num is 5
second.num <- 5

#takes the score and information matrix for meta-analysis
#the first 12 cols are basic SNPs information
#col start.ind:(start.ind+second.num-1) is the score
#col (start.ind+second.num):(start.ind+second.num^2-1) is the vector format of information matrix
start.ind = 13
icog_score_infor <- icog_result_shared_1p[,start.ind:(start.ind+second.num+second.num^2-1)]
onco_score_infor <- onco_result_shared_1p[,start.ind:(start.ind+second.num+second.num^2-1)]

icog_only_score_infor_temp <- icog_result_only_shared_1p[,start.ind:(start.ind+second.num+second.num^2-1)]
#for icog only SNPs, i created this icog_only_score_infor_sudo vector with score as 0 just for simplicity of the code
#it won't affect the result
icog_only_score_infor_sudo <- icog_only_score_infor_temp
icog_only_score_infor_sudo[] <- 0
#similar things for onco only SNPs
onco_only_score_infor_temp <- onco_result_only_shared_1p[,start.ind:(start.ind+second.num+second.num^2-1)]
onco_only_score_infor_sudo <- onco_only_score_infor_temp
onco_only_score_infor_sudo[] <- 0

#combine the icog_score_infor and onco_score_infor
#this is a long vector
#1:(second.num) is the score for icogs
#(second.num+1):(second.num+second.num^2) is the vector of information matrix for icogs
#(second.num+second.num^2+1):(second.num+second.num^2+second.num) is the score for oncoarray
#(second.num+second.num^2+second.num+1):(second.num+second.num^2+second.num+second.num^2) is the vector of information matrix for onco
shared_score_infor <- cbind(icog_score_infor,onco_score_infor)
icog_only_score_infor <- cbind(icog_only_score_infor_temp,icog_only_score_infor_sudo)
onco_only_score_infor <- cbind(onco_only_score_infor_temp,onco_only_score_infor_sudo)

#combine the score_infor for icog, onco shared SNP, icog only SNPs, onco only SNPs
icog_onco_score_infor <- rbind(shared_score_infor,icog_only_score_infor,onco_only_score_infor)
save(icog_onco_score_infor, file = "/data/zhangh24/breast_cancer_data_analysis/whole_genome_age/ICOG/ERPRHER2GRADE_fixed_baseline/result/icog_onco_score_infor.rdata")







#
# load("/data/zhangh24/breast_cancer_data_analysis/whole_genome_age/ICOG/Intrinsic_subtypes/result/meta_result_shared_1p_final.Rdata")
# 
# temp = meta_result_shared_1p[1028143,]
# scoreinfor_one =  meta_result_shared_1p[1028143,c(14:43)]
# score
# 
# 
# 







############baseline only





# load("/data/zhangh24/breast_cancer_data_analysis/whole_genome/ICOG/ERPRHER2GRADE_fixed_baseline/result/Icog_result_baseline.Rdata")
# load("/data/zhangh24/breast_cancer_data_analysis/whole_genome/ONCO/ERPRHER2GRADE_fixed_baseline/result/onco_result_baseline.Rdata")
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
# 
# 
# load("/data/zhangh24/match.Rdata")
# icog_result <- icog_result_baseline
# onco_result <- onco_result_baseline
# rm(icog_result_baseline)
# rm(onco_result_baseline)
# gc()
# 
# 
# 
# 
# idx = which(is.na(data$SNP.ICOGS)|is.na(data$SNP.ONCO)|is.na(data$var_name))
# 
# data_c = data[-idx,]
# 
# shared_rs_id = intersect(data_c$SNP.ICOGS,icog_result$rs_id)
# shared_rs_id2=intersect(data_c$SNP.ONCO,onco_result$rs_id)
# idx.icog_shared = which((icog_result$rs_id%in%shared_rs_id)==T)
# icog_result_shared = icog_result[idx.icog_shared,]
# idx.icog_match = match(shared_rs_id,icog_result_shared$rs_id)
# icog_result_shared = icog_result_shared[idx.icog_match,]
# 
# 
# 
# idx.onco_shared = which((onco_result$rs_id%in%shared_rs_id2)==T)
# onco_result_shared = onco_result[idx.onco_shared,]
# idx.onco_match = match(shared_rs_id2,onco_result_shared$rs_id)
# onco_result_shared = onco_result_shared[idx.onco_match,]
# 
# ####take out data_c
# idx.shared_data_c <- which((data_c$SNP.ICOGS%in%shared_rs_id)==T)
# 
# data_c_shared <- data_c[idx.shared_data_c,]
# idx.icog_match_data_c <- match(shared_rs_id,data_c_shared$SNP.ICOGS)
# data_c_shared <- data_c_shared[idx.icog_match_data_c,]
# 
# 
# #icog_result_shared <- icog_result_shared[,-ncol(icog_result_shared)]
# icog_result_shared <- cbind(icog_result_shared,data_c_shared)
# all.equal(icog_result_shared$rs_id,icog_result_shared$SNP.ICOGS)
# 
# 
# #onco_result_shared <- onco_result_shared[,-ncol(onco_result_shared)]
# onco_result_shared <- cbind(onco_result_shared,data_c_shared)
# all.equal(onco_result_shared$rs_id,onco_result_shared$SNP.ONCO)
# 
# save(icog_result_shared,file="/data/zhangh24/breast_cancer_data_analysis/whole_genome/ICOG/ERPRHER2GRADE_fixed_baseline/result/icog_result_shared_baseline.Rdata")
# save(onco_result_shared,file="/data/zhangh24/breast_cancer_data_analysis/whole_genome/ONCO/ERPRHER2GRADE_fixed_baseline/result/onco_result_shared_baseline.Rdata")
# 
# #load("/data/zhangh24/breast_cancer_data_analysis/whole_genome/ICOG/ERPRHER2_fixed/result/icog_result_shared.Rdata")
# #load("/data/zhangh24/breast_cancer_data_analysis/whole_genome/ONCO/ERPRHER2_fixed/result/onco_result_shared.Rdata")
# 
# 
# 
# idx.filter <- which(icog_result_shared$exp_freq_a1>=0.01&
#                       onco_result_shared$exp_freq_a1>=0.01&
#                       icog_result_shared$exp_freq_a1<=0.99&
#                       onco_result_shared$exp_freq_a1<=0.99)
# icog_result_shared_1p <- icog_result_shared[idx.filter,]
# onco_result_shared_1p <- onco_result_shared[idx.filter,]
# 
# 
# 
# 
# 
# save(icog_result_shared_1p,file="/data/zhangh24/breast_cancer_data_analysis/whole_genome/ICOG/ERPRHER2GRADE_fixed_baseline/result/icog_result_shared_1p_baseline.Rdata")
# save(onco_result_shared_1p,file="/data/zhangh24/breast_cancer_data_analysis/whole_genome/ONCO/ERPRHER2GRADE_fixed_baseline/result/onco_result_shared_1p_baseline.Rdata")
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
# 
# 
# 



















