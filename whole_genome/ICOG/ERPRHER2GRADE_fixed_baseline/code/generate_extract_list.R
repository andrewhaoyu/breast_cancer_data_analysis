load("/spin1/users/zhangh24/breast_cancer_data_analysis/whole_genome/ICOG/ERPRHER2GRADE_fixed_baseline/result/meta_result_shared_1p_filter_1M.Rdata")

idx <- which(meta_result_shared_1p_filter$pvalue<=(5*10^-6))













extract.list <- meta_result_shared_1p_filter[idx,]

save(extract.list,file = paste0("/spin1/users/zhangh24/breast_cancer_data_analysis/whole_genome/ICOG/ERPRHER2GRADE_fixed_baseline/result/extract_list.Rdata"))
