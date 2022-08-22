#recreate information matrix for icogs and oncoarray
#the info matrix for icogs and oncoarrary are created separately for auto and sex chromosomes
#combine them into one file
load("/data/zhangh24/breast_cancer_data_analysis/whole_genome/ICOG/ERPRHER2GRADE_fixed_baseline/result/icog_info.Rdata")
icog_info_auto = icog_info
load("/data/zhangh24/breast_cancer_data_analysis/whole_genome/ICOG/ERPRHER2GRADE_fixed_baseline/result/icog_info_sex.Rdata")
icog_info_sex = icog_info
icog_info_sex = icog_info_sex[,c(1,2,3,6:13)]
icog_info = rbind(icog_info_auto,icog_info_sex)
save(icog_info, file = "/data/zhangh24/breast_cancer_data_analysis/whole_genome/ICOG/ERPRHER2GRADE_fixed_baseline/result/icog_info_clean.Rdata")
