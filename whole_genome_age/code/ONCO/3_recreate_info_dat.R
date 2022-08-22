#recreate information matrix for icogs and oncoarray
#the info matrix for icogs and oncoarrary are created separately for auto and sex chromosomes
#combine them into one file
load("/data/zhangh24/breast_cancer_data_analysis/whole_genome/ONCO/ERPRHER2GRADE_fixed_baseline/result/onco_info.Rdata")
onco_info_auto = onco_info
# 
# setwd()
# filepath = "/data/NC_BW/icogs_onco/genotype/imputed2/onco_info_files/chr23"
# Files = dir(path = filepath, pattern = "OncoArray_chr23_euro15_phased_")
# Files_sub <- data.frame(chr=rep(1,length(Files)),p1=rep(0,length(Files)),p2=rep(0,length(Files)))
# for(i in 1:length(Files)){
#   temp <- gsub("OncoArray_chr23_euro15_phased_","",Files[i])
#   temp <- gsub(".txt_info","",temp)
#   p_temp <- strsplit(temp,"_")
#   p_temp <- unlist(p_temp)
#   p1 <- as.integer(p_temp[1])
#   p2 <- as.integer(p_temp[2])
#   Files_sub[i,] <- c(23,p1,p2)
# }
# idx <- order(Files_sub$chr,Files_sub$p1)
# 
# library(data.table)
# onco_info_list = list()
# for(i in 1:length(Files)){
#   data = as.data.frame(fread(Files[idx[i]]))
#   data$CHR = 23
#   onco_info_list[[i]] = data
# }
# onco_info = rbindlist(onco_info_list)
# save(onco_info, file = "/data/zhangh24/breast_cancer_data_analysis/whole_genome/ONCO/ERPRHER2GRADE_fixed_baseline/result/onco_info_sex.Rdata")



load("/data/zhangh24/breast_cancer_data_analysis/whole_genome/ONCO/ERPRHER2GRADE_fixed_baseline/result/onco_info_sex.Rdata")

onco_info_sex = onco_info
onco_info_sex = onco_info_sex[,c(1,2,3,6:13)]
onco_info = rbind(onco_info_auto,onco_info_sex)
save(onco_info, file = "/data/zhangh24/breast_cancer_data_analysis/whole_genome/ONCO/ERPRHER2GRADE_fixed_baseline/result/onco_info_clean.Rdata")
