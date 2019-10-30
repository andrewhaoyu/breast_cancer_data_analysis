load("/data/zhangh24/breast_cancer_data_analysis/whole_genome_age/ICOG/ERPRHER2GRADE_fixed_baseline/result/meta_result_shared_1p_filter_1M_Ju.Rdata")
meta_result_shared_1p_filter <- meta_result_shared_1p_filter_Ju
idx <- which(meta_result_shared_1p_filter$p.value<(5*10^-6))



extract.list <- meta_result_shared_1p_filter[idx,]

# new_filter <- read.csv("/data/zhangh24/breast_cancer_data_analysis/whole_genome/ICOG/ERPRHER2_fixed/result/Filter_based_on_Montse.csv",header=T,stringsAsFactors = F)
# new_filter[,2] <- as.numeric(gsub(",","",new_filter[,2]))
# 
# 
# idx <- NULL
# 
# idx<- which((meta_result_shared_1p_filter$position%in%new_filter[,2]==T)&
#               (meta_result_shared_1p_filter$CHR%in%new_filter[,3])==T)  
# 
# extract.list <- rbind(extract.list,meta_result_shared_1p_filter[idx,])



save(extract.list,file = paste0("/data/zhangh24/breast_cancer_data_analysis/whole_genome_age/ICOG/ERPRHER2GRADE_fixed_baseline/result/extract_list.Rdata"))











