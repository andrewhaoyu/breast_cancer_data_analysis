load("/data/zhangh24/breast_cancer_data_analysis/whole_genome_age/ICOG/ERPRHER2GRADE_fixed_baseline/result/icog_result_shared_1p.Rdata")
load("/data/zhangh24/breast_cancer_data_analysis/whole_genome_age/ONCO/ERPRHER2GRADE_fixed_baseline/result/onco_result_shared_1p.Rdata")
load("/data/zhangh24/breast_cancer_data_analysis/whole_genome_age/ICOG/ERPRHER2GRADE_fixed_baseline/result/icog_result_only_shared_1p.Rdata")
load("/data/zhangh24/breast_cancer_data_analysis/whole_genome_age/ONCO/ERPRHER2GRADE_fixed_baseline/result/onco_result_only_shared_1p.Rdata")

#takes the basic SNPs information: first 8 columns and the last column;
meta_result_shared_1p <- icog_result_shared_1p[,c(1:8,ncol(icog_result_shared_1p))]
meta_result_shared_1p_icog_only <- icog_result_only_shared_1p[,c(1:8,ncol(icog_result_shared_1p))]
meta_result_shared_1p_onco_only <- onco_result_only_shared_1p[,c(1:8,ncol(icog_result_shared_1p))]

#combine them together
meta_result_shared_1p_no_pvalue <- rbind(meta_result_shared_1p,meta_result_shared_1p_icog_only,meta_result_shared_1p_onco_only)


#save(meta_result_shared_1p_no_pvalue,file="/data/zhangh24/breast_cancer_data_analysis/whole_genome_age/ICOG/ERPRHER2GRADE_fixed_baseline/result/meta_result_shared_1p_no_pvalue.Rdata")
 
 #load the p-value results from the 1000 subjobs
 n <- nrow( meta_result_shared_1p_no_pvalue)
 
 p.value <- rep(0,n)
 total <- 0
 for(i1 in 1:1000){
   
   load(paste0("/data/zhangh24/breast_cancer_data_analysis/whole_genome_age/ICOG/ERPRHER2GRADE_fixed_baseline/result/p_value_sub",i1,".Rdata"))
   temp <- length(pvalue_sub)
   
   p.value[total+(1:temp)] <- pvalue_sub
   
   total <- total + temp
   
 }

#combine the SNPs information with the p-value
meta_result_shared_1p <- cbind(meta_result_shared_1p_no_pvalue,p.value)
#save the results
 save(meta_result_shared_1p,file=paste0("/data/zhangh24/breast_cancer_data_analysis/whole_genome_age/ICOG/ERPRHER2GRADE_fixed_baseline/result/meta_result_shared_1p.Rdata"))
 
# fine_mapping <- read.csv("/data/zhangh24/breast_cancer_data_analysis/data/fine_mapping_regions.csv",header= T)
# 
# 
# 
# 
# 
# idx_cut <- NULL
# start <- fine_mapping$start
# end <- fine_mapping$end
# CHR <- fine_mapping$V3
# 
# 
# for(i in 1:nrow(fine_mapping)){
#   print(i)
#   chr_temp <- CHR[i]
#   start_temp <- start[i]
#   end_temp <- end[i]
#   idx <- which(meta_result_shared_1p$CHR==chr_temp&meta_result_shared_1p$position>=start_temp&
#                  meta_result_shared_1p$position<=end_temp)
#   idx_cut <- c(idx_cut,idx)
# }
# ############duplicate variables won't mater
# idx_cut <- unique(idx_cut)
# meta_result_shared_1p_filter <- meta_result_shared_1p[-idx_cut,]
# 
# save(meta_result_shared_1p_filter,file="/data/zhangh24/breast_cancer_data_analysis/whole_genome_age/ICOG/ERPRHER2GRADE_fixed_baseline/result/meta_result_shared_1p_filter_1M.Rdata")
# 
# 
#  new_filter <- read.csv("/data/zhangh24/breast_cancer_data_analysis/whole_genome/ICOG/ERPRHER2_fixed/result/Filter_based_on_Montse.csv",header=T,stringsAsFactors = F)
#  new_filter <- new_filter[1:22,]
#  new_filter[,2] <- as.numeric(gsub(",","",new_filter[,2]))
# 
# idx_cut <- NULL
# 
# position.cut <- 500*10^3
# 
# for(i in 1:nrow(new_filter)){
#   print(i)
#   chr_temp <- new_filter[i,3]
#   position_temp <- new_filter[i,2]
#   position_low <- position_temp-position.cut
#   position_high <- position_temp+position.cut
#   idx <- which(meta_result_shared_1p_filter$CHR==chr_temp&meta_result_shared_1p_filter$position>position_low&
#                  meta_result_shared_1p_filter$position<position_high)
#   idx_cut <- c(idx_cut,idx)
# }
# idx_cut <- unique(idx_cut)
# meta_result_shared_1p_filter_Ju <- meta_result_shared_1p_filter[-idx_cut,]
# 
# save(meta_result_shared_1p_filter_Ju,file="/data/zhangh24/breast_cancer_data_analysis/whole_genome_age/ICOG/ERPRHER2GRADE_fixed_baseline/result/meta_result_shared_1p_filter_1M_Ju.Rdata")
# 
# 
# 
# 
# 
# # idx <- which.min(meta_result_shared_1p_filter_Ju$p.value)
# # 
# # meta_result_shared_1p_filter_Ju[idx,]
# 
# 
# # load("/data/zhangh24/breast_cancer_data_analysis/whole_genome_age/ICOG/ERPRHER2GRADE_fixed_baseline/result/meta_result_shared_1p_filter_1M_Ju.Rdata")
# # 
# meta_result_shared_1p_filter_Ju_fix <- meta_result_shared_1p_filter_Ju
# 
# idx <- which(meta_result_shared_1p_filter_Ju_fix$CHR==2&meta_result_shared_1p_filter_Ju_fix$position== 67902524)
# meta_result_shared_1p_filter_Ju_fix[idx,]
# idx <- which(meta_result_shared_1p_filter_Ju_fix$CHR==11&meta_result_shared_1p_filter_Ju_fix$position== 120233626)
# meta_result_shared_1p_filter_Ju_fix[idx,]
# idx <- which(meta_result_shared_1p_filter_Ju_fix$CHR==18&meta_result_shared_1p_filter_Ju_fix$position== 10354649)
# meta_result_shared_1p_filter_Ju_fix[idx,]
#  
# 
# load(paste0("/data/zhangh24/breast_cancer_data_analysis/whole_genome_age/ICOG/ERPRHER2GRADE_fixed_baseline/result/meta_result_shared_1p.Rdata"))
# meta_result_shared_1p <- meta_result_shared_1p[,c(2,11,3,15)]
# write.table(meta_result_shared_1p,file=paste0("/data/zhangh24/breast_cancer_data_analysis/whole_genome_age/ICOG/ERPRHER2GRADE_fixed_baseline/result/meta_result_shared_1p_fixed.txt"),row.names = F,
#             quote=F)
# 
# 
# 
# 
# 
# 
# 
# 
# 
