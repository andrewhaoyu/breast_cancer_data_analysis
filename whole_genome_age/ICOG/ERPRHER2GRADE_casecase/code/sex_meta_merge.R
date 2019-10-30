load("/data/zhangh24/breast_cancer_data_analysis/whole_genome_age/ICOG/ERPRHER2GRADE_casecase/result/icog_result_shared_1p_sex.Rdata")
load("/data/zhangh24/breast_cancer_data_analysis/whole_genome_age/ONCO/ERPRHER2GRADE_casecase/result/onco_result_shared_1p_sex.Rdata")


meta_result_shared_1p <- icog_result_shared_1p[,c(1:10,(ncol(icog_result_shared_1p)-3):ncol(icog_result_shared_1p))]


meta_result_shared_1p_no_pvalue <- meta_result_shared_1p

n <- nrow(meta_result_shared_1p_no_pvalue)

p.value <- rep(0,n)
total <- 0
for(i1 in 1:1000){
  
  load(paste0("/data/zhangh24/breast_cancer_data_analysis/whole_genome_age/ICOG/ERPRHER2GRADE_casecase/result/p_value_sub_sex",i1,".Rdata"))
  temp <- length(pvalue_sub)
  
  p.value[total+(1:temp)] <- pvalue_sub
  
  total <- total + temp
  
}


meta_result_shared_1p_sex <- cbind(meta_result_shared_1p_no_pvalue,p.value)
 #save(meta_result_shared_1p_sex,file=paste0("/data/zhangh24/breast_cancer_data_analysis/whole_genome_age/ICOG/ERPRHER2GRADE_casecase/result/meta_result_shared_1p_sex.Rdata"))




load(paste0("/data/zhangh24/breast_cancer_data_analysis/whole_genome_age/ICOG/ERPRHER2GRADE_casecase/result/meta_result_shared_1p.Rdata"))

head(meta_result_shared_1p)
#match the column name between the sex dataset and #22 auto chromosomes data
names.sex <- colnames(meta_result_shared_1p_sex)
names.standard <- colnames(meta_result_shared_1p)
#sex data 1:3,6:15
#standard data 1:8,11:15
meta_result_shared_1p_final = rbind(
  meta_result_shared_1p[,c(1:8,11:15)],
  meta_result_shared_1p_sex[,c(1:3,6:15)] )
meta_result_shared_1p <- meta_result_shared_1p_final

save(meta_result_shared_1p,file = paste0("/data/zhangh24/breast_cancer_data_analysis/whole_genome_age/ICOG/ERPRHER2GRADE_casecase/result/meta_result_shared_1p_final.Rdata"))


meta_result_shared_1p <- meta_result_shared_1p[,c(2,9,3,13)]
write.table(meta_result_shared_1p,file=paste0("/data/zhangh24/breast_cancer_data_analysis/whole_genome_age/ICOG/ERPRHER2GRADE_casecase/result/meta_result_shared_1p_mixed.txt"),row.names = F,quote = F)







# fine_mapping <- read.csv("/data/zhangh24/breast_cancer_data_analysis/data/fine_mapping_regions.csv",header= T)
# 
# idx <- which(meta_result_shared_1p$position== 47780223&meta_result_shared_1p$CHR==21)
# meta_result_shared_1p[idx,]
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
# save(meta_result_shared_1p_filter,file="/data/zhangh24/breast_cancer_data_analysis/whole_genome_age/ICOG/ERPRHER2GRADE_casecase/result/meta_result_shared_1p_filter_1M.Rdata")
# 
# 
# new_filter <- read.csv("/data/zhangh24/breast_cancer_data_analysis/whole_genome/ICOG/ERPRHER2_fixed/result/Filter_based_on_Montse.csv",header=T,stringsAsFactors = F)
# new_filter <- new_filter[1:22,]
# new_filter[,2] <- as.numeric(gsub(",","",new_filter[,2]))
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
# save(meta_result_shared_1p_filter_Ju,file="/data/zhangh24/breast_cancer_data_analysis/whole_genome_age/ICOG/ERPRHER2GRADE_casecase/result/meta_result_shared_1p_filter_1M_Ju.Rdata")
# 
# 
# 
# load(paste0("/data/zhangh24/breast_cancer_data_analysis/whole_genome_age/ICOG/ERPRHER2GRADE_casecase/result/meta_result_shared_1p.Rdata"))









# idx <- which.min(meta_result_shared_1p_filter_Ju$p.value)
# 
# meta_result_shared_1p_filter_Ju[idx,]





# load("/data/zhangh24/breast_cancer_data_analysis/whole_genome_age/ICOG/ERPRHER2GRADE_casecase/result/meta_result_shared_1p_filter_1M_Ju.Rdata")
# 
# idx <- which(meta_result_shared_1p_filter_Ju$CHR==1&meta_result_shared_1p_filter_Ju$position== 120485335)
# meta_result_shared_1p_filter_Ju[idx,]
# 
# idx <- which(meta_result_shared_1p_filter_Ju$CHR==2&meta_result_shared_1p_filter_Ju$position== 67902524)
# meta_result_shared_1p_filter_Ju[idx,]
# idx <- which(meta_result_shared_1p_filter_Ju$CHR==11&meta_result_shared_1p_filter_Ju$position== 120233626)
# meta_result_shared_1p_filter_Ju[idx,]
# idx <- which(meta_result_shared_1p_filter_Ju$CHR==18&meta_result_shared_1p_filter_Ju$position== 10354649)
# meta_result_shared_1p_filter_Ju[idx,]
# 
# 


























