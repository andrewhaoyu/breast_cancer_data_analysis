args = commandArgs(trailingOnly = T)
i1 = as.numeric(args[[1]])
load("/spin1/users/zhangh24/breast_cancer_data_analysis/whole_genome_age/ICOG/Intrinsic_subtypes/result/icog_result_shared_1p.Rdata")
load("/spin1/users/zhangh24/breast_cancer_data_analysis/whole_genome_age/ONCO/intrinsic_subtypes/result/onco_result_shared_1p.Rdata")
#load("/spin1/users/zhangh24/breast_cancer_data_analysis/whole_genome_age/ICOG/Intrinsic_subtypes/result/icog_result_shared.Rdata")
#load("/spin1/users/zhangh24/breast_cancer_data_analysis/whole_genome_age/ONCO/intrinsic_subtypes/result/onco_result_shared.Rdata")
load("/spin1/users/zhangh24/breast_cancer_data_analysis/whole_genome_age/ICOG/Intrinsic_subtypes/result/icog_result_only_shared_1p.Rdata")
load("/spin1/users/zhangh24/breast_cancer_data_analysis/whole_genome_age/ONCO/intrinsic_subtypes/result/onco_result_only_shared_1p.Rdata")

second.num <- 5

icog_result_shared_1p_casecase <- icog_result_shared_1p
onco_result_shared_1p_casecase <- onco_result_shared_1p
icog_result_only_shared_1p_casecase <- icog_result_only_shared_1p
onco_result_only_shared_1p_casecase <- onco_result_only_shared_1p


icog_score_infor_casecase <- icog_result_shared_1p_casecase[,11:(11+second.num+second.num^2-1)]
onco_score_infor_casecase <- onco_result_shared_1p_casecase[,11:(11+second.num+second.num^2-1)]

icog_score_infor_icog_only_casecase <- icog_result_only_shared_1p_casecase[,11:(11+second.num+second.num^2-1)]

onco_score_infor_icog_only_sudo_casecase <- icog_score_infor_icog_only_casecase
onco_score_infor_icog_only_sudo_casecase[] <- 0
for(i in 1:nrow(onco_score_infor_icog_only_sudo_casecase)){
  onco_score_infor_icog_only_sudo_casecase[i,(second.num+1):(second.num+second.num^2)] <- as.vector(diag(10000,second.num))
}
onco_score_infor_onco_only_casecase <- onco_result_only_shared_1p_casecase[,11:(11+second.num+second.num^2-1)]
icog_score_infor_onco_only_sudo_casecase <- onco_score_infor_onco_only_casecase
icog_score_infor_onco_only_sudo_casecase[] <- 0
for(i in 1:nrow(icog_score_infor_onco_only_sudo_casecase)){
  icog_score_infor_onco_only_sudo_casecase[i,(second.num+1):(second.num+second.num^2)] <- as.vector(diag(10000,second.num))
}


icog_onco_score_infor_casecase <- cbind(icog_score_infor_casecase,onco_score_infor_casecase)
icog_onco_score_infor_icog_only_casecase <- cbind(icog_score_infor_icog_only_casecase,onco_score_infor_icog_only_sudo_casecase)
icog_onco_score_infor_onco_only_casecase <- cbind(icog_score_infor_onco_only_sudo_casecase,onco_score_infor_onco_only_casecase)

icog_onco_score_infor_final_casecase <- rbind(icog_onco_score_infor_casecase,icog_onco_score_infor_icog_only_casecase,icog_onco_score_infor_onco_only_casecase)
icog_onco_score_infor_casecase <- icog_onco_score_infor_final_casecase

n <- nrow(icog_onco_score_infor_casecase)
# debug.one.line <- c(rep(0,second.num),as.vector(diag(second.num)),rep(0,second.num),as.vector(diag(second.num)))
# debug.idx <- c(9649548,9650051)
# for(k in 1:length(debug.idx)){
#   print(k)
#   icog_onco_score_infor_casecase[debug.idx[k],] <- debug.one.line  
# }
# 

icog_onco_score_infor <- icog_onco_score_infor_casecase



library(bc2)

rm(icog_result_shared_1p)
rm(onco_result_shared_1p)
rm(icog_result_shared_1p_casecase)
rm(onco_result_shared_1p_casecase)
gc()

#library(foreach)
#library(doParallel)
#no.cores <- 12
n <- nrow(icog_onco_score_infor)
#icog_onco_score_infor_temp <- icog_onco_score_infor[1:10^5,]
#n <- nrow(icog_onco_score_infor_temp)
#pvalue <- rep(0,n)
size <- 1000
#n <- nrow(icog_onco_score_infor)


# fixed.second.num <- 2
# random.second.num <- 3
#registerDoParallel(no.cores)
#pvalue <- rep(0,nrow(icog_onco_score_infor))


#pvalue.list <- foreach(i=1:size)%dopar%
#{
#print(i)
start.end <- startend(n,size,i1)
start <- start.end[1]
end <- start.end[2]
result_summary <- matrix(0,end-start+1,(second.num+second.num^2+1))
#pvalue_sub <- rep(0,end-start+1)
temp = 1
for(j in start:end){
  #print(j)
  icog_onco_score_infor_oneline <- icog_onco_score_infor[j,]
  result_temp <- MetaFixedPfunction_temp(icog_onco_score_infor_oneline,second.num)
  result_summary[j,1:second.num] <- as.vector(result_temp[[1]])
  result_summary[j,(second.num+1):(second.num+second.num^2)] <- as.vector(result_temp[[2]])
  result_summary[j,(second.num+second.num^2+1)] <- as.numeric(result_temp[[3]][11])
  temp = temp+1
}
#return(pvalue_sub)


#}
#stopImplicitCluster()

# pvalue <- rep(0,n)
# temp.total = 0
# for(i in 1:size){
#   print(i)
#   temp <- length(pvalue.list[[i]])
#   pvalue[(1:temp)+temp.total] <- pvalue.list[[i]]
#   temp.total <- temp+temp.total
#   
#   
# }
# 


#meta_result_shared_1p <- cbind(meta_result_shared_1p,pvalue)

save(result_summary,file=paste0("/spin1/users/zhangh24/breast_cancer_data_analysis/whole_genome_age/ICOG/Intrinsic_subtypes/result/reuslt_summary_sub",i1,".Rdata"))




#save(p_value_sub,file=paste0("/spin1/users/zhangh24/breast_cancer_data_analysis/whole_genome/ICOG/ERPRHER2GRADE_casecase/result/meta_result_shared_1p.Rdata"))

# known_snps <- read.csv("/spin1/users/zhangh24/breast_cancer_data_analysis/data/known_SNPs_info.csv",header= T)
# 
# 
# 
# 
# position.cut <- 10^6
# 
# 
# idx_cut <- NULL
# 
# 
# for(i in 1:nrow(known_snps)){
#   print(i)
#   chr_temp <- known_snps[i,3]
#   position_temp <- known_snps[i,4]
#   position_low <- position_temp-position.cut
#   position_high <- position_temp+position.cut
#   idx <- which(meta_result_shared_1p$CHR==chr_temp&meta_result_shared_1p$position>position_low&
#                  meta_result_shared_1p$position<position_high)
#   idx_cut <- c(idx_cut,idx)
# }
# ############duplicate variables won't mater
# idx_cut <- unique(idx_cut)
# meta_result_shared_1p_filter <- meta_result_shared_1p[-idx_cut,]
# 
# save(meta_result_shared_1p_filter,file="/spin1/users/zhangh24/breast_cancer_data_analysis/whole_genome/ICOG/ERPRHER2GRADE_casecase/result/meta_result_shared_1p_filter_1M.Rdata")
# 
# 
# 
# 
# new_filter <- read.csv("/spin1/users/zhangh24/breast_cancer_data_analysis/whole_genome/ICOG/ERPRHER2_fixed/result/Filter_based_on_Montse.csv",header=T,stringsAsFactors = F)
# new_filter[,2] <- as.numeric(gsub(",","",new_filter[,2]))
# 
# idx_cut <- NULL
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
# save(meta_result_shared_1p_filter_Ju,file="/spin1/users/zhangh24/breast_cancer_data_analysis/whole_genome/ICOG/ERPRHER2GRADE_casecase/result/meta_result_shared_1p_filter_1M_Ju.Rdata")
# 





















