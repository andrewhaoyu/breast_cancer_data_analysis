#goal: meta analysis using mixed-effect two-stage model
#since sex chromosome data is small, directly run on single computer
args = commandArgs(trailingOnly = T)
i1 = as.numeric(args[[1]])

load("/data/zhangh24/breast_cancer_data_analysis/whole_genome_age/ICOG/ERPRHER2GRADE_casecase/result/icog_result_shared_1p_sex.Rdata")
load("/data/zhangh24/breast_cancer_data_analysis/whole_genome_age/ONCO/ERPRHER2GRADE_casecase/result/onco_result_shared_1p_sex.Rdata")

second.num <- 3

icog_result_shared_1p_casecase <- icog_result_shared_1p
onco_result_shared_1p_casecase <- onco_result_shared_1p
# icog_result_only_shared_1p_casecase <- icog_result_only_shared_1p
# onco_result_only_shared_1p_casecase <- onco_result_only_shared_1p


icog_score_infor_casecase <- icog_result_shared_1p_casecase[,11:(11+second.num+second.num^2-1)]
onco_score_infor_casecase <- onco_result_shared_1p_casecase[,11:(11+second.num+second.num^2-1)]




icog_onco_score_infor_casecase <- cbind(icog_score_infor_casecase,onco_score_infor_casecase)

n <- nrow(icog_onco_score_infor_casecase)


load("/data/zhangh24/breast_cancer_data_analysis/whole_genome_age/ICOG/ERPRHER2GRADE_fixed_baseline/result/icog_result_shared_1p_sex.Rdata")
load("/data/zhangh24/breast_cancer_data_analysis/whole_genome_age/ONCO/ERPRHER2GRADE_fixed_baseline/result/onco_result_shared_1p_sex.Rdata")
library(bc2, lib.loc ="/home/zhangh24/R/x86_64-pc-linux-gnu-library/3.6/")


second.num <- 5

#####get the score and infor for baseline only

icog_score_infor <- icog_result_shared_1p[,11:(11+second.num+second.num^2-1)]
onco_score_infor <- onco_result_shared_1p[,11:(11+second.num+second.num^2-1)]
icog_score_infor <- icog_score_infor[,c(1,2,6,7,11,12)]
onco_score_infor <- onco_score_infor[,c(1,2,6,7,11,12)]


icog_onco_score_infor <- cbind(icog_score_infor,onco_score_infor)

second.num <- 2

n <- nrow(icog_onco_score_infor)






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

#n <- nrow(icog_onco_score_infor)

size <- 1000
#n <- nrow(icog_onco_score_infor)


fixed.second.num <- 2
random.second.num <- 3
#registerDoParallel(no.cores)
#pvalue <- rep(0,nrow(icog_onco_score_infor))


#pvalue.list <- foreach(i=1:size)%dopar%
#{
#print(i)
start.end <- startend(n,size,i1)
start <- start.end[1]
end <- start.end[2]
pvalue_sub <- rep(0,end-start+1)
temp = 1
for(j in start:end){
  print(j)
  icog_onco_score_infor_oneline <- icog_onco_score_infor[j,]
  icog_onco_score_infor_casecase_oneline <- icog_onco_score_infor_casecase[j,]
  pvalue_sub[temp] <- MetaMixedPfunction_temp(icog_onco_score_infor_oneline,icog_onco_score_infor_casecase_oneline,fixed.second.num,random.second.num)
  temp = temp+1
}


save(pvalue_sub,file=paste0("/data/zhangh24/breast_cancer_data_analysis/whole_genome_age/ICOG/ERPRHER2GRADE_casecase/result/p_value_sub_sex",i1,".Rdata"))






#p.value <- pvalue_sub



# save(pvalue_sub,file=paste0("/data/zhangh24/breast_cancer_data_analysis/whole_genome_age/ICOG/ERPRHER2GRADE_casecase/result/p_value_sub",i1,".Rdata"))




#save(p_value_sub,file=paste0("/data/zhangh24/breast_cancer_data_analysis/whole_genome/ICOG/ERPRHER2GRADE_casecase/result/meta_result_shared_1p.Rdata"))

# known_snps <- read.csv("/data/zhangh24/breast_cancer_data_analysis/data/known_SNPs_info.csv",header= T)
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
# save(meta_result_shared_1p_filter,file="/data/zhangh24/breast_cancer_data_analysis/whole_genome/ICOG/ERPRHER2GRADE_casecase/result/meta_result_shared_1p_filter_1M.Rdata")
# 
# 
# 
# 
# new_filter <- read.csv("/data/zhangh24/breast_cancer_data_analysis/whole_genome/ICOG/ERPRHER2_fixed/result/Filter_based_on_Montse.csv",header=T,stringsAsFactors = F)
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
# save(meta_result_shared_1p_filter_Ju,file="/data/zhangh24/breast_cancer_data_analysis/whole_genome/ICOG/ERPRHER2GRADE_casecase/result/meta_result_shared_1p_filter_1M_Ju.Rdata")
# 





















