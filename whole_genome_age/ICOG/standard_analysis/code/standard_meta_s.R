args = commandArgs(trailingOnly = T)
i1 = as.numeric(args[[1]])
second.num <- 6
setwd("/spin1/users/zhangh24/breast_cancer_data_analysis/")
load("./whole_genome_age/ICOG/standard_analysis/result/icog_result_shared_1p_s.Rdata")
load("./whole_genome_age/ONCO/standard_analysis/result/onco_result_shared_1p_s.Rdata")
load("./whole_genome_age/ONCO/standard_analysis/result/icog_result_only_shared_1p_s.Rdata")
load("./whole_genome_age/ONCO/standard_analysis/result/onco_result_only_shared_1p_s.Rdata")






library(bc2, lib.loc ="/home/zhangh24/R/x86_64-pc-linux-gnu-library/3.6/")


#####get the score and infor for baseline only
icog_score_infor <- icog_result_shared_1p[,11:(11+second.num+second.num-1)]
onco_score_infor <- onco_result_shared_1p[,11:(11+second.num+second.num-1)]



icog_score_infor_icog_only <- icog_result_only_shared_1p[,11:(11+second.num+second.num-1)]
onco_score_infor_icog_only_sudo <- icog_score_infor_icog_only
onco_score_infor_icog_only_sudo[] <- 0

onco_score_infor_onco_only <- onco_result_only_shared_1p[,11:(11+second.num+second.num-1)]


icog_score_infor_onco_only_sudo <- onco_score_infor_onco_only
icog_score_infor_onco_only_sudo[] <- 0



icog_onco_score_infor <- cbind(icog_score_infor,onco_score_infor)
icog_onco_score_infor_icog_only <- cbind(icog_score_infor_icog_only,onco_score_infor_icog_only_sudo)
icog_onco_score_infor_onco_only <- cbind(icog_score_infor_onco_only_sudo,onco_score_infor_onco_only)
colnames(icog_onco_score_infor) <- colnames(icog_onco_score_infor_icog_only)
icog_onco_score_infor_final <- rbind(icog_onco_score_infor,icog_onco_score_infor_icog_only,icog_onco_score_infor_onco_only)
icog_onco_score_infor <- icog_onco_score_infor_final






#library(foreach)
#library(doParallel)
#no.cores <- 12
n <- nrow(icog_onco_score_infor)
#icog_onco_score_infor_temp <- icog_onco_score_infor[1:10^5,]
#n <- nrow(icog_onco_score_infor_temp)
#pvalue <- rep(0,n)
size <- 1000
#n <- nrow(icog_onco_score_infor)


#registerDoParallel(no.cores)
#pvalue <- rep(0,nrow(icog_onco_score_infor))


#pvalue.list <- foreach(i=1:size)%dopar%
#{
#print(i)
start.end <- startend(n,size,i1)
start <- start.end[1]
end <- start.end[2]
#for easy coding
second.num = 6
result_summary <- matrix(0,end-start+1,(second.num+second.num+1))
temp = 1
for(j in start:end){
  #print(j)
  second.num.s = 6
  icog_onco_score_infor_oneline_s <- icog_onco_score_infor[j,]
  icog_infor_s <- icog_onco_score_infor_oneline_s[(second.num.s+1):(second.num.s+second.num.s)]
  onco_infor_s <- icog_onco_score_infor_oneline_s[(second.num.s+second.num.s+second.num.s+1):(second.num.s+second.num.s+second.num.s+second.num.s)]
  icog_onco_score_infor_oneline <- rep(0,4)
  
  for(l in 1:6){
    #for easy coding, still use second.num
    second.num = 1
    icog_onco_score_infor_oneline[1] <- as.numeric(icog_onco_score_infor_oneline_s[l])
    icog_onco_score_infor_oneline[2] <- as.numeric(icog_infor_s[l])
    icog_onco_score_infor_oneline[3] <- as.numeric(icog_onco_score_infor_oneline_s[l+second.num.s+second.num.s])
    icog_onco_score_infor_oneline[4] <- as.numeric(onco_infor_s[second.num.s])
    
    icog_infor <- icog_onco_score_infor_oneline[2]
    onco_infor <- icog_onco_score_infor_oneline[4]
    
    if(icog_infor==0){
      print(j)
      icog_onco_score_infor_oneline[(second.num+1):(second.num+second.num)] <- as.vector(diag(10000,second.num))
    }
    if(onco_infor==0){
      print(j)
      icog_onco_score_infor_oneline[(second.num+second.num+second.num+1):(second.num+second.num+second.num+second.num)] <- as.vector(diag(10000,second.num))
    }
    
    
    result_temp <- MetaFixedPfunction_temp(icog_onco_score_infor_oneline,second.num)
    result_summary[temp,l] <- as.vector(result_temp[[1]])
    result_summary[temp,second.num.s+l] <- as.vector(result_temp[[2]])
    result_summary[temp,(second.num.s+second.num.s+1)] <- as.numeric(result_temp[[3]][3])
   
    
    
  }
  
  temp = temp+1
  
  
  
 
}

save(result_summary,file=paste0("./whole_genome_age/ICOG/ERPRHER2GRADE_casecase/result/result_sub_s",i1,".Rdata"))




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





















