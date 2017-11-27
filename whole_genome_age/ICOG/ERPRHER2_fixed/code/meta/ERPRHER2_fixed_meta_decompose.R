load("/spin1/users/zhangh24/breast_cancer_data_analysis/whole_genome/ICOG/ERPRHER2_fixed/result/icog_result_shared_1p.Rdata")
load("/spin1/users/zhangh24/breast_cancer_data_analysis/whole_genome/ONCO/ERPRHER2_fixed/result/onco_result_shared_1p.Rdata")
library(bc2)
second.num <- 4


icog_score_infor <- icog_result_shared_1p[,11:(11+second.num+second.num^2-1)]
onco_score_infor <- onco_result_shared_1p[,11:(11+second.num+second.num^2-1)]




icog_onco_score_infor <- cbind(icog_score_infor,onco_score_infor)
meta_result_shared_1p <- icog_result_shared_1p[,c(1:10,31:34)]
rm(icog_result_shared_1p)
rm(onco_result_shared_1p)
gc()




num <- nrow(meta_result_shared_1p)
size = 1000
library(bc2)




library(foreach)
library(doParallel)
no.cores <- 20
registerDoParallel(no.cores)
foreach(i = 1:size)%dopar%
{
  print(i)
  start.end <- startend(num,size,i)  
  start <- start.end[1]
  end <- start.end[2]
  meta_result_shared_1p_sub <- meta_result_shared_1p[start:end,]
  save(meta_result_shared_1p_sub,
       file=paste0("/spin1/users/zhangh24/breast_cancer_data_analysis/whole_genome/ICOG/ERPRHER2_fixed/result/meta/meta_result_shared_1p_sub",i,".Rdata"))
  icog_onco_score_infor_sub <- icog_onco_score_infor[start:end,]
  save(icog_onco_score_infor_sub,
       file=paste0("/spin1/users/zhangh24/breast_cancer_data_analysis/whole_genome/ICOG/ERPRHER2_fixed/result/meta/icog_onco_score_infor_sub",i,".Rdata"))
 
}
stopImplicitCluster()
# for(i in 1:size){
#  
# }
# 

Filesdir <- "/spin1/users/zhangh24/breast_cancer_data_analysis/whole_genome/ICOG/ERPRHER2_fixed/result/meta/"
files <- dir(Filesdir,pattern="icog_result_shared_1p_")

files<- dir(Filesdir,pattern="icog_onco_score_infor_sub")












# #icog_onco_score_infor_temp <- icog_onco_score_infor[1:10000,]
# 
# library(foreach)
# library(doParallel)
# no.cores <- 30
# #pvalue <- rep(0,nrow(icog_onco_score_infor))
# registerDoParallel(no.cores)
# 
# pvalue <- foreach(i=1:nrow(icog_onco_score_infor),
#                   .combine=c)%dopar%
# {
#   icog_onco_score_infor_oneline <- icog_onco_score_infor[i,]
#   return(MetaPfunction(icog_onco_score_infor_oneline))
# }
# stopImplicitCluster()




# meta_result_shared_1p <- cbind(meta_result_shared_1p,pvalue)
# 
# 
# 
# 
# 
# save(meta_result_shared_1p,file=paste0("/spin1/users/zhangh24/breast_cancer_data_analysis/whole_genome/ICOG/ERPRHER2_fixed/result/meta_result_shared_1p.Rdata"))
# 
# known_snps <- read.csv("/spin1/users/zhangh24/breast_cancer_data_analysis/data/known_SNPs_info.csv",header= T)
# 
# 
# 
# 
# 
# 
# 
# idx_cut <- NULL
# 
# 
# for(i in 1:nrow(known_snps)){
#   print(i)
#   chr_temp <- known_snps[i,3]
#   position_temp <- known_snps[i,4]
#   position_low <- position_temp-0.5*10^6
#   position_high <- position_temp+0.5*10^6
#   idx <- which(meta_result_shared_1p$CHR==chr_temp&meta_result_shared_1p$position>position_low&
#                  meta_result_shared_1p$position<position_high)
#   idx_cut <- c(idx_cut,idx)
# }
# 
# meta_result_shared_1p_filter <- meta_result_shared_1p[-idx_cut,]
# 
# save(meta_result_shared_1p_filter,file="/spin1/users/zhangh24/breast_cancer_data_analysis/whole_genome/ICOG/ERPRHER2_fixed/result/meta_result_shared_1p_filter.Rdata")
# 
# 
# 
# 
# 




