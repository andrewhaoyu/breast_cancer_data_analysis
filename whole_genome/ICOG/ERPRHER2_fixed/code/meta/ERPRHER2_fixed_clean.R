load("/data/zhangh24/breast_cancer_data_analysis/whole_genome/ICOG/ERPRHER2_fixed/result/icog_result_shared_1p.Rdata")
load("/data/zhangh24/breast_cancer_data_analysis/whole_genome/ONCO/ERPRHER2_fixed/result/onco_result_shared_1p.Rdata")
library(bc2)
second.num <- 4


icog_score_infor <- icog_result_shared_1p[,11:(11+second.num+second.num^2-1)]
onco_score_infor <- onco_result_shared_1p[,11:(11+second.num+second.num^2-1)]




icog_onco_score_infor <- cbind(icog_score_infor,onco_score_infor)
meta_result_shared_1p <- icog_result_shared_1p[,c(1:10,31:34)]

rm(icog_result_shared_1p)
rm(onco_result_shared_1p)
gc()

library(foreach)
library(doParallel)
no.cores <- 30
n <- nrow(icog_onco_score_infor)
#icog_onco_score_infor_temp <- icog_onco_score_infor[1:10^5,]
#n <- nrow(icog_onco_score_infor_temp)
#pvalue <- rep(0,n)
size <- 1000
#n <- nrow(icog_onco_score_infor)


  
registerDoParallel(no.cores)
#pvalue <- rep(0,nrow(icog_onco_score_infor))
library(bc2)

pvalue.list <- foreach(i=1:size)%dopar%
{
  print(i)
  start.end <- startend(n,size,i)
  start <- start.end[1]
  end <- start.end[2]
  pvalue_sub <- rep(0,end-start+1)
  temp = 1
  for(j in start:end){
    icog_onco_score_infor_oneline <- icog_onco_score_infor[j,]
    pvalue_sub[temp] <- MetaPfunction(icog_onco_score_infor_oneline,second.num)
    temp = temp+1
  }
  return(pvalue_sub)
 
  
}
stopImplicitCluster()

pvalue <- rep(0,n)
temp.total = 0
for(i in 1:size){
  print(i)
  temp <- length(pvalue.list[[i]])
  pvalue[(1:temp)+temp.total] <- pvalue.list[[i]]
  temp.total <- temp+temp.total
  
  
}



meta_result_shared_1p <- cbind(meta_result_shared_1p,pvalue)





save(meta_result_shared_1p,file=paste0("/data/zhangh24/breast_cancer_data_analysis/whole_genome/ICOG/ERPRHER2_fixed/result/meta_result_shared_1p.Rdata"))

known_snps <- read.csv("/data/zhangh24/breast_cancer_data_analysis/data/known_SNPs_info.csv",header= T)



position.cut <- 10^6



idx_cut <- NULL


for(i in 1:nrow(known_snps)){
  print(i)
  chr_temp <- known_snps[i,3]
  position_temp <- known_snps[i,4]
  position_low <- position_temp-position.cut
  position_high <- position_temp+position.cut
  idx <- which(meta_result_shared_1p$CHR==chr_temp&meta_result_shared_1p$position>position_low&
                 meta_result_shared_1p$position<position_high)
  idx_cut <- c(idx_cut,idx)
}
############duplicate variables won't mater
idx_cut <- unique(idx_cut)
meta_result_shared_1p_filter <- meta_result_shared_1p[-idx_cut,]

save(meta_result_shared_1p_filter,file="/data/zhangh24/breast_cancer_data_analysis/whole_genome/ICOG/ERPRHER2_fixed/result/meta_result_shared_1p_filter_1M.Rdata")

idx <- which(meta_result_shared_1p_filter$pvalue<5E-06)
extract.list <- meta_result_shared_1p_filter[idx,]
save(extract.list,file="/data/zhangh24/breast_cancer_data_analysis/whole_genome/ICOG/ERPRHER2_fixed/result/extract_list.Rdata")





new_filter <- read.csv("/data/zhangh24/breast_cancer_data_analysis/whole_genome/ICOG/ERPRHER2_fixed/result/Filter_based_on_Montse.csv",header=T,stringsAsFactors = F)
new_filter[,2] <- as.numeric(gsub(",","",new_filter[,2]))

idx_cut <- NULL

for(i in 1:nrow(new_filter)){
  print(i)
  chr_temp <- new_filter[i,3]
  position_temp <- new_filter[i,2]
  position_low <- position_temp-position.cut
  position_high <- position_temp+position.cut
  idx <- which(meta_result_shared_1p_filter$CHR==chr_temp&meta_result_shared_1p_filter$position>position_low&
                 meta_result_shared_1p_filter$position<position_high)
  idx_cut <- c(idx_cut,idx)
}
idx_cut <- unique(idx_cut)
meta_result_shared_1p_filter_Ju <- meta_result_shared_1p_filter[-idx_cut,]

save(meta_result_shared_1p_filter_Ju,file="/data/zhangh24/breast_cancer_data_analysis/whole_genome/ICOG/ERPRHER2_fixed/result/meta_result_shared_1p_filter_1M_Ju.Rdata")












