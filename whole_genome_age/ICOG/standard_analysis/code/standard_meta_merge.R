second.num <- 1
setwd("/spin1/users/zhangh24/breast_cancer_data_analysis/")
load("./whole_genome_age/ICOG/standard_analysis/result/icog_result_shared_1p.Rdata")
load("./whole_genome_age/ONCO/standard_analysis/result/onco_result_shared_1p.Rdata")
load("./whole_genome_age/ONCO/standard_analysis/result/icog_result_only_shared_1p.Rdata")
load("./whole_genome_age/ONCO/standard_analysis/result/onco_result_only_shared_1p.Rdata")


meta_result_shared_1p <- icog_result_shared_1p[,c(1:10,(ncol(icog_result_shared_1p)-3):ncol(icog_result_shared_1p))]
meta_result_shared_1p_icog_only <- icog_result_only_shared_1p[,c(1:10,(ncol(icog_result_shared_1p)-3):ncol(icog_result_shared_1p))]
meta_result_shared_1p_onco_only <- onco_result_only_shared_1p[,c(1:10,(ncol(icog_result_shared_1p)-3):ncol(icog_result_shared_1p))]


meta_result_shared_1p_no_pvalue <- rbind(meta_result_shared_1p,meta_result_shared_1p_icog_only,meta_result_shared_1p_onco_only)


# save(meta_result_shared_1p_no_pvalue,file="/spin1/users/zhangh24/breast_cancer_data_analysis/whole_genome_age/ICOG/Intrinsic_subtypes/result/meta_result_shared_1p_no_pvalue.Rdata")


n <- nrow(meta_result_shared_1p_no_pvalue)

p.value <- rep(0,n)
total <- 0
log.odds <- matrix(0,n,second.num)
sigma <- matrix(0,n,second.num^2)
for(i1 in 1:1000){
  print(i1)  
  load(paste0("./whole_genome_age/ICOG/ERPRHER2GRADE_casecase/result/result_sub",i1,".Rdata"))
  temp <- nrow(result_summary)
  
  p.value[total+(1:temp)] <- result_summary[,3]
  log.odds[total+(1:temp),] <- result_summary[,1]
  sigma[total+(1:temp),] <- result_summary[,2]
  total <- total + temp
  
}


meta_result_shared_1p <- cbind(meta_result_shared_1p_no_pvalue,p.value,log.odds,sigma)
save(meta_result_shared_1p,file=paste0("./whole_genome_age/ICOG/standard_analysis/result/meta_result_shared_1p.Rdata"))





