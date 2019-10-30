#merge meta-analysis results for overall analysis and subtypes analysis without adjusting for country
second.num <- 6
setwd("/data/zhangh24/breast_cancer_data_analysis/")
load("./whole_genome_age/ICOG/standard_analysis/result/icog_result_shared_1p_s.Rdata")
load("./whole_genome_age/ONCO/standard_analysis/result/onco_result_shared_1p_s.Rdata")
load("./whole_genome_age/ONCO/standard_analysis/result/icog_result_only_shared_1p_s.Rdata")
load("./whole_genome_age/ONCO/standard_analysis/result/onco_result_only_shared_1p_s.Rdata")


meta_result_shared_1p <- icog_result_shared_1p[,c(1:10,(ncol(icog_result_shared_1p)-3):ncol(icog_result_shared_1p))]
meta_result_shared_1p_icog_only <- icog_result_only_shared_1p[,c(1:10,(ncol(icog_result_shared_1p)-3):ncol(icog_result_shared_1p))]
meta_result_shared_1p_onco_only <- onco_result_only_shared_1p[,c(1:10,(ncol(icog_result_shared_1p)-3):ncol(icog_result_shared_1p))]


meta_result_shared_1p_no_pvalue <- rbind(meta_result_shared_1p,meta_result_shared_1p_icog_only,meta_result_shared_1p_onco_only)


# save(meta_result_shared_1p_no_pvalue,file="/data/zhangh24/breast_cancer_data_analysis/whole_genome_age/ICOG/Intrinsic_subtypes/result/meta_result_shared_1p_no_pvalue.Rdata")


n <- nrow(meta_result_shared_1p_no_pvalue)


total <- 0
log.odds <- matrix(0,n,second.num)
sigma <- matrix(0,n,second.num)
p.value <- matrix(0,n,second.num)
for(i1 in 1:1000){
  print(i1)  
  load(paste0("./whole_genome_age/ICOG/ERPRHER2GRADE_casecase/result/result_sub_s",i1,".Rdata"))
  temp <- nrow(result_summary)
  log.odds.temp <- result_summary[,1:second.num]
  sig.temp <- result_summary[,(1+second.num):(2*second.num)]
  z = log.odds.temp/sqrt(sig.temp)
  p.temp = 2*pnorm(-abs(z),lower.tail = T)
  log.odds[total+(1:temp),] <- log.odds.temp
  sigma[total+(1:temp),] <- sig.temp
  p.value[total+(1:temp),] <- p.temp
  total <- total + temp
  
}


meta_result_shared_1p <- cbind(meta_result_shared_1p_no_pvalue,log.odds,sigma,p.value)
save(meta_result_shared_1p,file=paste0("./whole_genome_age/ICOG/standard_analysis/result/meta_result_shared_1p_s.Rdata"))





