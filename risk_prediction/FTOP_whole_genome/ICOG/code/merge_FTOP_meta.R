setwd("/data/zhangh24/breast_cancer_data_analysis/risk_prediction/FTOP_whole_genome")
load("./ICOG/result/icog_result_shared_1p.Rdata")
load("./ONCO/result/onco_result_shared_1p.Rdata")
size = 1000
total = 0
for(i1 in 1:size){
  
  load(paste0("./ICOG/result/p_value_sub",i1,".Rdata"))  
  total <- total + length(pvalue_sub)
}


p.value <- rep(0,total)

total <- 0
for(i1 in 1:size){
  load(paste0("./ICOG/result/p_value_sub",i1,".Rdata"))  
  temp <- length(pvalue_sub)

  p.value[total+(1:temp)] <- pvalue_sub
  total <- total + temp
}

#icog_result_shared_1p[,11] <- log.odds
#icog_result_shared_1p[,12] <- covar
meta_result_shared_1p_FTOP <- cbind(icog_result_shared_1p[c(1:10,41:44)],
                               p.value) 
#meta_result_shared_1p_FTOP <- meta_result_shared_1p
save(meta_result_shared_1p_FTOP,
     file="./ICOG/result/meta_result_shared_1p.Rdata")
