setwd('/data/zhangh24/breast_cancer_data_analysis/risk_prediction/standard_whole_genome/')
load("./ICOG/result/icog_result_shared_1p.Rdata")
load("./ONCO/result/onco_result_shared_1p.Rdata")
size = 1000
total = 0
for(i1 in 1:size){
  load(paste0("./ICOG/result/result_sub_sub",i1,".Rdata"))  
  total <- total + length(result_sub[[1]])
}

log.odds <- rep(0,total)
covar <- rep(0,total)
p.value <- rep(0,total)

total <- 0
for(i1 in 1:size){
  load(paste0("./ICOG/result/result_sub_sub",i1,".Rdata"))  
  temp <- length(result_sub[[1]])
  log.odds[total+(1:temp)] <- result_sub[[1]]
  covar[total+(1:temp)] <- result_sub[[2]]
  p.value[total+(1:temp)] <- result_sub[[3]]
  total <- total + temp
}

icog_result_shared_1p[,11] <- log.odds
icog_result_shared_1p[,12] <- covar
meta_result_shared_1p <- cbind(icog_result_shared_1p,p.value) 
save(meta_result_shared_1p,
     file="./ICOG/result/meta_result_shared_1p.Rdata")
