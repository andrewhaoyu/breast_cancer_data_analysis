setwd('/spin1/users/zhangh24/breast_cancer_data_analysis/risk_prediction')
load('./FTOP_whole_genome/ICOG/result/meta_result_shared_1p.Rdata')
meta_FTOP <- meta_result_shared_1p_FTOP
#load("./ICOG/result/icog_result_shared_1p.Rdata")
#load("./ONCO/result/onco_result_shared_1p.Rdata")
size = 1000
total = 0
for(i1 in 1:size){
  load(paste0("./EB_whole_genome/result/p_value_sub",i1,".Rdata"))
  #load(paste0("./ICOG/result/result_sub_sub",i1,".Rdata"))  
  total <- total + nrow(result_sub[[1]])
}
second.num <- 5
log.odds <- matrix(0,total,second.num)
#covar <- rep(0,total)
heter.var <- rep(0,total)

total <- 0
for(i1 in 1:size){
  load(paste0("./EB_whole_genome/result/p_value_sub",i1,".Rdata"))
  temp <- nrow(result_sub[[1]])
  log.odds[total+(1:temp)] <- result_sub[[1]]
  heter.var[total+(1:temp)] <- result_sub[[2]]
  total <- total + temp
}

icog_result_shared_1p[,11] <- log.odds
icog_result_shared_1p[,12] <- covar
meta_result_shared_1p <- cbind(icog_result_shared_1p,p.value) 
save(meta_result_shared_1p,
     file="./ICOG/result/meta_result_shared_1p.Rdata")
