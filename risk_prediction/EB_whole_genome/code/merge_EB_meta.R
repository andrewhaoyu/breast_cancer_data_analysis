setwd('/spin1/users/zhangh24/breast_cancer_data_analysis/risk_prediction')
load('./intrinsic_subtypes_whole_genome/ICOG/result/meta_result_shared_1p.Rdata')
meta_intrin <- meta_result_shared_1p
#load("./ICOG/result/icog_result_shared_1p.Rdata")
#load("./ONCO/result/onco_result_shared_1p.Rdata")
size = 1000
total = 0
for(i1 in 1:size){
  load(paste0("./EB_whole_genome/result/result_sub",i1,".Rdata"))
  #load(paste0("./ICOG/result/result_sub_sub",i1,".Rdata"))  
  total <- total + nrow(result_sub[[1]])
}
second.num <- 5
log.odds <- matrix(0,total,second.num)
#covar <- rep(0,total)
heter.var <- rep(0,total)

total <- 0
for(i1 in 1:size){
  load(paste0("./EB_whole_genome/result/result_sub",i1,".Rdata"))
  temp <- nrow(result_sub[[1]])
  log.odds[total+(1:temp),] <- result_sub[[1]]
  heter.var[total+(1:temp)] <- result_sub[[2]]
  total <- total + temp
}

meta_intrin[,11:15] <- log.odds
meta_eb <- meta_intrin[,1:15]
meta_eb <- cbind(meta_eb,heter.var)
# icog_result_shared_1p[,16:40] <- covar
# meta_result_shared_1p <- cbind(icog_result_shared_1p,p.value) 
save(meta_eb,
     file="./EB_whole_genome/result/meta_result_shared_1p.Rdata")
