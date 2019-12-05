args = commandArgs(trailingOnly = T)
i1 = as.numeric(args[[1]])

setwd('/data/zhangh24/breast_cancer_data_analysis/risk_prediction/FTOP_whole_genome/')
load("./ICOG/result/icog_result_shared_1p.Rdata")
load("./ONCO/result/onco_result_shared_1p.Rdata")
library(bc2, lib.loc ="/home/zhangh24/R/x86_64-pc-linux-gnu-library/3.6/")
second.num <- 5


icog_score_infor <- icog_result_shared_1p[,11:(11+second.num+second.num^2-1)]
onco_score_infor <- onco_result_shared_1p[,11:(11+second.num+second.num^2-1)]





icog_onco_score_infor <- cbind(icog_score_infor,onco_score_infor)

rm(icog_result_shared_1p)
rm(onco_result_shared_1p)
gc()
n <- nrow(icog_onco_score_infor)
size <- 1000
start.end <- startend(n,size,i1)
start <- start.end[1]
end <- start.end[2]
pvalue_sub <- rep(0,end-start+1)
temp = 1
for(j in start:end){
  print(j)
  
  icog_onco_score_infor_oneline <- icog_onco_score_infor[j,]
  pvalue_sub[temp] <- MetaPfunction(icog_onco_score_infor_oneline,second.num)
  temp = temp+1
}



save(pvalue_sub,file=paste0("./ICOG/result/p_value_sub",i1,".Rdata"))
















