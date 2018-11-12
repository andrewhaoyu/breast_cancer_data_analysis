args = commandArgs(trailingOnly = T)
i1 = as.numeric(args[[1]])

setwd('/spin1/users/zhangh24/breast_cancer_data_analysis/risk_prediction/intrinsic_subtypes_whole_genome/')

load("./ICOG/result/icog_result_shared_1p.Rdata")
load("./ONCO/result/onco_result_shared_1p.Rdata")
# load("/spin1/users/zhangh24/breast_cancer_data_analysis/whole_genome_age/ICOG/ERPRHER2GRADE_fixed_baseline/result/icog_result_only_shared_1p.Rdata")
# load("/spin1/users/zhangh24/breast_cancer_data_analysis/whole_genome_age/ONCO/ERPRHER2GRADE_fixed_baseline/result/onco_result_only_shared_1p.Rdata")
library(bc2)
second.num <-5


icog_score_infor <- icog_result_shared_1p[,11:(11+second.num+second.num^2-1)]
onco_score_infor <- onco_result_shared_1p[,11:(11+second.num+second.num^2-1)]

icog_onco_score_infor <- cbind(icog_score_infor,onco_score_infor)
n <- nrow(icog_onco_score_infor)
rm(icog_result_shared_1p)
rm(onco_result_shared_1p)
gc()
# 
size <- 1000
Metaoddsfunction <- function(icog_onco_score_infor_one,second.num){
  score.icog <- as.numeric(icog_onco_score_infor_one[1:(second.num)])
  infor.icog <- matrix(as.numeric(icog_onco_score_infor_one[(second.num+1):
                                                              (second.num+second.num^2) ]),
                       ncol = second.num)
  start <- second.num+second.num^2
  score.onco <- as.numeric(icog_onco_score_infor_one[(1+start):
                                                       (second.num+start)])
  infor.onco <- matrix(as.numeric(icog_onco_score_infor_one[(second.num+1+start):
                                                              (second.num+second.num^2+start) ]),ncol=second.num)
  
  
  meta.result <- LogoddsMetaAnalysis(score.icog,infor.icog,
                                     score.onco,infor.onco)
  score.meta <- as.vector(meta.result[[1]])
  infor.meta <- meta.result[[2]]
  p.meta <- DisplaySecondStageTestResult(score.meta,infor.meta)[[second.num*2+1]]
  return(list(log.odds=score.meta,
              covar = infor.meta,
              p = p.meta))
}


#registerDoParallel(no.cores)
#pvalue <- rep(0,nrow(icog_onco_score_infor))


#pvalue.list <- foreach(i=1:size)%dopar%
#{
#print(i)
start.end <- startend(n,size,i1)
start <- start.end[1]
end <- start.end[2]
logodds_sub <- matrix(0,end-start+1,second.num)
covar_sub <- matrix(0,end-start+1,second.num^2)
pvalue_sub <- rep(0,end-start+1)
temp = 1
for(j in start:end){
  print(j)
  
  icog_onco_score_infor_oneline <- icog_onco_score_infor[j,]
  result_temp <- Metaoddsfunction(icog_onco_score_infor_oneline,second.num)
  logodds_sub[temp,] <- result_temp[[1]]
  covar_sub[temp,] <- as.vector(result_temp[[2]])
  pvalue_sub[temp] <- result_temp[[3]]
  temp = temp+1
}


result_sub <- list(logodds_sub,
                   covar_sub,
                   as.numeric(pvalue_sub))
save(result_sub,file=paste0("./ICOG/result/result_sub",i1,".Rdata"))
