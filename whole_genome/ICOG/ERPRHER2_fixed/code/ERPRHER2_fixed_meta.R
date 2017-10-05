load("/spin1/users/zhangh24/breast_cancer_data_analysis/whole_genome/ICOG/ERPRHER2_fixed/result/icog_result_shared_1p.Rdata")
load("/spin1/users/zhangh24/breast_cancer_data_analysis/whole_genome/ONCO/ERPRHER2_fixed/result/onco_result_shared_1p.Rdata")
library(bc2)
second.num <- 4


icog_score_infor <- icog_result_shared_1p[,11:(11+second.num+second.num^2-1)]
onco_score_infor <- onco_result_shared_1p[,11:(11+second.num+second.num^2-1)]




icog_onco_score_infor <- cbind(icog_score_infor,onco_score_infor)
meta_result_shared_1p <- icog_result_shared_1p[,c(1:10,31:33)]
rm(icog_result_shared_1p)
rm(onco_result_shared_1p)
gc()

MetaPfunction <- function(icog_onco_score_infor_one){
  score.icog <- as.numeric(icog_onco_score_infor_one[1:(second.num)])
  infor.icog <- matrix(as.numeric(icog_onco_score_infor_one[(second.num+1):
                                                  (second.num+second.num^2) ]),
                       ncol = second.num)
  start <- second.num+second.num^2
  score.onco <- as.numeric(icog_onco_score_infor_one[(1+start):
                                            (second.num+start)])
  infor.onco <- matrix(as.numeric(icog_onco_score_infor_one[(second.num+1+start):
                                                   (second.num+second.num^2+start) ]),ncol=second.num)
  
  
  meta.result <- ScoreMetaAnalysis(score.icog,infor.icog,
                                         score.onco,infor.onco)
  score.meta <- t(meta.result[[1]])
  infor.meta <- meta.result[[2]]
  DisplayFixedScoreTestResult(score.meta,infor.meta)
}


#icog_onco_score_infor_temp <- icog_onco_score_infor[1:10000,]

library(foreach)
library(doParallel)
no.cores <- 20
registerDoParallel(no.cores)
pvalue <- foreach(i=1:nrow(icog_onco_score_infor),
                  .combine=c)%dopar%
{
  icog_onco_score_infor_oneline <- icog_onco_score_infor[i,]
  return(MetaPfunction(icog_onco_score_infor_oneline))
}
stopImplicitCluster()
















#p.value <- apply(icog_onco_score_infor_temp,1,MetaPfunction)


#meta_result_shared_1p_temp <- meta_result_shared_1p[1:10000,]

meta_result_shared_1p <- cbind(meta_result_shared_1p,pvalue)


save(meta_result_shared_1p,file=paste0("/spin1/users/zhangh24/breast_cancer_data_analysis/whole_genome/ICOG/ERPRHER2_fixed/result/meta_result_shared_1p.Rdata"))





