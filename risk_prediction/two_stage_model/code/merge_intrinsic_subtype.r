setwd("/spin1/users/zhangh24/breast_cancer_data_analysis/risk_prediction/two_stage_model/result/")
log.odds.meta.two.stage
for(i1 in 1:205){
print(i1)
    load(paste0("/spin1/users/zhangh24/breast_cancer_data_analysis/risk_prediction/two_stage_model/result/meta.result",i1,".Rdata"))
}