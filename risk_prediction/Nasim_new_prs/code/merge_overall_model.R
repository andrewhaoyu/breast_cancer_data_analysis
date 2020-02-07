setwd('/data/zhangh24/breast_cancer_data_analysis/')
logodds_overall <- NULL
var_overall <- NULL
for(i1 in 1:32){
  load(paste0("/data/zhangh24/breast_cancer_data_analysis/risk_prediction/Nasim_new_prs/result/overall_log_odds_",i1,".Rdata"))
  logodds_overall <- c(logodds_overall,result[[1]])
  var_overall <- c(var_overall,
                   result[[2]])

}

load(paste0("./risk_prediction/Nasim_new_prs/result/discover_snp_id.rdata"))

final_result <- data.frame(snp_id[,1:3],
                           logodds_overall,
                           var_overall,
                           stringsAsFactors = F)


save(final_result,file= "./risk_prediction/Nasim_new_prs/result/32_overall_logodds.Rdata")
