setwd('/data/zhangh24/breast_cancer_data_analysis/')
logodds_overall <- NULL
logodds_erpos <- NULL
logodds_erneg <- NULL
for(i1 in 1:313){
  load(paste0("/data/zhangh24/breast_cancer_data_analysis/risk_prediction/Nasim_prs/result/overall_log_odds_country_",i1,".Rdata"))
  logodds_overall <- c(logodds_overall,result[[1]])
  logodds_erpos <- c(logodds_erpos,result[[3]])
  logodds_erneg <- c(logodds_erneg,result[[5]])
  
}

load(paste0("./risk_prediction/Nasim_prs/result/nasim_snp_id.rdata"))

final_result <- data.frame(snp_id,
                           logodds_overall,
                           logodds_erpos,
                           logodds_erneg,
                           stringsAsFactors = F)


save(final_result,file= "./risk_prediction/Nasim_prs/result/313_overall_logodds_country.Rdata")

