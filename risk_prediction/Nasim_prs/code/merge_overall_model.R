setwd('/data/zhangh24/breast_cancer_data_analysis/')
logodds_overall <- NULL
var_overall <- NULL
logodds_erpos <- NULL
logodds_erneg <- NULL
for(i1 in 1:313){
  load(paste0("/data/zhangh24/breast_cancer_data_analysis/risk_prediction/Nasim_prs/result/overall_log_odds_",i1,".Rdata"))
  logodds_overall <- c(logodds_overall,result[[1]])
  var_overall <- c(var_overall,
                   result[[2]])
  logodds_erpos <- c(logodds_erpos,result[[3]])
  logodds_erneg <- c(logodds_erneg,result[[5]])
  
}

load(paste0("./risk_prediction/Nasim_prs/result/nasim_snp_id.rdata"))

final_result <- data.frame(snp_id,
                           logodds_overall,
                           logodds_erpos,
                           logodds_erneg,
                           stringsAsFactors = F)


save(final_result,file= "./risk_prediction/Nasim_prs/result/313_overall_logodds.Rdata")

#save the results for future LD pruning analysis due to 3 of the 313 SNPs having MAF lower than 0.01
final_result <- data.frame(snp_id,
                           logodds_overall,
                           var_overall,
                           
                           stringsAsFactors = F)
save(final_result,file= "./risk_prediction/Nasim_prs/result/313_overall_logodds_var.Rdata")
