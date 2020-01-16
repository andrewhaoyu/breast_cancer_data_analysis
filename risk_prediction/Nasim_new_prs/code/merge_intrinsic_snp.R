setwd('/data/zhangh24/breast_cancer_data_analysis/')
logodds <- NULL
for(i1 in 1:32){
  load( paste0("/data/zhangh24/breast_cancer_data_analysis/risk_prediction/Nasim_new_prs/result/intrinsic_subtype_logodds_dis_",i1,".Rdata"))
  logodds <- rbind(logodds,t(result[[1]]))
}
colnames(logodds) <- c("Luminal_A",
                       "Luminal_B",
                       "Luminal_B_HER2Neg",
                       "HER2_Enriched",
                       "TN")
load(paste0("./risk_prediction/Nasim_new_prs/result/discover_snp_id.rdata"))

final_result <- data.frame(snp_id[,1:3],
                           logodds,
                           stringsAsFactors = F)


save(final_result,file= "./risk_prediction/Nasim_new_prs/result/32_intrinsic_subtype_logodds.Rdata")

