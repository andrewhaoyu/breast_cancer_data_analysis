setwd('/data/zhangh24/breast_cancer_data_analysis/')
logodds <- NULL
var.logodds <- NULL
for(i1 in 1:313){
  load(paste0("/data/zhangh24/breast_cancer_data_analysis/risk_prediction/Nasim_prs/result/intrinsic_subtype_logodds",i1,".Rdata"))
  logodds <- rbind(logodds,t(result[[1]]))
  var.logodds <- rbind(var.logodds,
                       as.vector(result[[2]]))
}
colnames(logodds) <- c("Luminal_A",
                       "Luminal_B",
                       "Luminal_B_HER2Neg",
                       "HER2_Enriched",
                       "TN")
load(paste0("./risk_prediction/Nasim_prs/result/nasim_snp_id.rdata"))

final_result <- data.frame(snp_id,
                           logodds,
                           stringsAsFactors = F)


save(final_result,file= "./risk_prediction/Nasim_prs/result/313_intrinsic_subtype_logodds.Rdata")
#save this results for the LD pruning risk prediction
final_result <- data.frame(snp_id,
                           logodds,
                           var.logodds,
                           stringsAsFactors = F)
save(final_result,file= "./risk_prediction/Nasim_prs/result/313_intrinsic_subtype_logodds_var.Rdata")






#debug : figure out the difference between whole genome intrinsic and extracted snps results
load("./risk_prediction/intrinsic_subtypes_whole_genome/ICOG/result/whole_gonome.rdata")
load("./risk_prediction/intrinsic_subtypes_whole_genome/ICOG/result/icog_result_shared_1p.Rdata")
load("./risk_prediction/intrinsic_subtypes_whole_genome/ONCO/result/onco_result_shared_1p.Rdata")



library(dplyr)
temp_result <- left_join(final_result,whole_genome,
                         by="var_name")
temp_result = temp_result %>% 
  select(var_name,SNP.ICOGS.x,SNP.ONCO.x,
         Luminal_A.x,Luminal_B.x,Luminal_B_HER2Neg.x,HER2_Enriched.x,TN.x,
         Luminal_A.y,Luminal_B.y,Luminal_B_HER2Neg.y,HER2_Enriched.y,TN.y)

colnames(icog_result_shared_1p)
idx <- which(icog_result_shared_1p$var_name=="1_10566215_A_G")
icog_result_shared_1p[idx,]
idx <- which(onco_result_shared_1p$SNP.ONCO=="rs1416885:100010765:C:G")

onco_result_shared_1p[idx,]

