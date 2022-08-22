#add standard error to the supplementary table 19
load("./risk_prediction/Nasim_prs/result/313_intrinsic_subtype_logodds_var.Rdata")
final_result_all = final_result
load("./risk_prediction/Nasim_new_prs/result/32_intrinsic_subtype_logodds_var.Rdata")
load("./risk_prediction/Nasim_prs/result/independent_SNPs_idx_in_32_novel_SNP.rdata")
final_result_all<- rbind(final_result_all,final_result[idx,])

#
sd_mat = matrix(0,nrow(final_result_all),5)
for(l in 1:nrow(final_result_all)){
  sd_mat[l,] = sqrt(as.numeric(diag(matrix(final_result_all[l,9:33],5,5))))
}
final_result_all = final_result_all[,c(1:8)]
final_result_all <- cbind(final_result_all,sd_mat)
write.csv(final_result_all,file ="./risk_prediction/Nasim_new_prs/result/updated_sup_table19.csv" )
