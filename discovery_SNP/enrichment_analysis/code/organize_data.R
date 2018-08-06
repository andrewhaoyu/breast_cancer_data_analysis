
load(paste0("./discovery_SNP/CIMBA_BCAC_meta_analysis/result/CIMBA.BCAC.meta.result.Rdata"))
enrichment_result <- CIMBA.BCAC.meta.result[,c(9,18,10,2,3,4,12,22,26,
                                               27,31,47,51)]
colnames(enrichment_result)[c(1,4,5,6,7,8,9)] <- c("rs_id",
                                                "Reference_allele",
                                                "Effect_allele",
                                                "Freq_Reference_allele",
                                                "imputation_info",
                                                "Log_OR_Luminal_A",
                                                "Log_OR_Triple_Neg"
                                                
  
)



##############combine CIMBA and BCAC
combine.data <- merge(CIMBA,meta_result_shared_1p,by.x = "MarkerName",
                      by.y = "var_name")
CIMBA.BCAC.combine <- combine.data
##############save BCAC and CIMBA combine data for meta-analysis
save(CIMBA.BCAC.combine,file = paste0("./discovery_SNP/CIMBA_BCAC_meta_analysis/result/CIMBA.BCAC.combine.Rdata"))


