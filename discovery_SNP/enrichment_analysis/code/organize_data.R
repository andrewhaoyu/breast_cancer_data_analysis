setwd("/spin1/users/zhangh24/breast_cancer_data_analysis/")
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
enrichment_analysis_data <- enrichment_result
enrichment_analysis_data <- as.data.frame(enrichment_analysis_data)
#save(enrichment_analysis_data,file = "./discovery_SNP/enrichment_analysis/result/enrichment_analysis_data.Rdata")

head(enrichment_analysis_data)

Pfunction <- function(z){
  return(2*pnorm(-abs(z)))
}


n <- nrow(enrichment_analysis_data)
p_value2 <- p_value1 <- rep(0,n)

# for(i in 1:n){
#   if(i%%1000==0){
#     print(i)
#   }
  p_value1 <- Pfunction(enrichment_analysis_data[,8]/sqrt(enrichment_analysis_data[,10]))
  p_value2 <- Pfunction(enrichment_analysis_data[,9]/sqrt(enrichment_analysis_data[,13]))
#}

enrichment_analysis_data <- data.frame(enrichment_analysis_data,p_value1,p_value2)
colnames(enrichment_analysis_data)[c(14,15)] <- c("Luminal_A_p_value",
                                                  "triple_negative_p_value")
#save(enrichment_analysis_data,file = "./discovery_SNP/enrichment_analysis/result/enrichment_analysis_data.Rdata")
#load("./discovery_SNP/enrichment_analysis/result/enrichment_analysis_data.Rdata")

n <- nrow(enrichment_analysis_data)
z_LuA <- rep(0,n)
z_tn <- rep(0,n)
n_LuA <- rep(0,n)
n_tn <- rep(0,n)
z_LuA <- enrichment_analysis_data[,8]/sqrt(enrichment_analysis_data[,10])
z_tn <- enrichment_analysis_data[,9]/sqrt(enrichment_analysis_data[,13])
  n_LuA <- 1/(enrichment_analysis_data[,10]*2*enrichment_analysis_data[,6]*(1-enrichment_analysis_data[,6]))
  n_tn <- 1/(enrichment_analysis_data[,13]*2*enrichment_analysis_data[,6]*(1-enrichment_analysis_data[,6]))
  enrichment_analysis_data <- cbind(enrichment_analysis_data,z_LuA,z_tn,n_LuA,n_tn)
  head(enrichment_analysis_data)
  colnames(enrichment_analysis_data)[16:19] <- c("Z_Luminal_A",
                                                 "Z_triple_neg",
                                                 "sample_size_Luminal_A",
                                                 "sample_size_triple_neg")
  save(enrichment_analysis_data,file="./discovery_SNP/enrichment_analysis/result/enrichment_analysis_data.Rdata")
#}

# matrix(as.numeric(enrichment_analysis_data[1,10:13]),2,2)
# 
# ##############combine CIMBA and BCAC
# combine.data <- merge(CIMBA,meta_result_shared_1p,by.x = "MarkerName",
#                       by.y = "var_name")
# CIMBA.BCAC.combine <- combine.data
# ##############save BCAC and CIMBA combine data for meta-analysis
# save(CIMBA.BCAC.combine,file = paste0("./discovery_SNP/CIMBA_BCAC_meta_analysis/result/CIMBA.BCAC.combine.Rdata"))
# 
# 
