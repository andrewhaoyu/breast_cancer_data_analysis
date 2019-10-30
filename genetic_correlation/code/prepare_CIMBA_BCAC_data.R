setwd('/data/zhangh24/breast_cancer_data_analysis/')
library(data.table)
###########load CIMBA data
load(paste0("./discovery_SNP/CIMBA_BCAC_meta_analysis/result/CIMBA.clean.Rdata"))
CIMBA <- CIMBA.clean
###########load BCAC intrinsic subtype data
load(paste0("/data/zhangh24/breast_cancer_data_analysis/whole_genome_age/ICOG/Intrinsic_subtypes/result/meta_result_shared_1p_082119.Rdata"))


N.TN <- 1/(meta_result_shared_1p[,45]*2*meta_result_shared_1p[,4]*(1-meta_result_shared_1p[,4]))
median(N.TN)
N.TN.CIMBA <- 1/(CIMBA[,6]^2*2*CIMBA[,4]*(1-CIMBA[,4]))
head(N.TN.CIMBA)
median(N.TN.CIMBA)

##########filter out the SNP with imputation quality more than 0.3
##########data for genetic correlation of BCAC
idx <- which(meta_result_shared_1p[,5]>=0.3)
meta_result_shared_1p <- meta_result_shared_1p[idx,]
colnames(meta_result_shared_1p)
snp.split <- strsplit(meta_result_shared_1p[,12],
                      split = "_")
n <- nrow(meta_result_shared_1p)
references_allele <- rep("c",n)
effect_allele <- rep("c",n)
for(i in 1:n){
  if(i%%10000==0){
    print(i)  
  }
  references_allele[i] <- snp.split[[i]][3]
  effect_allele[i] <- snp.split[[i]][4]
}
meta_result_shared_1p$references_allele <- references_allele
meta_result_shared_1p$effect_allele <- effect_allele
colnames(meta_result_shared_1p)
#try <- (meta_result_shared_1p[,c(1,7,13,19,25)+20]),2,sd)
z.stat <- meta_result_shared_1p[,16:20]/sqrt(meta_result_shared_1p[,c(1,7,13,19,25)+20])
p.value <- z.stat
for(i in 1:5){
  p.value[,i] <- 2*pnorm(abs(z.stat[,i]),lower.tail = F)  
}


N <- 1/(meta_result_shared_1p[,c(1,7,13,19,25)+20]*2*meta_result_shared_1p[,4]*(1-meta_result_shared_1p[,4]))
meta_result_shared_1p <- meta_result_shared_1p[,c(2,11,3,4,46,47,16:20,c(1,7,13,19,25)+20)]
BCAC_subtypes_result <- data.frame(meta_result_shared_1p,
                                    z.stat,
                                    N,
                                    p.value)
subtypes <- c("Luminial_A","Luminal_B",
                                    "Luminal_B_HER2Neg",
                                    "HER2_Enriched",
                                    "Triple_Negative") %>% %>% 
colnames(BCAC_subtypes_result)[7:31] <- c(paste0("log_or_",subtypes),
                                           paste0("var_",subtypes),
                                           paste0("z_stat_",subtypes),
                                           paste0("effect_sample_size_",subtypes),
                                           paste0("p_value_",subtypes))
colnames(BCAC_subtypes_result)
save(BCAC_subtypes_result,file= paste0("./genetic_correlation/result/BCAC_subtypes_result_082119.Rdata"))
head(meta_result_shared_1p)

#############prepare data for CIMBA to genetic correlation analysis
z.stat <- CIMBA[,5]/CIMBA[,6]
N <- 1/(CIMBA[,6]^2*2*CIMBA[,4]*(1-CIMBA[,4]))
p.value <- 2*pnorm(abs(z.stat),lower.tail = F)
CIMBA.result <- data.frame(CIMBA[,1:5],
                           CIMBA[,6]^2,
                           z.stat,
                           N,
                           p.value)

colnames(CIMBA.result) <- c("MarkerName",
                            "ReferenceAllele",
                            "EffectAllele",
                            "Freq_effect_allele",
                            "LogOR",
                            "Var",
                            "Zstat",
                            "EffectSampleSize",
                            "pvalue"
                            )
save(CIMBA.result,file = "./genetic_correlation/result/CIMBA_result.Rdata")


###############prepare data for breast cancer standard to genetic correlation analysis
load("/data/zhangh24/breast_cancer/standard_gwas/case_control/meta/case_control_meta_result_final.Rdata")
load("/data/zhangh24/breast_cancer/standard_gwas/ER-_control/meta/ER-_control_meta_result_final.Rdata")
load("/data/zhangh24/breast_cancer/standard_gwas/ER+_control/meta/ER+_control_meta_result_final.Rdata")


standard.analysis <- data.frame(case_control_meta_result_final,
                                ERNeg_control_meta_result_final,
                                ERPos_control_meta_result_final)



standard.analysis.m <- merge(standard.analysis,
                                          BCAC_subtypes_result,
                                          by.x= "rs_id",
                                          by.y = "rs_id")
standard.analysis.result <- standard.analysis.m[,c(1:5,9:10,14:15,18:20)]

idx <- which(is.na(standard.analysis.result[,5]))
standard.analysis.result <- standard.analysis.result[-idx,]
z.stat <- cbind(standard.analysis.result[,4]/standard.analysis.result[,5],
                standard.analysis.result[,6]/standard.analysis.result[,7],
                standard.analysis.result[,8]/standard.analysis.result[,9])
p.value <- z.stat
for(i in 1:3){
  p.value[,i] <- 2*pnorm(abs(z.stat[,i]),lower.tail=F)
}
N <- cbind((1/(standard.analysis.result[,5]^2*2*standard.analysis.result[,10]*(1-standard.analysis.result[,10]))),
           (1/(standard.analysis.result[,7]^2*2*standard.analysis.result[,10]*(1-standard.analysis.result[,10]))),
           (1/(standard.analysis.result[,9]^2*2*standard.analysis.result[,10]*(1-standard.analysis.result[,10]))))


apply(N,2,median)


head(standard.analysis.result)
standard.analysis.result.c <- data.frame(standard.analysis.result[,c(1,2,3,10,11,12)],
                                         standard.analysis.result[,4],
                                         standard.analysis.result[,5]^2,
                                         standard.analysis.result[,6],
                                         standard.analysis.result[,7]^2,
                                         standard.analysis.result[,8],
                                         standard.analysis.result[,9]^2,
                                         z.stat,
                                         N,
                                         p.value
                                         )
colnames(standard.analysis.result.c) <- c("rs_id",
                                          "CHR",
                                          "position",
                                          "freq_a1",
                                          "reference_allele",
                                          "effect_allele",
                                          "overall_breast_cancer_log_or",
                                          "overall_braest_cancer_var",
                                          "ER_neg_log_or",
                                          "ER_neg_var",
                                          "ER_pos_log_or",
                                          "ER_pos_var",
                                          "z_stat_overall_breast_cancer",
                                          "z_stat_ER_neg",
                                          "z_stat_ER_pos",
                                          "Effect_sample_size_overall_breast_cancer",
                                          "Effect_sample_size_ER_neg",
                                          "Effect_sample_size_ER_pos",
                                          "p_value_overall_breast_cancer",
                                          "p_value_ER_neg",
                                          "p_value_ER_pos")
standard.analysis.result <- standard.analysis.result.c
save(standard.analysis.result,file="./genetic_correlation/result/standard_analysis_result.Rdata")

####################CIMBA and BCAC meta analysis genetic correlation
load(paste0("./discovery_SNP/CIMBA_BCAC_meta_analysis/result/CIMBA.BCAC.meta.result.Rdata"))
head(CIMBA.BCAC.meta.result)


idx <- which(CIMBA.BCAC.meta.result[,11]>=0.3)
CIMBA.BCAC.meta.result <- CIMBA.BCAC.meta.result[idx,]
CIMBA.BCAC.meta.result <- as.data.frame(CIMBA.BCAC.meta.result)
#try <- (meta_result_shared_1p[,c(1,7,13,19,25)+20]),2,sd)
z.stat <- CIMBA.BCAC.meta.result[,21:25]/sqrt(CIMBA.BCAC.meta.result[,c(1,7,13,19,25)+25])
p.value <- z.stat
for(i in 1:5){
  p.value[,i] <- 2*pnorm(abs(z.stat[,i]),lower.tail = F)  
}


N <- 1/(CIMBA.BCAC.meta.result[,c(1,7,13,19,25)+25]*2*CIMBA.BCAC.meta.result[,4]*(1-CIMBA.BCAC.meta.result[,4]))


CIMBA.BCAC.meta.result.m<- CIMBA.BCAC.meta.result[,c(8,17,9,4,2,3,21:25,c(1,7,13,19,25)+25)]
head(CIMBA.BCAC.meta.result.m)
CIMBA.BCAC.meta.result <- data.frame(CIMBA.BCAC.meta.result.m,
                                   z.stat,
                                   N,
                                   p.value)
colnames(CIMBA.BCAC.meta.result)
subtypes <- c("Luminial_A","Luminal_B",
              "Luminal_B_HER2Neg",
              "HER2_Enriched",
              "Triple_Negative")
colnames(CIMBA.BCAC.meta.result)[7:31] <- c(paste0("log_or_",subtypes),
                                          paste0("var_",subtypes),
                                          paste0("z_stat_",subtypes),
                                          paste0("effect_sample_size_",subtypes),
                                          paste0("p_value_",subtypes))
colnames(CIMBA.BCAC.meta.result)
save(CIMBA.BCAC.meta.result,file= paste0("./genetic_correlation/result/CIMBA.BCAC.meta.result.Rdata"))
head(meta_result_shared_1p)







#
# new.data <- merge(CIMBA.BCAC.meta.result,
#                   enrichment_analysis_data,
#                   by.x ="rs_id",
#                   by.y = "rs_id")
# all.equal(new.data$log_or_Triple_Negative,new.data$Log_OR_Triple_Neg)

#save()




################change CIMBA MarkerName into BCAC rs_id style
#CIMBA_rs_id <- gsub("_",":",CIMBA$MarkerName)
#head(CIMBA_rs_id)
#CIMBA$rs_id <- CIMBA_rs_id
# head(CIMBA)
# colnames(meta_result_shared_1p)[21:45] <- paste0("var",c(1:25))
# ##############combine CIMBA and BCAC
# combine.data <- merge(CIMBA,meta_result_shared_1p,by.x = "MarkerName",
#                       by.y = "var_name")
# CIMBA.BCAC.combine <- combine.data
# ##############save BCAC and CIMBA combine data for meta-analysis
# save(CIMBA.BCAC.combine,file = paste0("./discovery_SNP/CIMBA_BCAC_meta_analysis/result/CIMBA.BCAC.combine.Rdata"))
