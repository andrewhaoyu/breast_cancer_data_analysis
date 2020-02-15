library(data.table)
#data <- as.data.frame(fread("./data/oncoarray_bcac_public_release_oct17.txt"))

#load results of standard analysis
library(data.table)
setwd('/data/zhangh24/breast_cancer_data_analysis/')

all.snp <- as.data.frame(fread("./discovery_SNP/result/ResultsMeta_GWAS_iCOGs_Onco_filter_R2_MAF.txt",header=T))
#load icogs intrinsic subtypes result
load("/data/zhangh24/breast_cancer_data_analysis/whole_genome_age/ICOG/Intrinsic_subtypes/result/icog_result_shared_1p_082119.Rdata")
load("/data/zhangh24/breast_cancer_data_analysis/whole_genome_age/ICOG/Intrinsic_subtypes/result/icog_result_only_shared_1p_082119.Rdata")
icog_result <- rbind(icog_result_shared_1p,icog_result_only_shared_1p)

load("/data/zhangh24/breast_cancer_data_analysis/whole_genome_age/ICOG/Intrinsic_subtypes/result/icog_result_shared_1p_sex.Rdata")
icog_result_sex <- icog_result_shared_1p
#some column may not match
#but it doesn't matter since we only select the matched column latter
colnames(icog_result_sex) = colnames(icog_result)
icog_result = rbind(icog_result,
                    icog_result_sex )
#load oncoarray intrinsic subtypes result
load("/data/zhangh24/breast_cancer_data_analysis/whole_genome_age/ONCO/intrinsic_subtypes/result/onco_result_shared_1p_082119.Rdata")
load("/data/zhangh24/breast_cancer_data_analysis/whole_genome_age/ONCO/intrinsic_subtypes/result/onco_result_only_shared_1p_082119.Rdata")
onco_result <- rbind(onco_result_shared_1p,onco_result_only_shared_1p)
load("/data/zhangh24/breast_cancer_data_analysis/whole_genome_age/ONCO/intrinsic_subtypes/result/onco_result_shared_1p_sex.Rdata")
onco_result_sex = onco_result_shared_1p
colnames(onco_result_sex) = colnames(onco_result)
onco_result = rbind(onco_result,onco_result_sex)
#load meta-analysis result
load("/data/zhangh24/breast_cancer_data_analysis/whole_genome_age/ICOG/Intrinsic_subtypes/result/meta_result_shared_1p_082119.Rdata")
meta_result = meta_result_shared_1p
load("/data/zhangh24/breast_cancer_data_analysis/whole_genome_age/ICOG/Intrinsic_subtypes/result/meta_result_shared_1p_final.Rdata")
meta_result_sex = meta_result_shared_1p
colnames(meta_result_sex) = colnames(meta_result)
meta_result = rbind(meta_result,meta_result_sex)
#load FTOP result
load("/data/zhangh24/breast_cancer_data_analysis/whole_genome_age/ICOG/Intrinsic_subtypes/result/meta_result_shared_1p_final.Rdata")
library(dplyr)
ftop = meta_result_shared_1p %>% 
  rename(FTOP_p_value = p.value) %>% 
  select(var_name,FTOP_p_value)
#load MTOP
load("/data/zhangh24/breast_cancer_data_analysis/whole_genome_age/ICOG/ERPRHER2GRADE_casecase/result/meta_result_shared_1p_final.Rdata")
mtop = meta_result_shared_1p %>% 
  rename(MTOP_p_value = p.value) %>% 
  select(var_name,MTOP_p_value)



#clean icog_result
icog_result_clean = icog_result[,c(11:15,c(1,7,13,19,25)+15,42)]
icog_result_select = icog_result_clean %>% 
  mutate(Luminal_A_sd = sqrt(X1.1),
         Luminal_B_sd = sqrt(X7),
         Luminal_B_HER2Neg_sd = sqrt(X13),
         HER2_Enriched_sd = sqrt(X19),
         Triple_Neg_sd = sqrt(X25)) %>% 
  rename(luminal_A_log_or = X1,
         Luminal_B_log_or = X2,
         Luminal_B_HER2Neg_log_or = X3,
         HER2_Enriched_log_or = X4,
         Triple_Neg_log_or = X5) %>% 
    select(var_name,luminal_A_log_or,Luminal_A_sd,
           Luminal_B_log_or,Luminal_B_sd,
           Luminal_B_HER2Neg_log_or,Luminal_B_HER2Neg_sd,
           HER2_Enriched_log_or,HER2_Enriched_sd,
           Triple_Neg_log_or,Triple_Neg_sd)
#clean oncoarray result 
onco_result_clean = onco_result[,c(11:15,c(1,7,13,19,25)+15,42)]
onco_result_select = onco_result_clean %>% 
  mutate(Luminal_A_sd = sqrt(X1.1),
         Luminal_B_sd = sqrt(X7),
         Luminal_B_HER2Neg_sd = sqrt(X13),
         HER2_Enriched_sd = sqrt(X19),
         Triple_Neg_sd = sqrt(X25)) %>% 
  rename(luminal_A_log_or = X1,
         Luminal_B_log_or = X2,
         Luminal_B_HER2Neg_log_or = X3,
         HER2_Enriched_log_or = X4,
         Triple_Neg_log_or = X5) %>% 
  select(var_name,luminal_A_log_or,Luminal_A_sd,
         Luminal_B_log_or,Luminal_B_sd,
         Luminal_B_HER2Neg_log_or,Luminal_B_HER2Neg_sd,
         HER2_Enriched_log_or,HER2_Enriched_sd,
         Triple_Neg_log_or,Triple_Neg_sd)

colnames(icog_result_select)[c(2:11)] = 
  paste0(colnames(icog_result_select)[c(2:11)],"_iCOGS")
colnames(onco_result_select)[c(2:11)] = paste0(colnames(onco_result_select)[c(2:11)],"_ONCO")
#join icogs and oncoarray result
icog_onco_select = full_join(icog_result_select,
                        onco_result_select,by="var_name")
#clean meta-analysis result

meta_result_clean = meta_result[,c(
                                   16:20,c(1,7,13,19,25)+20,12)]
colnames(meta_result_clean)[c(1:10)] <- paste0("X",colnames(meta_result_clean)[c(1:10)])
meta_result_select = meta_result_clean %>% 
  mutate(Luminal_A_sd = sqrt(X1.1),
         Luminal_B_sd = sqrt(X7),
         Luminal_B_HER2Neg_sd = sqrt(X13),
         HER2_Enriched_sd = sqrt(X19),
         Triple_Neg_sd = sqrt(X25)) %>% 
  rename(luminal_A_log_or = X1,
         Luminal_B_log_or = X2,
         Luminal_B_HER2Neg_log_or = X3,
         HER2_Enriched_log_or = X4,
         Triple_Neg_log_or = X5) %>% 
  select(var_name,luminal_A_log_or,Luminal_A_sd,
         Luminal_B_log_or,Luminal_B_sd,
         Luminal_B_HER2Neg_log_or,Luminal_B_HER2Neg_sd,
         HER2_Enriched_log_or,HER2_Enriched_sd,
         Triple_Neg_log_or,Triple_Neg_sd)
colnames(meta_result_select)[c(2:11)] = paste0(colnames(meta_result_select)[c(2:11)],"_meta")

#join meta analysis and ftop and mtop
meta_result_select = full_join(meta_result_select,ftop,by="var_name") %>% 
  rename(FTOP_p_value_meta = FTOP_p_value)
meta_result_select = full_join(meta_result_select,mtop,by="var_name") %>% 
  rename(MTOP_p_value_meta = MTOP_p_value)
#join icog onco meta analysis result
icog_onco_meta_select = full_join(icog_onco_select,
                                  meta_result_select,
                                  by="var_name")
#join the overall analysis information, with icog_onco_meta_select
icog_onco_meta_all = left_join(icog_onco_meta_select,
                               all.snp,by="var_name")

icog_onco_meta_final = icog_onco_meta_all[,c(1,45:53,2:11,61:69,12:21,77:78,22:33)]
write.table(icog_onco_meta_final,file = "")