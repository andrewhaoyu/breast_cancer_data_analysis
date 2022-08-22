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

#load meta-analysis result for sex chromosome
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
write.table(icog_onco_meta_final,file = "./discovery_SNP/prepare_summary_level_statistics/result/icogs_onco_meta_intrinsic_subtypes_summary_level_statistics.txt",col.names = T,row.names = F,quote=F)

#load cimba bcac meta-analysis result

meta_result <- as.data.frame(fread("./data/brca1_bcac_tn_meta.txt",header=T))
colnames(meta_result)[c(2,3)] <- c("eff_allele","ref_allele")
meta_result <- meta_result[,-c(16,17)]
library(tidyr)
meta_result_new = meta_result %>% 
  separate(MarkerName,
           c("CHR","position",
             "Allele1","Allele2"),sep="_",remove=F)
idx <- which(meta_result_new$Allele2!=
               toupper(meta_result_new$eff_allele))
meta_result_new$Effect[idx] <- -meta_result_new$Effect[idx]
meta_result_new$Freq1[idx] <- 1-meta_result_new$Freq1[idx]
colnames(meta_result_new)[8] <- c("eff_freq")
meta_result_clean = meta_result_new[,c(1:5,8,9,12:19)]
jdx <- which(meta_result_clean$MarkerName=="17_7571752_T_G")
write.table(meta_result_clean,file ="./discovery_SNP/prepare_summary_level_statistics/result/CIMBA_BRCA1_BCAC_TN_meta_summary_level_statistics.txt",quote=F,row.names = F,col.names = F)

meta_result_clean <- as.data.frame(fread("./discovery_SNP/prepare_summary_level_statistics/result/CIMBA_BRCA1_BCAC_TN_meta_summary_level_statistics.txt"))
#load the 1kg data

data.list = list()
for(i in 1:22){
  leg <- as.data.frame(fread(paste0("/data/zhangh24/KG.impute2/1000GP_Phase3/1000GP_Phase3_chr",i,".legend"),header=T))
  n <- nrow(leg)
  CHR = rep(i,n)
  temp.data = cbind(leg[,1:4],CHR)
  data.list[[i]] = temp.data
}

kg <- rbindlist(data.list)
library(tidyr)
kg = kg %>% separate(id,
              c("rsid","position_new",
                "Allele1_new","Allele2_new"),sep=":",remove=F) %>% 
  select(rsid,CHR,position,a0,a1) %>% 
  mutate(chr.pos = paste0(CHR,"_",position))

idx <- which(kg$chr.pos=="5_7352452")

kg_new = kg %>% select(rsid,chr.pos)
#keep 1kg unique
kg_new <- distinct(kg_new,chr.pos,.keep_all = TRUE)

meta_result_clean = meta_result_clean %>% 
  mutate(chr.pos = paste0(CHR,"_",position))

meta_result_new = left_join(meta_result_clean,kg_new,by="chr.pos")

idx <- which(meta_result_new$chr.pos=="5_7352452")
meta_result_new[idx,]
  meta_result_new_id = meta_result_new %>% 
  select(MarkerName,rsid)
write.table(meta_result_new_id,file ="./discovery_SNP/prepare_summary_level_statistics/result/CIMBA_BRCA1_BCAC_TN_meta_summary_level_id_information.txt",quote=F,row.names = F,col.names = F)

library(data.table)
data <- as.data.frame(fread("./discovery_SNP/prepare_summary_level_statistics/result/icogs_onco_meta_intrinsic_subtypes_summary_level_statistics.txt"))
idx <- which(data$var_name=="5_162902516_T_C" )
data[idx,]

data <- data[idx,c(1,2,3,4,5,6,7,8,9,42:51)]
p <- rep(0,5)
for(k in 1:5){
  z = data[,9+2*k-1]/data[,9+2*k]
  p[k] <- 2*pnorm(-abs(z))
}
Luminal_A_p <- p[1]
Luminal_B_p <- p[2]
Luminal_B_HER2Neg_p <- p[3]
HER2_Enriched_p <- p[4]
Triple_Neg_p <- p[5]
data_new <- cbind(data[,c(1:11)],
              Luminal_A_p,
              data[,c(12:13)],
              Luminal_B_p,
              data[,c(14:15)],
              Luminal_B_HER2Neg_p,
              data[,c(16:17)],
              HER2_Enriched_p,
              data[,c(18:19)],
              Triple_Neg_p)
write.csv(data_new,file = "./discovery_SNP/prepare_summary_level_statistics/result/requested_snp_by_collaborator.csv")








#load the results locally
setwd("/Users/zhangh24/GoogleDrive/breast_cancer_data_analysis/discovery_SNP/prepare_summary_level_statistics/result")
library(data.table)
overall_result <- as.data.frame(fread("icogs_onco_gwas_meta_overall_breast_cancer_summary_level_statistics.txt"))
#overall_result <- overall_result[,-c(8:12)]
write.table(overall_result,file = "icogs_onco_gwas_meta_overall_breast_cancer_summary_level_statistics.txt",row.names = F,col.names = T,quote=F)


setwd("/Users/zhangh24/GoogleDrive/breast_cancer_data_analysis/discovery_SNP/prepare_summary_level_statistics/result")
library(data.table)
subtypes_result <- as.data.frame(fread("icogs_onco_meta_intrinsic_subtypes_summary_level_statistics.txt"))
names <- colnames(subtypes_result)
new.names <- gsub("_sd_","_se_",names)
new.names <- gsub("luminal","Luminal",new.names)
colnames(subtypes_result) = new.names
#overall_result <- overall_result[,-c(8:12)]
write.table(subtypes_result,file = "icogs_onco_meta_intrinsic_subtypes_summary_level_statistics.txt",row.names = F,col.names = T,quote=F)
