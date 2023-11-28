#prepare the known SNP data for icogs and oncoarray
#load the orignal extract snps
library(data.table)
library(dplyr)
discovery_snp = fread("/data/zhangh24/breast_cancer_data_analysis/data/discovery_snp_summary_new.csv")
#Originally, there are 35 discovery SNPs
known_snp = fread("/data/zhangh24/breast_cancer_data_analysis/data/210_known_discovery_snp_paper_order.csv")
known_snp = known_snp %>% 
  mutate(chr.pos = paste0(CHR,":",position))
#210th SNP is the OncoArray only SNP from original 178 known SNP
discovery_snp_paper_order = known_snp[178:209,]
#there are three snps removed in the published paper due to conditional analysis p-value <1E-06
idx <- which(discovery_snp$Pos%in%known_snp$position==F)
remove_snp_infor = discovery_snp[idx,]
discovery_snp_clean = discovery_snp[-idx,]
#prepare idx_match to match the genotype data to the same order as 210_known_snp
idx_match = match(discovery_snp_paper_order$chr.pos, discovery_snp_clean$chr.pos)
discovery_snp_mathch = discovery_snp_clean[idx_match,]
#sig_snp data were prepared for PRS analysis
#the original 178 known SNPs are coded in rsid prepared by Bill
#the later 32 discovery SNPs are coded in SNP.ICOGs or SNP.Onco prepared by Haoyu
#we need to align the 32 discovery SNPs and rename column name with rsid
data1 = as.data.frame(fread("/data/zhangh24/breast_cancer_data_analysis/data/sig_snp_icog_prs.csv"))
#when three SNPs should be removed due to 
name1 = colnames(data1)
idx <- which(name1%in%remove_snp_infor$SNP.ICOGS)
data1_clean = data1[,-idx]
name1 = colnames(data1_clean)
#test whether the later 32 SNPs are the same order as discovery_snp_clean
all.equal(name1[195:226],discovery_snp_clean$SNP.ICOGS)
#test whether the matched discovery snp are the same as the required order
all.equal(discovery_snp_mathch$chr.pos, discovery_snp_paper_order$chr.pos)
required_name_order = c(name1[18:194],discovery_snp_mathch$SNP.ICOGS)
#use data1_clean to obtain the name for 210 discovery snp
#load the original dataset
data1_known = as.data.frame(fread("/data/zhangh24/breast_cancer_data_analysis/data/iCOGS_euro_v10_10232017.csv",header=T))
data1_discovery <- as.data.frame(fread("/data/zhangh24/breast_cancer_data_analysis/discovery_SNP/result/discovery_icog_data.csv",header=T))
data1_com = cbind(data1_known,
                  data1_discovery)
#find the SNP data and age
data1_clean_new = data1_com[,required_name_order]
colnames(data1_clean_new) = as.character(known_snp$Best.published.SNP[1:209])
library(dplyr)
data1_infor = data1_com %>% 
  select(SG_ID, Behaviour1, Morphologygroup1_derived, ER_status1, PR_status1,
         HER2_status1, Grade1, age, pc1, 
         pc2, pc3, pc4, pc5, 
         pc6, pc7, pc8, pc9, pc10)
data1_clean = cbind(data1_infor,data1_clean_new)


write.csv(data1_clean, file = "/data/zhangh24/breast_cancer_data_analysis/data/210_known_snp_icog_genotype_data.csv", 
          quote = F, row.names = F)
write.csv(data1_clean, file = "/data/NC_BW/210_known_snp_icog_genotype_data.csv", 
          quote = F, row.names = F)

setwd("/data/zhangh24/breast_cancer_data_analysis/")
data2 = as.data.frame(fread("/data/zhangh24/breast_cancer_data_analysis/data/sig_snp_onco_prs.csv"))
name2 = colnames(data2)
idx <- which(name2%in%remove_snp_infor$SNP.ONCO)
data2_clean = data2[,-idx]
colnames(data2_clean)[228] = "rs554219"
name2 = colnames(data2_clean)
#test whether the later 32 SNPs are the same order as discovery_snp_clean
all.equal(name2[195:226],discovery_snp_clean$SNP.ONCO)
#test whether the matched discovery snp are the same as the required order
all.equal(discovery_snp_mathch$chr.pos, discovery_snp_paper_order$chr.pos)
required_name_order = c(name2[18:194],discovery_snp_mathch$SNP.ONCO,name2[228])
#use data1_clean to obtain the name for 210 discovery snp
#load the original dataset
#use data1_clean to obtain the name for 210 discovery snp
#load the original dataset
data2_known <- fread("./data/Onco_euro_v10_10232017.csv",header=T)
data2_known <- as.data.frame(data2_known)
data2_discovery <- as.data.frame(fread("/data/zhangh24/breast_cancer_data_analysis/discovery_SNP/result/discovery_onco_data.csv",header=T))
data2_com = cbind(data2_known,
                  data2_discovery)
data2_clean_new = data2_com[,required_name_order]
colnames(data2_clean_new) =as.character(known_snp$Best.published.SNP[1:210])
data2_infor = data2_com %>% 
  select(Onc_ID, Behaviour1, Morphologygroup1_derived, 
         ER_status1, PR_status1,
         HER2_status1, Grade1, age, PC_1, 
         PC_2, PC_3, PC_4, PC_5, 
         PC_6, PC_7, PC_8, PC_9, PC_10) %>% 
  rename(pc1 = PC_1,
         pc2 = PC_2,
         pc3 = PC_3,
         pc4 = PC_4,
         pc5 = PC_5,
         pc6 = PC_6,
         pc7 = PC_7,
         pc8 = PC_8,
         pc9 = PC_9,
         pc10 = PC_10)
data2_clean = cbind(data2_infor,data2_clean_new)


write.csv(data2_clean, file = "/data/zhangh24/breast_cancer_data_analysis/data/210_known_snp_onco_genotype_data.csv", 
          quote = F, row.names = F)

write.csv(data2_clean, file = "/data/NC_BW/210_known_snp_onco_genotype_data.csv", 
          quote = F, row.names = F)
#idx <- which(name1%in%remove_snp_infor$SNP.ICOGS)