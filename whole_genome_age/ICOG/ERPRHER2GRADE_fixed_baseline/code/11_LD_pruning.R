#Remove all the SNPs with LD2 > 0.1 with 210 known SNPs
#we will use the control group from oncoarray data to calculate the R2 with known SNPs since OncoArray is larger
setwd("/data/zhangh24/breast_cancer_data_analysis/")
library(data.table)
#set to your direcotry
data2 <- fread("./data/210_known_snp_onco_genotype_data.csv",header=T)
data2 <- as.data.frame(data2)
#take the 210 known snps value
known_snp_genotype <- data2[,c(19:228)]
#find the control group
idx.control <- which(data2$Behaviour1==0)
#find the control value
known_snp_genotype_control <- known_snp_genotype[idx.control,]
#load the merged extracted SNPs results (504 SNPs) from your directory
load("/data/zhangh24/breast_cancer_data_analysis/whole_genome_age/ONCO/ERPRHER2GRADE_fixed_baseline/result/extract_result_shared.Rdata")
# extract.list <- extract.list[-c(1700,1701,1702,1703,1705,1706,1707),]
# save(extract.list,file="/data/zhangh24/breast_cancer_data_analysis/whole_genome_age/ICOG/ERPRHER2_fixed/result/extract_list.Rdata")
#try <- duplicated(extract.list)
#take the SNP ID
extract_snp_id = extract.result[[1]]
#take the SNPs value
extract.result.onco <- extract.result[[2]]
#take the control group
extract.result.onco.control <- extract.result.onco[idx.control,]

#this is conceptual code;
#you need to modify it based on your setting
remove_idx = NULL
for(k in 1:ncol(extract.result.onco.control)){
  temp_snp = extract.result.onco.control[,k]
  LD = cor(temp_snp, known_snp_genotype_control)^2
  if(max(LD)>0.1){
    remove_idx = c(remove_idx,k)
  }
}
#load data_clean_select_sig from investigate_result.r
#remove the SNPs from the data_clean_select_sig
remove_id = extract_snp_id[remove_idx]
idx = which(data_clean_select_sig$SNP.ONCO%in%remove_id)
data_clean_significant = data_clean_select_sig[-idx,]

#data_clean significant are the SNPs removing +-500kb of known SNPs and LD >0.1 with any known SNPs
load("/data/NC_BW/HZ_SF/data_clean_significant.rdata")

#only keep top signal in one region
FilterSignal <- function(CHR, position, p_value){
  picked_ind = NULL
  #create temporary data for removing
  CHR_temp = CHR
  position_temp = position
  p_value_temp = p_value
  while(length(CHR_temp)!=0){
    #find the snp with smallest p-value
    idx = which.min(p_value_temp)
    #remove all SNPs within +-500kb region
    CHR_pick = CHR_temp[idx]
    position_pick = position_temp[idx]
    #find the original index in the data
    jdx = which(CHR==CHR_pick&position==position_pick)
    picked_ind = c(picked_ind,jdx)
    remove_idx = which(CHR_temp==CHR_pick&
                         (position_temp<=position_pick+500000)&
                         (position_temp>=position_pick-500000))
    
    CHR_temp = CHR_temp[-remove_idx]
    position_temp = position_temp[-remove_idx]
    p_value_temp = p_value_temp[-remove_idx] 
  }
  return(picked_ind)
}


CHR = data_clean_significant$CHR
position = data_clean_significant$position
p_value = data_clean_significant$p.acat

picked_ind = FilterSignal(CHR, position, p_value)

#find the top results in each region
data_clean_select_sig = data_clean_significant[picked_ind,]
dim(data_clean_select_sig)
library(data.table)
#check whether the SNPs are within +-2Mb within any known SNPs
known_snp = fread("/data/zhangh24/breast_cancer_data_analysis/data/210_known_discovery_snp_paper_order.csv")
# #load overall GWAS result to find CHR and position information
# data = fread("/data/zhangh24/breast_cancer_data_analysis/discovery_SNP/prepare_summary_level_statistics/result/icogs_onco_gwas_meta_overall_breast_cancer_summary_level_statistics.txt")
# library(dplyr)
# data = data %>% 
#   mutate(chr.pos = paste0(chr.Onco,":",Position.Onco)) %>% 
#   select(var_name, chr.pos, SNP.iCOGs, chr.iCOGs, Position.iCOGs,
#          SNP.Onco, chr.Onco, Position.Onco)
# known_snp = known_snp %>% 
#   mutate(chr.pos = paste0(CHR,":",position))
# known_snp_update = left_join(known_snp,data,by="chr.pos")
#   
# idx <- which(duplicated(known_snp_update$chr.pos))
# jdx <- which(duplicated(known_snp$Best.published.SNP))

n.known = nrow(known_snp)
idx = NULL
for(i in 1:n.known){
  CHRi = known_snp$CHR[i]
  lb = known_snp$position[i]-2*10^6
  ub = known_snp$position[i]+2*10^6
  temp =  which((data_clean_select_sig$CHR==CHRi) & (data_clean_select_sig$position>=lb) & (data_clean_select_sig$position<=ub))
  idx = c(idx,temp)
  cat(i,length(idx),"\n")
}

idx = unique(idx)
length(idx)  


#10 SNPs within +-2MB
data_clean_select_1 = data_clean_select_sig[idx,]
#6 SNPs with clean signals data_clean_select_2
data_clean_select_2 = data_clean_select_sig[-idx,]
#find corresponding known SNP within +-2MB
for(k in 1:nrow(data_clean_select_1)){
  chr.temp = data_clean_select_1[k,"CHR"]
  position.temp = data_clean_select_1[k,"position"]
  idx <- which(known_snp$CHR==chr.temp&
                 (position.temp>= known_snp$position-2*10^6)&
                 position.temp<= known_snp$position+2*10^6)
  print(length(idx))
  known_snp[idx,]
}




