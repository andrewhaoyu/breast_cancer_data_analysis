#-------------------------------------------------------------------
# Update Date: 01/21/2020
# Create Date: 11/22/2018
# Goal: prepare different different prs files for plink
# Author: Haoyu Zhang
#-------------------------------------------------------------------
#---------------------------------------#---------------------------------------
# setwd('/data/zhangh24/breast_cancer_data_analysis/risk_prediction')
# 
# #load("./EB_whole_genome/result/whole_gonome.rdata")
# load("./intrinsic_subtypes_whole_genome/ICOG/result/whole_gonome.rdata")
# idx <- which(whole_genome$SNP.ONCO=="rs56097627:110198129:CAAA:C")
# whole_genome[idx,c(14,19:23)]
# #LD.assoc <- read.table("/spin1/users/zhangh24/BCAC/impute_plink_onco/LD_assoc",header=T)
# #idx <- which(LD.assoc$SNP=="chr1_121280613_A_G")
# #LD.assoc[idx,]
# #######prs file need A1 to be effect_allele
# #######we need to reverse all the log odds ratio since we put A2 as effect allele
# #whole_genome = whole_genome %>% mutate(p.min = pmin(p.value,FTOP_result))
# #head(whole_genome)
# #save(whole_genome,file = "./EB_whole_genome/result/whole_gonome.rdata")
# #LD clumping based on the min-p value of two-stage model and standard analysis
# library(data.table)
# 
# library(dplyr)
# clump.snp <- as.data.frame(fread("/data/zhangh24/BCAC/impute_plink_onco/clump_snp",header=F))
# clump.snp <- clump.snp %>% filter(clump.snp!="SNP"&
#                                     clump.snp!="")
# dim(clump.snp)
# dim(whole_genome)
# colnames(clump.snp) <- c("SNP.ONCO")
# #check duplicated
# idx <- which(duplicated(clump.snp$SNP.ONCO)==T)
# length(idx)
# whole_genome_clump <- left_join(clump.snp,whole_genome)
# #load in 313 SNPs
# setwd('/data/zhangh24/breast_cancer_data_analysis')
# snp <- read.csv("./data/Nasim_313_SNPs_infor.csv",header=T)
# library(dplyr)
# snp = snp %>% mutate(chr.pos = paste0(CHR,":",position))
# #find out all the snps within 313 SNPs
# idx <- which(whole_genome_clump$var_name%in%
#                snp$SNPa)
# #put the 313 SNPs p-value as 0 into the list so it can be selected
# whole_genome_clump$p.min[idx] = 0
# #idx <- which(whole_genome_clump$SNP.ONCO=="chr1_121280613_A_G")
# #whole_genome_clump[idx,]
# #save(whole_genome_clump,file = "./EB_whole_genome/result/whole_genome_clump.rdata")
# #No need to rerun the previous code again
# 
# #create prs files based on different p-threshold
# setwd('/data/zhangh24/breast_cancer_data_analysis/risk_prediction')
# library(dplyr)
# library(data.table)
# #load("./EB_whole_genome/result/whole_genome_clump.rdata")
# dim(whole_genome_clump)
# #method <- c("standard","two-stage","eb")
# pthres <- c(5E-08,1E-07,5E-07,1E-06,5E-06,1E-05,5E-05,
#             1E-04,5E-04,1E-03,5E-03,1E-02)
# n.pthres <- length(pthres)
# 
# #create the prs file for two-stage and eb
# subtypes <- c("Luminal_A",
#               "Luminal_B",
#               "Luminal_B_HER2Neg",
#               "HER2_Enriched",
#               "TN")
# select.names <- subtypes
# score <- whole_genome_clump %>%  select(select.names)
# whole_genome_clump_new <- whole_genome_clump %>% mutate(SNP=SNP.ONCO) %>% select(SNP,effect_allele,p.min) %>% 
#   cbind(score)
# idx <- which(whole_genome_clump_new$SNP=="chr1_100880328_A_T")
# 
# 
# 
# idx <- which(whole_genome_clump$CHR==1)
# head(whole_genome_clump[idx,])
# 
# for(i in 1:n.pthres){
#   for(j in 1:length(select.names)){
#     prs <-  whole_genome_clump_new %>%
#       filter(p.min<=pthres[i]) %>% 
#       select(SNP,effect_allele,select.names[j])
#     colnames(prs) <- c("SNP","effect_allele","beta")
#     dim(prs)
#     write.table(prs,file = paste0("/data/zhangh24/breast_cancer_data_analysis/risk_prediction/subtypes_prs/result/",select.names[j],"_prs_pvaluecut_",i,".file"),row.names=F,col.names=T,quote=F)
#     
#   }
#   
# }
# 
# 
# 
# 
# 
# 
# 
# 
# idx <- which(prs$SNP=="chr1_100880328_A_T")
# prs[idx,]
# #
# i <- 8
# j <- 1
# 
# prs.temp <-  whole_genome_clump_new %>%
#   filter(p.min<=pthres[i]) %>%
#   select(SNP,effect_allele,select.names[j])
# colnames(prs.temp) <- c("SNP","effect_allele","beta")
# dim(prs.temp)
# n.snp <- 1701
# prs.temp <- prs.temp[c(1:n.snp),,drop=F]
# #prs.temp <- prs.temp[c(n.snp),,drop=F]
# write.table(prs.temp,file = paste0("/data/zhangh24/breast_cancer_data_analysis/risk_prediction/subtypes_prs/result/test.file"),row.names=F,col.names=T,quote=F)
# beta = prs[,3]
# 
# knwon_snp_new_update <- known_snp_new[complete.cases(known_snp_new[,2]),]
# prs.temp <-  knwon_snp_new_update %>%
#   select(SNP,effect_allele,select.names[j])
# colnames(prs.temp) <- c("SNP","effect_allele","beta")
# dim(prs.temp)
# n.snp <- 1
# prs.temp <- prs.temp[c(1:n.snp),,drop=F]
# #prs.temp <- prs.temp[c(n.snp),,drop=F]
# write.table(prs.temp,file = paste0("/data/zhangh24/breast_cancer_data_analysis/risk_prediction/subtypes_prs/result/test.file"),row.names=F,col.names=T,quote=F)
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# dosage.assoc <- read.table("/data/zhangh24/breast_cancer_data_analysis/risk_prediction/subtypes_prs/result/dosage_test_out.assoc.dosage",header=T)
# 
# idx <- which(dosage.assoc$SNP=="rs370540207:121450795:T:G")
# # 
# # 
# # 
# # 
# # 
# # 
# # 
# # 
# # 
# # 
# 
# 
# 
# 
# #create the file for standard logisic regression
# n.pthres <- length(pthres)
# whole_genome_clump_new <- whole_genome_clump %>% mutate(SNP=SNP.ONCO,
#                               beta = -score
# )
# for(i in 1:n.pthres){
#   prs <-  whole_genome_clump_new %>%
#     filter(p.min<=pthres[i]) %>% 
#     select(SNP,referece_allele,beta)
#   colnames(prs) <- c("SNP","effect_allele","beta")
#   dim(prs)
#   write.table(prs,file = paste0("/data/zhangh24/BCAC/prs_file/standard_prs_",i,".file"),row.names=F,col.names=T,quote=F)
#   
# }
# 
# #create the prs file for two-stage and eb
# subtypes <- c("Luminal_A",
#               "Luminal_B",
#               "Luminal_B_HER2Neg",
#               "HER2_Enriched",
#               "TN")
# select.names <- c(subtypes,paste0("eb_",subtypes))
# 
# score <- whole_genome_clump %>%  select(select.names)
# #reverse the odds ratio, since dosage use the first allele as reference
# score <- -score
# whole_genome_clump_new <- whole_genome_clump %>% mutate(SNP=SNP.ONCO) %>% select(SNP,referece_allele,p.min) %>% 
#   cbind(score)
# 
# for(i in 1:n.pthres){
#   for(j in 1:length(select.names)){
#     prs <-  whole_genome_clump_new %>%
#       filter(p.min<=pthres[i]) %>% 
#       select(SNP,referece_allele,select.names[j])
#     colnames(prs) <- c("SNP","effect_allele","beta")
#     dim(prs)
#     write.table(prs,file = paste0("/data/zhangh24/BCAC/prs_file/",select.names[j],"_prs_",i,".file"),row.names=F,col.names=T,quote=F)
#     
#   }
#   
# }
# 








#
#-------------------------------------------------------------------
# Update Date: 01/23/2020
# Create Date: 11/22/2018
# Goal: prepare different different prs files for plink
# Author: Haoyu Zhang
#-------------------------------------------------------------------
#---------------------------------------#---------------------------------------
setwd('/data/zhangh24/breast_cancer_data_analysis/risk_prediction')

#load("./EB_whole_genome/result/whole_gonome.rdata")
load("./intrinsic_subtypes_whole_genome/ICOG/result/whole_genome_threeadd.rdata")
idx <- which(whole_genome$SNP.ONCO=="rs56097627:110198129:CAAA:C")
whole_genome[idx,c(14,19:23)]
#LD.assoc <- read.table("/spin1/users/zhangh24/BCAC/impute_plink_onco/LD_assoc",header=T)
#idx <- which(LD.assoc$SNP=="chr1_121280613_A_G")
#LD.assoc[idx,]
#######prs file need A1 to be effect_allele
#######we need to reverse all the log odds ratio since we put A2 as effect allele
#whole_genome = whole_genome %>% mutate(p.min = pmin(p.value,FTOP_result))
#head(whole_genome)
#save(whole_genome,file = "./EB_whole_genome/result/whole_gonome.rdata")
#LD clumping based on the min-p value of two-stage model and standard analysis
library(data.table)

library(dplyr)
clump.snp <- as.data.frame(fread("/data/zhangh24/BCAC/impute_plink_onco/clump_snp_121019",header=F))
clump.snp <- clump.snp %>% filter(clump.snp!="SNP"&
                                    clump.snp!="")
dim(clump.snp)
dim(whole_genome)
colnames(clump.snp) <- c("SNP.ONCO")
#check duplicated
idx <- which(duplicated(clump.snp$SNP.ONCO)==T)
length(idx)
whole_genome_clump <- left_join(clump.snp,whole_genome)




setwd('/data/zhangh24/breast_cancer_data_analysis/')
load("./data/Nasim_313SNPs_complete_information.Rdata")
#check whether the 313 SNPs are in the list
head(snp.new)
idx <- which(snp.new$var_name%in%
               whole_genome_clump$var_name==F)
length(idx)

#idx <- which(whole_genome_clump$SNP.ONCO=="chr1_121280613_A_G")
#whole_genome_clump[idx,]
#save(whole_genome_clump,file = "./EB_whole_genome/result/whole_genome_clump.rdata")
#No need to rerun the previous code again

#create prs files based on different p-threshold
setwd('/data/zhangh24/breast_cancer_data_analysis/risk_prediction')
library(dplyr)
library(data.table)
#load("./EB_whole_genome/result/whole_genome_clump.rdata")
dim(whole_genome_clump)
#method <- c("standard","two-stage","eb")
 # pthres <- c(1E-30,5E-08,1E-07,5E-07,1E-06,5E-06,1E-05,5E-05,
 #             1E-04,5E-04,1E-03,5E-03)
 pthres <- c(1E-30,1E-10,5E-08,1E-07,5E-07,1E-06,5E-06,1E-05,5E-05,1E-04,1E-03,1E-02)
 
 n.pthres <- length(pthres)

#create the prs file for two-stage and eb
subtypes <- c("Luminal_A",
              "Luminal_B",
              "Luminal_B_HER2Neg",
              "HER2_Enriched",
              "TN")
select.names <- subtypes
score <- whole_genome_clump %>%  select(select.names)
# whole_genome_clump_new <- whole_genome_clump %>% mutate(SNP=SNP.ONCO) %>% select(SNP,effect_allele,p.min) %>% 
#   cbind(score)
whole_genome_clump_new <- whole_genome_clump %>% mutate(SNP=SNP.ONCO) %>% 
select(SNP,effect_allele,p.min,CHR,position) %>% 
  cbind(score)

idx <- which(whole_genome_clump$SNP.ONCO=="chr22_28195386_A_C")
whole_genome_clump[idx,]


idx <- which(whole_genome_clump$CHR==1)
head(whole_genome_clump[idx,])

for(i in 1:n.pthres){
  #keep the 313 SNPs, then for all other SNPs, make sure only the top one within 500kb are kept
  prs <-  whole_genome_clump_new %>%
    filter(p.min<=pthres[i]) 
  idx <- which(prs$SNP%in%snp.new$SNP.ONCO==F)
  if(length(idx)!=0){
    prs.no.nasim <- prs[idx,,drop=F]
    #remove SNP  rs76858104 which is in LD with Nasim SNP rs11249433
    #remove SNP that are not indepedent with prevoius known SNPs
    qdx <- which(prs.no.nasim$p.min<=8.2E-12)
    prs.no.nasim = prs.no.nasim[-qdx,,drop=F]
    snp.keep <- NULL
    for(k in 1:nrow(prs.no.nasim)){
      #check whether there are any SNPs that are within 500kb of the top SNPs
      jdx <- which(prs.no.nasim$CHR==prs.no.nasim$CHR[k]&
                     prs.no.nasim$position>=prs.no.nasim$position[k]-500000&
                     prs.no.nasim$position<=prs.no.nasim$position[k]+500000&           prs.no.nasim$SNP!=prs.no.nasim$SNP[k])
      
        if(sum(prs.no.nasim$p.min[k]<prs.no.nasim$p.min[jdx])==length(jdx)){
          snp.keep = c(snp.keep,prs.no.nasim$SNP[k])
        }
      
      
    }
    #only keep the 313 SNPs and snps that have no nearby snps in +-500kb
    kdx <- which(prs$SNP%in%snp.new$SNP.ONCO|
                   prs$SNP%in%snp.keep)
    prs = prs[kdx,,drop=F]
  }
  
  for(j in 1:length(select.names)){
    
    prs.new <-  prs %>%
      select(SNP,effect_allele,select.names[j])
    colnames(prs.new) <- c("SNP","effect_allele","beta")
      
    
    # prs <-  whole_genome_clump_new %>%
    #   filter(p.min<=pthres[i]) %>% 
    #   select(SNP,effect_allele,select.names[j])
    # colnames(prs) <- c("SNP","effect_allele","beta")
    dim(prs.new)
    write.table(prs.new,file = paste0("/data/zhangh24/breast_cancer_data_analysis/risk_prediction/subtypes_prs/result/",select.names[j],"_prs_pvaluecut_",i,"_121019.file"),row.names=F,col.names=T,quote=F)
    
  }
  
}

idx <- which(prs$SNP%in%
               snp.new$SNP.ONCO==F)

prs[idx,]







