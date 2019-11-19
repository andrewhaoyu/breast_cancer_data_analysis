#-------------------------------------------------------------------
# Update Date: 11/25/2018
# Create Date: 11/22/2018
# Goal: prepare different different prs files for plink
# Author: Haoyu Zhang
#-------------------------------------------------------------------
#---------------------------------------#---------------------------------------
setwd('/data/zhangh24/breast_cancer_data_analysis/risk_prediction')

#load("./EB_whole_genome/result/whole_gonome.rdata")
load("./intrinsic_subtypes_whole_genome/ICOG/result/whole_gonome.rdata")
#idx <- which(whole_genome$SNP.ONCO=="chr1_121280613_A_G")
#whole_genome[idx,]
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
clump.snp <- as.data.frame(fread("/data/zhangh24/BCAC/impute_plink_onco/clump_snp",header=F))
clump.snp <- clump.snp %>% filter(clump.snp!="SNP"&
                                    clump.snp!="")
dim(clump.snp)
dim(whole_genome)
colnames(clump.snp) <- c("SNP.ONCO")
#check duplicated
idx <- which(duplicated(clump.snp$SNP.ONCO)==T)
length(idx)
whole_genome_clump <- left_join(clump.snp,whole_genome)
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
pthres <- c(5E-08,1E-07,5E-07,1E-06,5E-06,1E-05,5E-05,
            1E-04,5E-04,1E-03,5E-03,1E-02)
n.pthres <- length(pthres)


subtypes <- c("stan_logodds")
select.names <- subtypes
score <- whole_genome_clump %>%  select(select.names)
whole_genome_clump_new <- whole_genome_clump %>% mutate(SNP=SNP.ONCO) %>% select(SNP,effect_allele,p.min) %>% 
  cbind(score)

for(i in 1:n.pthres){
  for(j in 1:length(select.names)){
    prs <-  whole_genome_clump_new %>%
      filter(p.min<=pthres[i]) %>% 
      select(SNP,effect_allele,select.names[j])
    colnames(prs) <- c("SNP","effect_allele","beta")
    dim(prs)
    write.table(prs,file = paste0("/data/zhangh24/breast_cancer_data_analysis/risk_prediction/subtypes_prs/result/",select.names[j],"_prs_pvaluecut_",i,".file"),row.names=F,col.names=T,quote=F)
    
  }
  
}



#
i <- 8
j <- 1

prs.temp <-  whole_genome_clump_new %>%
  filter(p.min<=pthres[i]) %>%
  select(SNP,effect_allele,select.names[j])
colnames(prs.temp) <- c("SNP","effect_allele","beta")
dim(prs.temp)
n.snp <- 1701
prs.temp <- prs.temp[c(1:n.snp),,drop=F]
#prs.temp <- prs.temp[c(n.snp),,drop=F]
write.table(prs.temp,file = paste0("/data/zhangh24/breast_cancer_data_analysis/risk_prediction/subtypes_prs/result/test.file"),row.names=F,col.names=T,quote=F)
beta = prs[,3]

knwon_snp_new_update <- known_snp_new[complete.cases(known_snp_new[,2]),]
prs.temp <-  knwon_snp_new_update %>%
  select(SNP,effect_allele,select.names[j])
colnames(prs.temp) <- c("SNP","effect_allele","beta")
dim(prs.temp)
n.snp <- 1
prs.temp <- prs.temp[c(1:n.snp),,drop=F]
#prs.temp <- prs.temp[c(n.snp),,drop=F]
write.table(prs.temp,file = paste0("/data/zhangh24/breast_cancer_data_analysis/risk_prediction/subtypes_prs/result/test.file"),row.names=F,col.names=T,quote=F)














dosage.assoc <- read.table("/data/zhangh24/breast_cancer_data_analysis/risk_prediction/subtypes_prs/result/dosage_test_out.assoc.dosage",header=T)

idx <- which(dosage.assoc$SNP=="rs370540207:121450795:T:G")
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




#create the file for standard logisic regression
n.pthres <- length(pthres)
whole_genome_clump_new <- whole_genome_clump %>% mutate(SNP=SNP.ONCO,
                              beta = -score
)
for(i in 1:n.pthres){
  prs <-  whole_genome_clump_new %>%
    filter(p.min<=pthres[i]) %>% 
    select(SNP,referece_allele,beta)
  colnames(prs) <- c("SNP","effect_allele","beta")
  dim(prs)
  write.table(prs,file = paste0("/data/zhangh24/BCAC/prs_file/standard_prs_",i,".file"),row.names=F,col.names=T,quote=F)
  
}

#create the prs file for two-stage and eb
subtypes <- c("Luminal_A",
              "Luminal_B",
              "Luminal_B_HER2Neg",
              "HER2_Enriched",
              "TN")
select.names <- c(subtypes,paste0("eb_",subtypes))

score <- whole_genome_clump %>%  select(select.names)
#reverse the odds ratio, since dosage use the first allele as reference
score <- -score
whole_genome_clump_new <- whole_genome_clump %>% mutate(SNP=SNP.ONCO) %>% select(SNP,referece_allele,p.min) %>% 
  cbind(score)

for(i in 1:n.pthres){
  for(j in 1:length(select.names)){
    prs <-  whole_genome_clump_new %>%
      filter(p.min<=pthres[i]) %>% 
      select(SNP,referece_allele,select.names[j])
    colnames(prs) <- c("SNP","effect_allele","beta")
    dim(prs)
    write.table(prs,file = paste0("/data/zhangh24/BCAC/prs_file/",select.names[j],"_prs_",i,".file"),row.names=F,col.names=T,quote=F)
    
  }
  
}


