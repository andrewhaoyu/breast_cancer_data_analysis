#-------------------------------------------------------------------
# Update Date: 11/23/2018
# Create Date: 11/22/2018
# Goal: prepare different different prs files for plink
# Author: Haoyu Zhang
#-------------------------------------------------------------------
#---------------------------------------#---------------------------------------
setwd('/spin1/users/zhangh24/breast_cancer_data_analysis/risk_prediction')
load("./EB_whole_genome/result/whole_gonome.rdata")
#######prs file need A1 to be effect_allele
#######we need to reverse all the log odds ratio since we put A2 as effect allele
library(dplyr)
library(data.table)
#whole_genome = whole_genome %>% mutate(p.min = pmin(p.value,FTOP_result))
#head(whole_genome)
#save(whole_genome,file = "./EB_whole_genome/result/whole_gonome.rdata")
clump.snp <- as.data.frame(fread("/spin1/users/zhangh24/BCAC/impute_plink_onco/clump_snp",header=F))
colnames(clump.snp) <- "SNP.ONCO"
dim(clump.snp)
dim(whole_genome)
whole_genome_clump <- left_join(clump.snp,whole_genome)
save(whole_genome_clump,file = "./EB_whole_genome/result/whole_genome_clump.rdata")
#No need to rerun the previous code again

#create prs files based on different p-threshold
load("./EB_whole_genome/result/whole_genome_clump.rdata")
dim(whole_genome_clump)
method <- c("standard","two-stage","eb")
pthres <- c(5E-08,1E-07,5E-07,1E-06,5E-06,1E-05,5E-05,
            1E-04,5E-04,1E-03,5E-03,1E-02)
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
  write.table(prs,file = paste0("/spin1/users/zhangh24/BCAC/prs_file/standard_prs_",i,".file"),row.names=F,col.names=T,quote=F)
  
}
#create the prs file for two-stage and eb
subtypes <- c("Luminial_A",
              "Luminal_B",
              "Luminal_B_HER2Neg",
              "HER2_Enriched",
              "TN")
select.names <- c(subtypes,paste0("eb_",subtypes))

try <- whole_genome_clump %>%  select(subtypes)

whole_genome_clump %>% mutate(SNP=SNP.ONCO,
                              beta = -score
)




prs <- whole_genome_clump %>% mutate(SNP=SNP.ONCO,
                               beta = -score
                               ) %>%
                     filter(p.min<=1E-06) %>% 
  select(SNP,referece_allele,beta)

