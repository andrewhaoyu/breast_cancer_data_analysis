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

#create the prs file for two-stage and eb
subtypes <- c("Luminal_A",
              "Luminal_B",
              "Luminal_B_HER2Neg",
              "HER2_Enriched",
              "TN")
select.names <- subtypes
score <- whole_genome_clump %>%  select(select.names)
whole_genome_clump_new <- whole_genome_clump %>% mutate(SNP=SNP.ONCO) %>% select(SNP,effect_allele,stan_p,FTOP_result) %>% 
  cbind(score)


EBEst <- function(sigma_p,beta,sigma_d){
  beta_est <- solve(solve(sigma_p)+solve(sigma_d))%*%solve(sigma_d)%*%beta
  return(beta_est)
}




for(i in 1:n.pthres){
  for(k in 1:n.pthres){
    logodds <- whole_genome_clump %>% 
      filter(stan_p<=pthres[i]|
               FTOP_result<=pthres[k]) %>% select(select.names)
    sigma_d_all <- whole_genome_clump %>% 
      filter(stan_p<=pthres[i]|
               FTOP_result<=pthres[k]) %>% select(24:48)
    
    sigma_p <- cov(logodds)
    eb_logodds <- logodds
    for(z in 1:nrow(eb_logodds)){
      beta = as.numeric(logodds[z,])
      sigma_d = matrix(as.numeric(sigma_d_all[z,]),5,5)
      eb_logodds[z,] <- EBEst(sigma_p,beta,sigma_d)
    }
    
    for(j in 1:length(select.names)){
   
  
      
      
      
      
        info <-  whole_genome_clump_new %>%
        filter(stan_p<=pthres[i]|
                 FTOP_result<=pthres[k]) %>% 
        select(SNP,effect_allele)
        prs_score <- eb_logodds[,j]
        prs <- cbind(info,prs_score)
      colnames(prs) <- c("SNP","effect_allele","beta")
      dim(prs)
      write.table(prs,file = paste0("/data/zhangh24/breast_cancer_data_analysis/risk_prediction/subtypes_prs/result/",select.names[j],"_prs_pvaluecut_",i,"_",k,"_eb_121019.file"),row.names=F,col.names=T,quote=F)
      
    }
  }
}






