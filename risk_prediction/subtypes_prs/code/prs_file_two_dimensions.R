#-------------------------------------------------------------------
# Update Date: 11/25/2018
# Create Date: 11/22/2018
# Goal: prepare different different prs files for plink
# Author: Haoyu Zhang
#-------------------------------------------------------------------
#---------------------------------------#---------------------------------------
setwd('/data/zhangh24/breast_cancer_data_analysis/risk_prediction')

#load("./EB_whole_genome/result/whole_gonome.rdata")
load("./intrinsic_subtypes_whole_genome/ICOG/result/whole_genome_threeadd.rdata")
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
load("./data/Nasim_330SNPs_complete_information.Rdata")
#check whether the 330 SNPs are in the list
head(snp.new)
idx <- which(snp.new$var_name%in%
               whole_genome_clump$var_name==F)
length(idx)
#put the 330 SNPs stan_p as 0 for later filtering
head(snp.new)
idx <- which(whole_genome_clump$var_name%in%snp.new$var_name==T)
length(idx)
whole_genome_clump$stan_p[idx] = 0

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
# pthres <- c(5E-08,1E-07,5E-07,1E-06,5E-06,1E-05,5E-05,
#             1E-04,5E-04,1E-03,5E-03,1E-02)
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
whole_genome_clump_new <- whole_genome_clump %>% mutate(SNP=SNP.ONCO) %>% select(SNP,effect_allele,stan_p,FTOP_result,p.min,CHR,position) %>% 
  cbind(score)

for(i in 1:n.pthres){
  for(k in 1:n.pthres){
    
    
    #keep the 330 SNPs, then for all other SNPs, make sure only the top one within 500kb are kept
    prs <-  whole_genome_clump_new %>%
      filter(stan_p<=pthres[i]|
               FTOP_result<=pthres[k]) 
    idx <- which(prs$SNP%in%snp.new$SNP.ONCO==F)
    if(length(idx)!=0){
      prs.no.nasim <- prs[idx,,drop=F]
      #remove SNP  rs76858104 which is in LD with Nasim SNP rs11249433
      #remove SNP that are not indepedent with prevoius known SNPs
      qdx <- which(prs.no.nasim$p.min<=8.2E-12)
      prs.no.nasim = prs.no.nasim[-qdx,,drop=F]
      if(length(qdx)<length(idx)){
        snp.keep <- NULL
        for(l in 1:nrow(prs.no.nasim)){
          #check whether there are any SNPs that are within 500kb of the top SNPs
          jdx <- which(prs.no.nasim$CHR==prs.no.nasim$CHR[l]&
                         prs.no.nasim$position>=prs.no.nasim$position[l]-500000&
                         prs.no.nasim$position<=prs.no.nasim$position[l]+500000&           prs.no.nasim$SNP!=prs.no.nasim$SNP[l])
          if(sum(prs.no.nasim$p.min[l]<prs.no.nasim$p.min[jdx])==length(jdx)&
             #remove SNP  rs76858104 which is in LD with Nasim SNP rs11249433
             prs.no.nasim$p.min[l]>=3E-12){
            snp.keep = c(snp.keep,prs.no.nasim$SNP[l])
          }
          
          
        }
        #only keep the 313 SNPs and snps that have no nearby snps in +-500kb
        kdx <- which(prs$SNP%in%snp.new$SNP.ONCO|
                       prs$SNP%in%snp.keep)
        prs = prs[kdx,,drop=F]
      }  
      }
      
    
    
    
    
    for(j in 1:length(select.names)){
      prs.new <-  prs %>%
        select(SNP,effect_allele,select.names[j])
      # prs <-  whole_genome_clump_new %>%
      # filter(stan_p<=pthres[i]|
      #          FTOP_result<=pthres[k]) %>% 
      # select(SNP,effect_allele,select.names[j])
    colnames(prs.new) <- c("SNP","effect_allele","beta")
    dim(prs.new)
    write.table(prs.new,file = paste0("/data/zhangh24/breast_cancer_data_analysis/risk_prediction/subtypes_prs/result/",select.names[j],"_prs_pvaluecut_",i,"_",k,"_121019.file"),row.names=F,col.names=T,quote=F)
    
  }
  }
}






