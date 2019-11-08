#purpose debug
#figure out where is wrong
#why the prs continues increase as more snps get in



setwd('/data/zhangh24/breast_cancer_data_analysis/risk_prediction')

#load("./EB_whole_genome/result/whole_gonome.rdata")
load("./intrinsic_subtypes_whole_genome/ICOG/result/whole_gonome.rdata")

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

#create the prs file for two-stage and eb
subtypes <- c("Luminal_A",
              "Luminal_B",
              "Luminal_B_HER2Neg",
              "HER2_Enriched",
              "TN")
select.names <- subtypes
score <- whole_genome_clump %>%  select(select.names)
#reverse the odds ratio, since dosage use the first allele as reference
score <- -score
whole_genome_clump_new <- whole_genome_clump %>% mutate(SNP=SNP.ONCO) %>% select(SNP,reference_allele,p.min) %>% 
  cbind(score)


#generate extract snps list
i <- j <- 1
    prs <-  whole_genome_clump_new %>%
      filter(p.min<=pthres[i]) %>% 
      select(SNP)
    colnames(prs) <- c("SNP")

  write.table(prs,file = paste0("/data/zhangh24/breast_cancer_data_analysis/risk_prediction/subtypes_prs/result/extrac_snp_list"),row.names=F,col.names=F,quote=F)
  
  
  
 #generate map files 
  load("/data/zhangh24/breast_cancer_data_analysis/whole_genome/ONCO/ERPRHER2GRADE_fixed_baseline/result/onco_info.Rdata")
  #map file have four columns, no header
  #chromosome (1-22, X, Y or 0 if unplaced)
  #rs# or snp identifier
  #Genetic distance (morgans)
  #Base-pair position (bp units)
  onco_info$genetic_distance <- 0
  map_file <- onco_info %>% 
    select(CHR,rs_id,genetic_distance,position)
  write.table(map_file,
              file = "/data/zhangh24/breast_cancer_data_analysis/risk_prediction/subtypes_prs/result/map.txt",
              row.names = F,
              col.names=F,
              quote = F)
  
  extract_rs_id <- rep("c",nrow(prs))
  for(i in 1:nrow(prs)){
    extract_rs_id[i] <- as.character(prs[i,1])
  }
  temp <- which(extract_rs_id%in%
                  as.character(map_file$rs_id))
  
  
  
  prs.code <- rep("c",length(select.names))
  prs.code <- data.frame(prs.code,stringsAsFactors=F)
  temp <- 1
  
      prs.code <- paste0("/data/zhangh24/plink --dosage /data/zhangh24/BCAC/impute_onco_dosage/dosage_all noheader skip0=1 skip1=1 format=1 --map  /data/zhangh24/breast_cancer_data_analysis/risk_prediction/subtypes_prs/result/map.txt --fam /data/zhangh24/BCAC/impute_onco/onco_plink.fam --extract /data/zhangh24/breast_cancer_data_analysis/risk_prediction/subtypes_prs/result/extrac_snp_list --out /data/zhangh24/breast_cancer_data_analysis/risk_prediction/subtypes_prs/result/extract_snps_genotype")
    
      
      "/data/zhangh24/qctool_v1.4-linux-x86_64/qctool -g chr22_test.gen -incl-rsids /data/zhangh24/breast_cancer_data_analysis/risk_prediction/subtypes_prs/result/extrac_snp_list.txt -og /data/zhangh24/BCAC/impute_onco/extracted_snp_test"
      
      
      "/data/zhangh24/software/qctool/qctool -g /data/zhangh24/BCAC/impute_onco/onco_all.gen.gz -incl-rsids /data/zhangh24/breast_cancer_data_analysis/risk_prediction/subtypes_prs/result/extrac_snp_list.txt -og /data/zhangh24/BCAC/impute_onco/extracted_snp_test"   
      /spin1/users/zhangh24/software/qctool/qctool
      
      
 