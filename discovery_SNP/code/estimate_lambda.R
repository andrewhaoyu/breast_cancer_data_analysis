#Goal: estimate lambda for qq-plot
library(data.table)
setwd("/dcl01/chatterj/data/hzhang1/breast_cancer_data_analysis")

Calculatelambda <- function(x,stat_type){
  if(stat_type!="CHISQ4"&
     stat_type!="CHISQ5"){
    if (stat_type == "Z")
      z = x
    
    if (stat_type == "CHISQ")
      z = sqrt(x)
    
    if (stat_type == "PVAL")
      z = qnorm(x / 2)
    lambda = round(median(z^2) / qchisq(0.5,1), 3)
    return(lambda)
    
  }else{
    if(stat_type =="CHISQ5"){ z = qchisq(x,5,lower.tail = F)
    lambda = round(median(z) / qchisq(0.5,5), 3)}
    if(stat_type =="CHISQ4"){  z = qchisq(x,4,lower.tail = F)
    lambda = round(median(z) / qchisq(0.5,4), 3)}
    return(lambda)
  }
  
}
FilterSNP <- function(gwas_result,fine_mapping){
  idx_cut <- NULL
  start <- fine_mapping$start
  end <- fine_mapping$end
  CHR <- fine_mapping$CHR
  
  #fine_mapping
  for(i in 1:nrow(fine_mapping)){
    print(i)
    chr_temp <- CHR[i]
    start_temp <- start[i]
    end_temp <- end[i]
    idx <- which(gwas_result$CHR==chr_temp&gwas_result$BP>=start_temp&gwas_result$BP<=end_temp)
    idx_cut <- c(idx_cut,idx)
  }
  ############duplicate variables remove by unqiue function
  idx_cut <- unique(idx_cut)
  gwas_result_filter <- gwas_result[-idx_cut,]
  return(gwas_result_filter)
}












data <- as.data.frame(fread("./discovery_SNP/result/ResultsMeta_GWAS_iCOGs_Onco_filter_R2_MAF.txt"))
library(dplyr)
gwas_result <- data %>%
  #filter(p.meta<=10^-4) %>% 
  select(c(SNP.Onco,chr.Onco,Position.Onco,p.meta))
colnames(gwas_result) <- c("SNP",
                           "CHR",
                           "BP",
                           "P")

#standard analysis
Calculatelambda(gwas_result$P,"PVAL")
#standard analysis after filtering
fine_mapping <- read.csv("./data/filter_regions_standard.csv",header= T)
gwas_result_filter <- FilterSNP(gwas_result,fine_mapping)
Calculatelambda(gwas_result_filter$P,"PVAL")




#fixed-effet two-stage
meta_result_shared_1p <- as.data.frame(fread("/dcl01/chatterj/data/hzhang1/breast_intrinsic/meta_result_shared_1p_fixed.txt",header=T))
ftop <- meta_result_shared_1p
ftop.p <- ftop$p.value
subtypes.p <- ftop.p
ftop$subtypes.p <- ftop$p.value
ftop =  ftop %>% 
  mutate(subtypes.p.new = ifelse((is.nan(subtypes.p)==T)|(subtypes.p==0),1E-20,subtypes.p))
#remove SNPs with missing id
idx <- which(is.na(subtypes.p))
ftop <- ftop[-idx,]
subtypes_gwas_result <- ftop %>%
  select(rs_id,CHR,position,subtypes.p.new)
colnames(subtypes_gwas_result) <- c("SNP",
                                    "CHR",
                                    "BP",
                                    "P")
#fixed-effect two-stage polytomous model
Calculatelambda(subtypes_gwas_result$P,"CHISQ5")
fine_mapping <- read.csv("./data/filter_regions_subtypes.csv",header= T)
subtypes_gwas_result_filter <- FilterSNP(subtypes_gwas_result,fine_mapping) 
#fixed-effect two-stage polytomous model after filtering
Calculatelambda(subtypes_gwas_result_filter$P,"CHISQ5")











#mixed-effet two-stage
meta_result_shared_1p <- as.data.frame(fread("/dcl01/chatterj/data/hzhang1/breast_intrinsic/meta_result_shared_1p_mixed.txt",header=T))
mtop <- meta_result_shared_1p
mtop.p <- mtop$p.value
#just for easy coding, put mtop as ftop
ftop = mtop
ftop$subtypes.p <- ftop$p.value
ftop =  ftop %>% 
  mutate(subtypes.p.new = ifelse((is.nan(subtypes.p)==T)|(subtypes.p==0),1E-20,subtypes.p))
#remove SNPs with missing id
idx <- which(is.na(subtypes.p))
ftop <- ftop[-idx,]
subtypes_gwas_result <- ftop %>%
  select(rs_id,CHR,position,subtypes.p.new)
colnames(subtypes_gwas_result) <- c("SNP",
                                    "CHR",
                                    "BP",
                                    "P")
#mixed-effect two-stage polytomous model analysis
Calculatelambda(subtypes_gwas_result$P,"CHISQ4")
#standard analysis after filtering
fine_mapping <- read.csv("./data/filter_regions_subtypes.csv",header= T)
subtypes_gwas_result_filter <- FilterSNP(subtypes_gwas_result,fine_mapping) 

Calculatelambda(subtypes_gwas_result_filter$P,"CHISQ4")





#CIMBA and BCAC meta-analysis
cimba_result_all <- as.data.frame(fread("/dcl01/chatterj/data/hzhang1/breast_intrinsic/brca1_bcac_tn_meta.txt",header = T))

colnames(cimba_result_all)[10] <- "P"
cimba_result = cimba_result_all %>% 
  filter(Freq1>=0.008&
           Freq1<=0.992&
           CHR!=23) %>% 
  select(MarkerName,CHR,
         position,P)

colnames(cimba_result) <- c("SNP",
                            "CHR",
                            "BP",
                            "P")
idx <- which(is.na(cimba_result$P))
length(idx)
#CIMBA and BCAC TN analysis
#cimba_result$P
Calculatelambda(cimba_result$P,"PVAL")
#CIMBA and BCAC TN analysis after filtering


fine_mapping <- read.csv("./data/filter_regions_cimba.csv",header= T)
cimba_result_filter <- FilterSNP(cimba_result,fine_mapping) 
Calculatelambda(cimba_result_filter$P,"PVAL")
