setwd('/data/zhangh24/breast_cancer_data_analysis/')
library(data.table)
###########load CIMBA data
CIMBA <- fread('./data/brca1_bc.txt',header=T)
###########load BCAC intrinsic subtype data
load(paste0("/data/zhangh24/breast_cancer_data_analysis/whole_genome_age/ICOG/Intrinsic_subtypes/result/meta_result_shared_1p.Rdata"))
###############Transform the allele in different coding
TranformAllele <- function(allele){
  if(allele=="A"){
    return("T")
  }else if(allele=="T"){
    return("A")
  }else if(allele=="C"){
    return("G")
  }else if(allele=="G"){
    return("C")
  }else{
    return("-")
  }
}

################change CIMBA MarkerName into BCAC rs_id style
#CIMBA_rs_id <- gsub("_",":",CIMBA$MarkerName)
#head(CIMBA_rs_id)
#CIMBA$rs_id <- CIMBA_rs_id
head(CIMBA)
colnames(meta_result_shared_1p)[21:45] <- paste0("var",c(1:25))
############clean CIMBA data to the right order
Allele.CIMBA.effect <- toupper(CIMBA$Allele1)
Allele.CIMBA.ref <- toupper(CIMBA$Allele2)
snp.split <- strsplit(CIMBA$MarkerName,"_")
n <- nrow(CIMBA)
Allele1.CIMBA <- rep("c",n)
Allele2.CIMBA <- rep("c",n)
effect <- CIMBA$Effect
notice <- NULL
Freq.update <- CIMBA$Freq1
for(i in 1:n){
  if(i%%10000==0){
    print(i)
  }
  Allele1.CIMBA[i] <- snp.split[[i]][3]
  Allele2.CIMBA[i] <- snp.split[[i]][4]
  if(snp.split[[i]][3]==Allele.CIMBA.effect[i]){
    effect[i] = -effect[i]
    Freq.update[i] = 1-CIMBA$Freq1[i]
  }else if(snp.split[[i]][4]==Allele.CIMBA.effect[i]){
    
  }else{
    if(snp.split[[i]][3]==TranformAllele(Allele.CIMBA.effect[i])){
      effect[i] = -effect[i]
      Freq.update[i] = 1-CIMBA$Freq1[i]
    }else if(snp.split[[i]][4]==TranformAllele(Allele.CIMBA.effect[i])){
      
    }else{
      notice <- c(notice,i)  
    }
  }
}
#CIMBA$Allele2[notice]

CIMBA.clean <- data.frame(CIMBA$MarkerName,
                          Allele1.CIMBA,
                          Allele2.CIMBA,
                          Freq.update,
                          effect,
                          CIMBA$StdErr,
                          stringsAsFactors = F
                          )
colnames(CIMBA.clean) <- c("MarkerName",
                           "Reference_Allele",
                           "Effect_Allele",
                           "Freq_Effect_Allele",
                           "Effect",
                           "StdErr")

save(CIMBA.clean, file = paste0("./discovery_SNP/CIMBA_BCAC_meta_analysis/result/CIMBA.clean.Rdata"))
##############combine CIMBA and BCAC
combine.data <- merge(CIMBA.clean,meta_result_shared_1p,by.x = "MarkerName",
                      by.y = "var_name")
freq.diff = combine.data$Freq_Effect_Allele-combine.data$exp_freq_a1
###########remove the difference freq SNPs
idx <- which(freq.diff>=0.3)

CIMBA.BCAC.combine <- combine.data[-idx,]
CIMBA.BCAC.combine <- as.data.frame(CIMBA.BCAC.combine)
##############save BCAC and CIMBA combine data for meta-analysis
save(CIMBA.BCAC.combine,file = paste0("./discovery_SNP/CIMBA_BCAC_meta_analysis/result/CIMBA.BCAC.combine.Rdata"))
