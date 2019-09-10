#-------------------------------------------------------------------
# Goal: estimate heritability for breast cancer overall risk using summary level statistics from icogs, onco and gwas
#-------------------------------------------------------------------
#---------------------------------------#---------------------------------------

SigmaEst <- function(data){
  #first column is log odds ratio
  #second column is variance
  #third column is allele frequency
  idx <- which(is.na(data[,1]))
  if(length(idx!=0)){
    data <- data[-idx,]  
  }
  
  beta <- data[,1]
  sigma <- data[,2]
  p <- data[,3]
  
  n <- length(p)
  result <- sum(2*p*(1-p)*(beta^2-sigma))
  return(result)
}

LD_pruning = function(sig_SNPs,LD2){
  sig_SNPs_temp =sig_SNPs
  filter_result = NULL
  LD2.temp =LD2
  temp.ind = 1
  while(nrow(sig_SNPs_temp)!=0){
    
    idx = which.min(sig_SNPs_temp$p.value)
    filter_result = rbind(filter_result,sig_SNPs_temp[idx,])
    LD2.single = LD2.temp[idx,]
    idx.cut = which( LD2.single>=0.1)
    position.range <- 500*10^3
    filter_result_position = sig_SNPs_temp$position[idx]
    filter_CHR = sig_SNPs_temp$CHR[idx]
    idx.cut2 <- which((sig_SNPs_temp$position>=filter_result_position-position.range)&(sig_SNPs_temp$position<=filter_result_position+position.range)&(sig_SNPs_temp$CHR==filter_CHR ))
    idx.cut <- c(idx.cut,idx.cut2)
    idx.cut <- unique(idx.cut)
    LD2.temp = LD2.temp[-idx.cut,-idx.cut]
    LD2.temp = as.matrix(LD2.temp)
    sig_SNPs_temp = sig_SNPs_temp[-idx.cut,]
    temp.ind = temp.ind+1
  }
  return(filter_result)
}



library(data.table)
setwd("/spin1/users/zhangh24/breast_cancer_data_analysis/")

#load summary level statistics for overall
#Haoyupdate added the meta-analysis between iCOGs and OncoArray
standard_result <- as.data.frame(fread("./discovery_SNP/result/ResultsMeta_GWAS_iCOGs_Onco_filter_R2_MAF_Haoyuupdate.txt"))
colnames(standard_result)[1] <- "MarkerName"
#load CIMBA data
##cimba only have markername
#need to match cimba data to bcac data
#the effect alleles of cimba data is not matching with bcac
#but it doesn't affect the heritability estimate

CIMBA <- as.data.frame(fread("./data/brca1_bc.txt",header=T))
colnames(CIMBA)[5] = "log_or_CIMBA"
CIMBA = CIMBA %>% mutate(var_CIMBA=StdErr^2)

standard_result <- left_join(standard_result,
                             CIMBA,
                             by="MarkerName")

#generate the summary level statistics based on iCOGS and ONCOArray
# Twometa <- function(beta1,var1,beta2,var2){
#   var_meta <- 1/(1/var1+1/var2)
#   beta_meta <- (var_meta)*(beta1/var1+
#                                beta2/var2)
#   return(list(beta_meta,var_meta))
# }
# library(dplyr)
# 
# 
# standard_result = standard_result %>% 
#   mutate(BCAC_meta_beta = Twometa(beta.iCOGs,SE.iCOGs^2,beta.Onco,SE.Onco^2)[[1]],
#          BCAC_meta_var = Twometa(beta.iCOGs,SE.iCOGs^2,beta.Onco,SE.Onco^2)[[2]],
#          Z = BCAC_meta_beta/sqrt(BCAC_meta_var),
#          P_BCAC_meta = 2*pnorm(-abs(Z))
#          )
# write.table(standard_result,file = "./discovery_SNP/result/ResultsMeta_GWAS_iCOGs_Onco_filter_R2_MAF_Haoyuupdate.txt",row.names = F,
# col.names=T,quote=F)

#load summary level statistics for subtypes
# load("/dcl01/chatterj/data/hzhang1/breast_intrinsic/whole_genome_breast_cancer_results/BCAC_subtypes_result.Rdata")







##################discovery snp were ordered based on the order they are extracted
discovery_snp <- read.csv("./data/discovery_snp_summary_new.csv",header=T,stringsAsFactors = F)
##create subtypes summary level statistics for 32 new discovery SNPs
log_or <- matrix(0,35,5)
colnames(log_or) <- paste0("log_or_",c("Luminial A","Luminal B",
                                       "Luminal B HER2Neg ",
                                       "HER2 Enriched",
                                       "Triple Negative"))
var_or <- matrix(0,35,5)
colnames(var_or) <- paste0("var_",c("Luminial A","Luminal B",
                                    "Luminal B HER2Neg ",
                                    "HER2 Enriched",
                                    "Triple Negative"))
for(i1 in 1:35){
  load(paste0("./discovery_SNP/additive_model/result/intrinsic_subtype_logodds",i1,".Rdata"))
  log_or[i1,] <- result[[1]]
  var_or[i1,] <- diag(result[[2]])
}
discovery_snp <- cbind(discovery_snp,log_or,var_or)


discovery_snp_paper_order <- read.csv("./data/discovery_snp_paper_order.csv",header=T,stringsAsFactors = F)
chr.pos.paper <- paste0(discovery_snp_paper_order$CHR,":",discovery_snp_paper_order$position)
library(dplyr)

chr.pos <- paste0(discovery_snp$CHR.x,":",discovery_snp$position)

idx.match <- match(chr.pos.paper,
                   chr.pos)
discovery_snp_new <- discovery_snp[idx.match,]
colnames(discovery_snp_new)[17] <- "MarkerName"


discovery_snp_subtypes <- left_join(discovery_snp_new,
                                    standard_result,
                                    by="MarkerName")
##create standard summary level statistics for 32 new discovery SNPs
standard_result = standard_result %>% 
  mutate(chr.pos = paste0(chr.Onco,":",Position.Onco))
discovery_snp_standard <- left_join(discovery_snp_new,standard_result,by="chr.pos")
##create subtypes summary level statistics for 178 known SNPs
log_or <- matrix(0,178,5)
colnames(log_or) <- paste0("log_or_",c("Luminial A","Luminal B",
                                       "Luminal B HER2Neg ",
                                       "HER2 Enriched",
                                       "Triple Negative"))
var_or <- matrix(0,178,5)
colnames(var_or) <- paste0("var_",c("Luminial A","Luminal B",
                                    "Luminal B HER2Neg ",
                                    "HER2 Enriched",
                                    "Triple Negative"))
for(i1 in 1:178){
  load(paste0("./known_SNPs/known_SNPs_analysis_G_revised/intrinsic_subtypes_pc_additive/result/heter_result_origin",i1,".Rdata"))
  log_or[i1,] <- heter.result[[1]]
  var_or[i1,] <- diag(heter.result[[2]])
}





fine_map <- read.csv("./data/fine_mapping_annotated_clean.csv")
fine_map <- fine_map %>% 
  mutate(chr.pos = paste0(CHR,":",position))
fine_map <- cbind(fine_map,log_or,var_or)
known_snp_subtypes <- left_join(fine_map,standard_result,
                                by="chr.pos")


#two snps duplicated
known_snp_subtypes <- known_snp_subtypes[-c(62,120),]
colnames(known_snp_subtypes)[84] <- "log_or_CIMBA"
#create standard summary level statistics for 178 known SNPs
known_snp_standard <- left_join(fine_map,standard_result,
                                by="chr.pos")
known_snp_standard <- known_snp_standard[-c(62,120),]

#LD pruning for known SNPs
data2 <- as.data.frame(fread("./data/Onco_euro_v10_10232017.csv",header=T))
names2 = colnames(data2)
idxi1 = which(names2=="rs554219")
x.test.all.mis2 <- data2[,c(27:203,idxi1)]

#p.meta is based on icogs,onco,gwas meta-analysis
p.value.known <- known_snp_standard$p.meta

idx.control <- which(data2$Behaviour1==0)

LD2 <- cor(x.test.all.mis2[idx.control,])^2
known.snp.infor.p <- known_snp_standard %>% 
  select(Best.published.SNP,CHR,position,p.meta)

colnames(known.snp.infor.p) <- c("SNP","CHR","position",
                                 "p.value")

known.snp.infor.pruned <- LD_pruning(known.snp.infor.p,LD2)

idx.fil <- which(known.snp.infor.p$SNP%in%known.snp.infor.pruned$SNP)
#only keep SNPs after LD pruning
known_snp_subtypes <- known_snp_subtypes[idx.fil,]
known_snp_standard <- known_snp_standard[idx.fil,]

known_snp_standard_cal = known_snp_standard %>% 
  select(Beta.meta,var.meta,EAFcontrols.Onco)

#heritability for known snps
SigmaEst(known_snp_standard_cal)
#heritability for discovery snps
discovery_snp_cal <- discovery_snp_standard %>% 
  select(Beta.meta,var.meta,EAFcontrols.Onco)
SigmaEst(discovery_snp_cal)
SigmaEst(known_snp_standard_cal)+SigmaEst(discovery_snp_cal)

