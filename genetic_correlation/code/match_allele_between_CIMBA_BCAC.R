#Goal: match the allele between CIMBA and BCAC

setwd("/dcl01/chatterj/data/hzhang1/breast_intrinsic/")

#load BCAC data
load("./two_stage_model_results/BCAC.meta.result.Rdata")
#merge BCAC results as snp.infor
snp.infor <- cbind(BCAC.meta.result[[1]],
                   BCAC.meta.result[[2]],
                   BCAC.meta.result[[3]],
                   BCAC.meta.result[[4]])
snp.infor$chr.pos <- paste0(snp.infor$CHR.x.x,"_",snp.infor$BP.x)
library(dplyr)
library(data.table)
#load in all SNPs information
all.snp <- fread("../breast_cancer_data_analysis/discovery_SNP/result/ResultsMeta_GWAS_iCOGs_Onco_filter_R2_MAF.txt",header=T)

all.snp = all.snp %>% 
  mutate(chr.pos=paste0(chr.Onco,"_",Position.Onco))
#merge the two datasets together
snp.infor.merge <- merge(snp.infor,all.snp,
                         by="chr.pos")
#fill in the SNPs with missing alleles
class(snp.infor.merge$alleles3)
snp.infor.merge$alleles3 <- as.character(snp.infor.merge$alleles3)
idx <- which(is.na(snp.infor.merge$alleles3))
snp.infor.merge$alleles3[idx] <- snp.infor.merge$Baseline.Meta[idx]
snp.infor.merge$alleles4 <- as.character(snp.infor.merge$alleles4)
idx <- which(is.na(snp.infor.merge$alleles4))
snp.infor.merge$alleles4[idx] <- snp.infor.merge$Effect.Meta[idx]


#load updated BCAC data whole genome 082119
load("/dcl01/chatterj/data/hzhang1/breast_intrinsic/whole_genome_breast_cancer_results/meta_result_shared_1p_082119.Rdata")
#merge the two datasets
colnames(snp.infor.merge)[12:36] <- paste0("cov",c(1:25))
snp.infor.merge.update <- merge(snp.infor.merge,meta_result_shared_1p,by="var_name")

#reorder the log odds ratio and covariance matrix based on snp.infor.merge
#meta_result_shared_1p used the order c("Luminial_A","Luminal_B","Luminal_B_HER2Neg","HER2_Enriched","Triple_Negative")
#snp.infor.merge used the order c("Triple_Negative","Luminial_A","HER2_Enriched","Luminal_B","Luminal_B_HER2Neg")
log.odds <- snp.infor.merge.update[,105:109][,c(5,1,4,2,3)]
log.odds.sigma <- matrix(0,
                         nrow(log.odds),
                         ncol(log.odds)^2)
for(l in 1:nrow(log.odds)){
  if(l%%1000==0){
    print(l)  
  }
  
  temp.mat <- matrix(as.numeric(snp.infor.merge.update[l,c(110:134)]),ncol(log.odds),ncol(log.odds))
  temp.mat <- temp.mat[c(5,1,4,2,3),c(5,1,4,2,3)]
  log.odds.sigma[l,] <- as.vector(temp.mat)
}
snp.infor.merge.update[,c(8:12)] <- log.odds
snp.infor.merge.update[,c(13:37)] <- log.odds.sigma
snp.infor.merge.temp <- snp.infor.merge.update[,c(1:90)]
# names.ori <- colnames(snp.infor.merge)
# names.temp <- colnames(snp.infor.merge.temp)
# idx.match <- match(names.ori,names.temp)


 %>% %>% 
head(CIMBA.result)
str.temp <- strsplit(CIMBA.result$MarkerName,"_")
n <- nrow(CIMBA.result)
CHR <- rep(0,n)
position <- rep(0,n)
for(i in 1:n){
  CHR[i] <- as.numeric(str.temp[[i]][1])
  position[i] <- as.numeric(str.temp[[i]][2])
}
CIMBA.result$CHR <- CHR
CIMBA.result$Position <- position
#CIMBA.result = CIMBA.result %>% mutate(chr.pos = paste0(CHR,"_",Position))

#library(data.table)
#merge CIMBA and BCAC data together
snp.infor.merge.update <- merge(snp.infor.merge.temp,CIMBA.result,by.x = "var_name",by.y = "MarkerName")

idx <- which(as.character(snp.infor.merge.update$alleles4)!=snp.infor.merge.update$EffectAllele)

snp.infor.merge.update <- snp.infor.merge.update[-idx,]

BCAC.meta.result.new <- list()
BCAC.meta.result.new[[1]] <- snp.infor.merge.update[,3:7]
colnames(BCAC.meta.result.new[[1]]) <- c("SNP",
                                         "CHR",
                                         "Position",
                                         "Reference_allele",
                                         "Effect_allele")
BCAC.meta.result.new[[2]] <- snp.infor.merge.update[,c(8:12,94)]
colnames(BCAC.meta.result.new[[2]]) <- c(colnames(BCAC.meta.result[[2]]),"CIMBA_BRCA1")
head(BCAC.meta.result.new[[1]])
head(BCAC.meta.result.new[[2]])
BCAC.meta.result.new[[3]] <- snp.infor.merge.update[,c(13:37,95)]
colnames(BCAC.meta.result.new[[3]])[26] <- "CIMBA_BRCA1_var"
BCAC.meta.result.new[[4]] <- snp.infor.merge.update[,c(38:39,93)]
colnames(BCAC.meta.result.new[[4]])[3] <- "freq_CIMBA"
save(BCAC.meta.result.new,file="./two_stage_model_results/BCAC_CIMBABRCA1_082119.Rdata")
# idx <- which(snp.infor.merge.update$alleles3!=snp.infor.merge.update$Allele1)
# snp.infor.merge.update[idx,]

#temp = merge(snp.infor,CIMBA.result,by.x="chr.pos")

#idx <- which(temp$alleles3!=temp$Allele1)
