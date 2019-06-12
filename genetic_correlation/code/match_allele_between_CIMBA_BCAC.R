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
all.snp <- fread("./ResultsMeta_GWAS_iCOGs_Onco_filter_R2_MAF.txt",header=T)

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


#load CIMBA data
load("./CIMBA/CIMBA.result.Rdata")
head(CIMBA.result)
CIMBA.result = CIMBA.result %>% mutate(chr.pos = paste0(CHR,"_",Position))

#library(data.table)
#merge CIMBA and BCAC data together
snp.infor.merge.update <- merge(snp.infor.merge,CIMBA.result,by.x = "chr.pos")


BCAC.meta.result.new <- list()
BCAC.meta.result.new[[1]] <- snp.infor.merge.update[,2:6]
colnames(BCAC.meta.result.new[[1]]) <- c("SNP",
                                         "CHR",
                                         "Position",
                                         "Reference_allele",
                                         "Effect_allele")
BCAC.meta.result.new[[2]] <- snp.infor.merge.update[,c(12:16,96)]
colnames(BCAC.meta.result.new[[2]]) <- c(colnames(BCAC.meta.result[[2]]),"CIMBA_BRCA1")
head(BCAC.meta.result.new[[1]])
head(BCAC.meta.result.new[[2]])
BCAC.meta.result.new[[3]] <- snp.infor.merge.update[,c(12:36,97)]
colnames(BCAC.meta.result.new[[3]])[26] <- "CIMBA_BRCA1_var"
BCAC.meta.result.new[[4]] <- snp.infor.merge.update[,c(37:38)]
save(BCAC.meta.result.new,file="./two_stage_model_results/BCAC_CIMBABRCA1.Rdata")
# idx <- which(snp.infor.merge.update$alleles3!=snp.infor.merge.update$Allele1)
# snp.infor.merge.update[idx,]

#temp = merge(snp.infor,CIMBA.result,by.x="chr.pos")

#idx <- which(temp$alleles3!=temp$Allele1)
