##########merge the CIMBA and BCAC meta analysis result together
setwd('/spin1/users/zhangh24/breast_cancer_data_analysis/')
load(paste0("./discovery_SNP/CIMBA_BCAC_meta_analysis/result/CIMBA.BCAC.combine.Rdata"))
total <- nrow(CIMBA.BCAC.combine)
meta.result <- matrix(0,total,30)
total <- 0
for(i1 in 1:1000){
  print(i1)
  load(paste0("./discovery_SNP/CIMBA_BCAC_meta_analysis/result/meta.result.sub",i1,".Rdata"))
  temp <- nrow(meta.result.sub)
  meta.result[total+(1:temp),] <- meta.result.sub
  total <- total+temp
}

CIMBA.BCAC.combine[,c(21:50)] <- meta.result
CIMBA.BCAC.meta.result <- CIMBA.BCAC.combine

save(CIMBA.BCAC.meta.result,file=paste0("./discovery_SNP/CIMBA_BCAC_meta_analysis/result/CIMBA.BCAC.meta.result.Rdata"))


############organize the data into genetic correlation analysis format
load(paste0("./discovery_SNP/CIMBA_BCAC_meta_analysis/result/CIMBA.BCAC.meta.result.Rdata"))

load(paste0("./genetic_correlation/ICOG/result/meta.result.completeglm.Rdata"))
load(paste0("./genetic_correlation/ICOG/result/ICOG.result.Rdata"))
meta.result <- meta.result.completeglm
chr.pos <- paste0(meta.result[[1]][,2],"_",
                  meta.result[[1]][,3],"_",
                  meta.result[[1]][,4],"_",
                  meta.result[[1]][,5])
##########Take out the CIMBA BCAC SNPs for genetic correlation analysis
snp.infor <- meta.result[[1]]
snp.infor$chr.pos <- chr.pos

CIMBA.BCAC.meta.temp <- merge(snp.infor,CIMBA.BCAC.meta.result,
                              by.x = "chr.pos",
                              by.y = "MarkerName")
snp.infor.meta <- CIMBA.BCAC.meta.temp[,c(2,3,4,5,6)]
colnames(CIMBA.BCAC.meta.temp)
colnames(snp.infor.meta) <- colnames(meta.result[[1]])
colnames(snp.infor.meta)[4:5] <- c("Reference_allele",
                                   "Effect_allele")
log.odds.meta <- CIMBA.BCAC.meta.temp[,c(26:30)]
colnames(log.odds.meta) <- colnames(meta.result[[2]])
head(log.odds.meta)
sigma.meta <- CIMBA.BCAC.meta.temp[,c(31:55)]
CIMBA.BCAC.meta <- list(snp.infor.meta,
                       log.odds.meta,
                       sigma.meta)
save(CIMBA.BCAC.meta,file = 
       paste0("./genetic_correlation/ICOG/result/CIMBA.BCAC.meta.Rdata"))

load(paste0("./genetic_correlation/ICOG/result/meta.result.Rdata"))
BCAC.meta.result <- meta.result
save(BCAC.meta.result,file = paste0("./genetic_correlation/ICOG/result/BCAC.meta.result.Rdata"))
#head(CIMBA.BCAC.meta.result)


#######Take out the CIMBA SNPs for genetic correlation analysis
load(paste0("./discovery_SNP/CIMBA_BCAC_meta_analysis/result/CIMBA.BCAC.combine.Rdata"))
load(paste0("./genetic_correlation/ICOG/result/meta.result.completeglm.Rdata"))
snp.infor <- meta.result.completeglm[[1]]
chr.pos <- paste0(meta.result.completeglm[[1]][,2],"_",
                  meta.result.completeglm[[1]][,3],"_",
                  meta.result.completeglm[[1]][,4],"_",
                  meta.result.completeglm[[1]][,5])
snp.infor$chr.pos <- chr.pos
CIMBA.temp <- merge(snp.infor,CIMBA.BCAC.combine,
                    by.x = "chr.pos",
                    by.y = "MarkerName")
CIMBA.result <- CIMBA.temp[,c(2,3,4,5,6,10,11)]
CIMBA.result[,7] <- CIMBA.result[,7]^2
colnames(CIMBA.result) <- c("SNP",
                            "CHR",
                            "Position",
                            "Allele1",
                            "Allele2",
                            "Triple_negative_log_odds_ratio",
                            "Triple_negative_var")
save(CIMBA.result,file = paste0("./genetic_correlation/CIMBA/result/CIMBA.result.Rdata"))
head(CIMBA.temp)
