setwd('/spin1/users/zhangh24/breast_cancer_data_analysis/')
result <- NULL
for(i1 in 1:35){
  load(paste0("./discovery_SNP/additive_model/result/intrinsic_subtype_hr_",i1,".Rdata"))
  result <- rbind(result,test.result.second.wald)
}
#SNP <- c(colnames(icog.julie),colnames(discovery.snp.icog)[1:18])

##################discovery snp were ordered based on the order they are extracted
discovery_snp <- read.csv("/spin1/users/zhangh24/breast_cancer_data_analysis/data/discovery_snp_summary_new.csv",header=T)
#SNP <- discovery_snp$SNP.ICOGS

################match the discovery snps to the order in the paper
discovery_snp_paper_order <- read.csv("./data/discovery_snp_paper_order.csv",header=T)
chr.pos.paper <- paste0(discovery_snp_paper_order$CHR,":",discovery_snp_paper_order$position)


chr.pos <- paste0(discovery_snp$CHR.x,":",discovery_snp$position)

idx.match <- match(chr.pos.paper,
                   chr.pos)
discovery_snp_new <- discovery_snp[idx.match,]
result <- result[idx.match,]
cbind(discovery_snp_new,result)
# colnames(result)[c(2*(1:n.subtypes))-1] <- paste0("logodds_",c("Luminial_A","Luminal_B",
#                                                                    "Luminal_B_HER2-",
#                                                                    "HER2_Enriched",
#                                                                    "Triple_Negative"))
# colnames(result)[c(2*(1:n.subtypes))] <- paste0("p value",c("Luminial A","Luminal B",
#                                                              "Luminal B HER2-",
#                                                              "HER2 Enriched",
#                                                              "Triple Negative"))
#   
#   









n <- nrow(discovery_snp_new)
major.minor <- rep("c",n)
allele.temp <- strsplit(as.character(discovery_snp_new$var_name),"_")
for(i in 1:n){
  if(discovery_snp_new$exp_freq_a1[i]>0.5){
    major.minor[i] <- paste0(allele.temp[[i]][4],"/",allele.temp[[i]][3])  
  }else{
    major.minor[i] <- paste0(allele.temp[[i]][3],"/",allele.temp[[i]][4])  
  }
}
##########order the alleles in terms of major versus minor
discovery_snp_new$marjor.minor <- major.minor




load(paste0("./discovery_SNP/CIMBA_BCAC_meta_analysis/result/CIMBA.BCAC.meta.result.Rdata"))
idx <- which(CIMBA.BCAC.meta.result$CHR==17&
               CIMBA.BCAC.meta.result$position==7571752)
idx <- which(CIMBA.BCAC.meta.result$CHR==17&
               CIMBA.BCAC.meta.result$position==43681771)
CIMBA.BCAC.meta.result[idx,]
head(CIMBA.BCAC.meta.result)
idx.fil <- which(CIMBA.BCAC.meta.result$MarkerName%in%
                   discovery_snp_new$var_name)
idx.match <- match(discovery_snp_new$var_name,CIMBA.BCAC.meta.result$MarkerName[idx.fil])

-0.3495
0.1033
-0.3790843
0.00445735

var_new = (1/0.1033^2 + 1/0.00445735)^-1
beta_new = var_new*(1/0.1033^2*-0.3495+1/0.00445735*-0.3790843)
CI95withP(beta_new,sqrt(var_new))
CIMBA.BCAC.meta.result.new <- CIMBA.BCAC.meta.result[idx.fil[idx.match],]
-0.042938989
0.006703173
CI95withP(-0.042938989,0.006703173)
CI95withP <- function(effect,std){
  z <- effect/std
  p <- 2*pnorm(abs(z),lower.tail = F)
  effect.low <- effect-1.96*std
  effect.high <- effect+1.96*std
  effect.result <- paste0(exp(effect),"(",
                          exp(effect.low),"-",
                          exp(effect.high),")")
  return(list(p,effect.result))
}
#CI95withP(0.06653665,0.011429894)
############Organize CIMBA ORs with per copy of minor alleles
CIMBA.p <- rep(0,n)
CIMBA.OR <- rep("c",n)
for(i in 1:n){
  if(discovery_snp_new$exp_freq_a1[i]>0.5){
    effect <-  -CIMBA.BCAC.meta.result.new$Effect[i]
    std <- CIMBA.BCAC.meta.result.new$StdErr[i]
  }else{
    effect <-  CIMBA.BCAC.meta.result.new$Effect[i]
    std <- CIMBA.BCAC.meta.result.new$StdErr[i]
  }
  result.temp <- CI95withP(effect,std)
  CIMBA.p[i] <- result.temp[[1]]
  CIMBA.OR[i] <- result.temp[[2]]
}

final.result <- cbind(discovery_snp_new$rs_id,
                      result,
                      CIMBA.OR,
                      CIMBA.p)


colnames(final.result) <- c("rs_id",
                            "Luminial A OR",
                            "Luminial A P",
                            "Luminal B OR",
                            "Luminal B P",
                            "Luminal B HER2Neg OR",
                            "Luminal B HER2Neg P",
                            "HER2 Enriched OR",
                            "HER2 Enriched P",
                            "Triple Negative OR",
                            "Triple Negative P",
                            "global association test",
                            "global heterogeneity test",
                            "CIMBA OR",
                            "CIMBA P")




write.csv(final.result,file= "./discovery_SNP/additive_model/result/intrinsic_subtype_final.csv")
