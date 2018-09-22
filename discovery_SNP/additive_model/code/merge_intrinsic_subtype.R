result <- NULL
for(i1 in 1:35){
  load(paste0("./discovery_SNP/additive_model/result/intrinsic_subtype_",i1,".Rdata"))
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


head(CIMBA.BCAC.meta.result)
idx.fil <- which(CIMBA.BCAC.meta.result$MarkerName%in%
                   discovery_snp_new$var_name)
idx.match <- match(discovery_snp_new$var_name,CIMBA.BCAC.meta.result$MarkerName[idx.fil])




CIMBA.BCAC.meta.result.new <- CIMBA.BCAC.meta.result[idx.fil[idx.match],]


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
                            
                              "HER2 Enriched",
                              "Triple Negative"))




write.csv(final.result,file= "./discovery_SNP/additive_model/result/intrinsic_subtype_final.csv")
