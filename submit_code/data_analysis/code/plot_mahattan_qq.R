args <- commandArgs(trailingOnly = T)
i1 <- as.numeric(args[[1]])
setwd("/data/zhangh24/breast_cancer_data_analysis/poland/result/whole_genome/")
library(qqman)
load("/data/zhangh24/breast_cancer_data_analysis/poland/result/whole_genome/onco_result_fixed_5p.Rdata")
load("/data/zhangh24/breast_cancer_data_analysis/poland/result/whole_genome/onco_result_casecase_5p.Rdata")
load("/data/zhangh24/breast_cancer_data_analysis/poland/result/whole_genome/onco_result_standard_5p.Rdata")
if(i1 == 1){
  idx <- which(onco_result_standard_5p$p_value1<=5E-08)
  onco_result_standard_5p[idx,]
  idx2 <- which(onco_result_casecase_5p$P<=5E-08)
  onco_result_casecase_5p[idx2,]
  onco_result_fixed_5p
  gwas_result1 <- data.frame(SNP=as.character(onco_result_fixed_5p$rs_id),CHR=onco_result_fixed_5p$CHR,BP=onco_result_fixed_5p$position,
                             P=onco_result_fixed_5p$P,stringsAsFactors =F)
  
  png(paste0("./man_fix.png"),width = 7.635,height =4.7175,units = "in",res = 600)
  manhattan(gwas_result1,suggestiveline = F, cex.axis = 2,main=paste0("Global Test Association Manhattan Plot"),ymax=12)
  dev.off()
  png(paste0("./qq_fix.png"),width = 7.635,height =4.7175,units = "in",res = 600)
  qq(gwas_result1$P,main=paste0("Global Test Association Manhattan QQ Plot"),ylim=c(0,12))
  dev.off()
  
}else if(i1==2){
  gwas_result2 <- data.frame(SNP=as.character(onco_result_casecase_5p$rs_id),CHR=onco_result_casecase_5p$CHR,BP=onco_result_casecase_5p$position,
                             P=onco_result_casecase_5p$P,stringsAsFactors =F)
  
  png(paste0("./man_random.png"),width = 7.635,height =4.7175,units = "in",res = 600)
  manhattan(gwas_result2,suggestiveline = F, cex.axis = 2,main=paste0("Global Test Association Manhattan Plot"),ymax=12)
  dev.off()
  png(paste0("./qq_random.png"),width = 7.635,height =4.7175,units = "in",res = 600)
  qq(gwas_result2$P,main=paste0("Global Test Association Manhattan QQ Plot"),ylim=c(0,12))
  dev.off()
  
}else if(i1==3){
  gwas_result3 <- data.frame(SNP=as.character(onco_result_standard_5p$rs_id),CHR=onco_result_standard_5p$CHR,BP=onco_result_standard_5p$position,
                             P=onco_result_standard_5p$p_value1,stringsAsFactors =F)
  png(paste0("./man_standard.png"),width = 7.635,height =4.7175,units = "in",res = 600)
  manhattan(gwas_result3,suggestiveline = F, cex.axis = 2,main=paste0("Global Test Association Manhattan Plot"), ylim=c(0,12))
  dev.off()
  png(paste0("./qq_standard.png"),width = 7.635,height =4.7175,units = "in",res = 600)
  qq(gwas_result3$P,main=paste0("Global Test Association Manhattan QQ Plot"),ylim=c(0,12))
  dev.off()
  
  
}else if(i1==4){
  # idx <- which(onco_result_standard_5p$exp_freq_a1>=0.2&
  #                onco_result_standard_5p$exp_freq_a1<=0.8)
  # onco_result_standard_20p <- onco_result_standard_5p[idx,]
  
  p_value2 <- onco_result_standard_5p$p_value2
  p_value2[is.nan(p_value2)] <- 1
  # too significant due the the insufficient 
  idx <- which(p_value2 <= 1E-10)
  p_value2[idx] <- 1
  gwas_result4 <- data.frame(SNP=as.character(onco_result_standard_5p$rs_id),CHR=onco_result_standard_5p$CHR,BP=onco_result_standard_5p$position,
                             P=p_value2,stringsAsFactors =F)
  idx <- which(p_value2==1)
  gwas_result4 <- gwas_result4[-idx,]
  png(paste0("./man_poly.png"),width = 7.635,height =4.7175,units = "in",res = 600)
  manhattan(gwas_result4,suggestiveline = F, cex.axis = 2,main=paste0("Global Test Association Manhattan Plot"), ylim=c(0,12))
  dev.off()
  png(paste0("./qq_poly.png"),width = 7.635,height =4.7175,units = "in",res = 600)
  qq(gwas_result4$P,main=paste0("Global Test Association Manhattan QQ Plot"),ylim=c(0,12))
  dev.off()
  
  
}











