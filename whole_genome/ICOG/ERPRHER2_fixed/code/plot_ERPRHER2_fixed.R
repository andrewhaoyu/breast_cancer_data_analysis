args <- commandArgs(trailingOnly=T)
args <- as.numeric(args[[1]])
i1 <- args
setwd("/spin1/users/zhangh24/breast_cancer_data_analysis/whole_genome/ICOG/ERPRHER2_fixed/")
library(qqman)
if(i1==1){
  load("/spin1/users/zhangh24/breast_cancer_data_analysis/whole_genome/ICOG/ERPRHER2_fixed/result/meta_result_shared_1p.Rdata")
  pvalue <- meta_result_shared_1p[,15]
  idx <- which(pvalue<=1.0e-20)
  pvalue[idx] <- 1.0e-20
  gwas_result <- data.frame(SNP=as.character(meta_result_shared_1p$rs_id),CHR=meta_result_shared_1p$CHR,BP=meta_result_shared_1p$position,
                            P=pvalue,stringsAsFactors =F)
  png(paste0("./result/plot/man.png"),width = 7.635,height =4.7175,units = "in",res = 600)
  manhattan(gwas_result,suggestiveline = F, cex.axis = 2,main=paste0(" Global Test Association Manhattan Plot"),ylim=c(0,20))
  dev.off()
  
}else if(i1==2){
  load("/spin1/users/zhangh24/breast_cancer_data_analysis/whole_genome/ICOG/ERPRHER2_fixed/result/meta_result_shared_1p_filter_1M.Rdata")
  pvalue <- meta_result_shared_1p_filter[,15]
  idx <- which(pvalue<=1.0e-20)
  pvalue[idx] <- 1.0e-20
  gwas_result <- data.frame(SNP=as.character(meta_result_shared_1p_filter$rs_id),CHR=meta_result_shared_1p_filter$CHR,BP=meta_result_shared_1p_filter$position,
                            P=pvalue,stringsAsFactors =F)
  png(paste0("./result/plot/man_filter.png"),width = 7.635,height =4.7175,units = "in",res = 600)
  manhattan(gwas_result,suggestiveline = F, cex.axis = 2,main=paste0(" Global Test Association Manhattan Plot"),ylim=c(0,20))
  dev.off()
}else{
  load("/spin1/users/zhangh24/breast_cancer_data_analysis/whole_genome/ICOG/ERPRHER2_fixed/result/meta_result_shared_1p_filter_1M_Ju.Rdata")
  meta_result_shared_1p_filter <- meta_result_shared_1p_filter_Ju
  
  pvalue <- meta_result_shared_1p_filter[,15]
  idx <- which(pvalue<=1.0e-20)
  pvalue[idx] <- 1.0e-20
  gwas_result <- data.frame(SNP=as.character(meta_result_shared_1p_filter$rs_id),CHR=meta_result_shared_1p_filter$CHR,BP=meta_result_shared_1p_filter$position,
                            P=pvalue,stringsAsFactors =F)
  

  
  png(paste0("./result/plot/man_filter_Montse.png"),width = 7.635,height =4.7175,units = "in",res = 600)
  manhattan(gwas_result,suggestiveline = F, cex.axis = 2,main=paste0(" Global Test Association Manhattan Plot"),ylim=c(0,20))
  dev.off()
  
}

position.cut = 0.5*10^6


idx.cut.more <- which(gwas_result$P<=5E-08&gwas_result$CHR==5&(gwas_result[,3]>=(46146805-position.cut)&gwas_result[,3]<=(46146805+position.cut)))
idx <- which(gwas_result$P<=5E-08&gwas_result$CHR==1)
idx.cut.more <- c(idx,idx.cut.more)
gwas_result <- gwas_result[-idx.cut.more,]
idx <- which(gwas_result$CHR==15&gwas_result$BP==67408298)
gwas_result[idx,4] <- 5E-08


png(paste0("./result/plot/man_filter_Montse_dirty.png"),width = 7.635,height =4.7175,units = "in",res = 600)
manhattan(gwas_result,suggestiveline = F, cex.axis = 2,main=paste0(" Global Test Association Manhattan Plot"),ylim=c(0,20))
dev.off()












