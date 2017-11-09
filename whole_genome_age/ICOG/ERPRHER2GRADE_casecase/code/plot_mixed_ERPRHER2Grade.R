args <- commandArgs(trailingOnly=T)
args <- as.numeric(args[[1]])
i1 <- args
setwd("/spin1/users/zhangh24/breast_cancer_data_analysis/whole_genome_age/ICOG/ERPRHER2GRADE_casecase/")
library(qqman)
if(i1==1){
  load("/spin1/users/zhangh24/breast_cancer_data_analysis/whole_genome_age/ICOG/ERPRHER2GRADE_casecase/result/meta_result_shared_1p.Rdata")
  pvalue <- meta_result_shared_1p[,15]
  idx <- which(pvalue<=1.0e-20)
  pvalue[idx] <- 1.0e-20
  gwas_result <- data.frame(SNP=as.character(meta_result_shared_1p$rs_id),CHR=meta_result_shared_1p$CHR,BP=meta_result_shared_1p$position,
                            P=pvalue,stringsAsFactors =F)
  png(paste0("./result/plot/man.png"),width = 7.635,height =4.7175,units = "in",res = 600)
  manhattan(gwas_result,suggestiveline = F, cex.axis = 2,main=paste0(" Global Test Association Manhattan Plot"),ylim=c(0,20))
  dev.off()
  
}else if(i1==2){
  load("/spin1/users/zhangh24/breast_cancer_data_analysis/whole_genome_age/ICOG/ERPRHER2GRADE_casecase/result/meta_result_shared_1p_filter_1M.Rdata")
  #load("/spin1/users/zhangh24/breast_cancer_data_analysis/whole_genome_age/ICOG/ERPRHER2GRADE_fixed_baseline/result/extract.list.ld")
  pvalue <- meta_result_shared_1p_filter[,15]
  idx <- which(pvalue<=1.0e-20)
  pvalue[idx] <- 1.0e-20
  gwas_result <- data.frame(SNP=as.character(meta_result_shared_1p_filter$rs_id),CHR=meta_result_shared_1p_filter$CHR,BP=meta_result_shared_1p_filter$position,
                            P=pvalue,stringsAsFactors =F)
  #idx.cut <- which((gwas_result$BP%in%extract.list.ld$position)&(gwas_result$CHR%in%extract.list.ld$CHR))
  #gwas_result <- gwas_result[-idx.cut,]
  png(paste0("./result/plot/man_filter.png"),width = 7.635,height =4.7175,units = "in",res = 600)
  manhattan(gwas_result,suggestiveline = F, cex.axis = 2,main=paste0(" Global Test Association Manhattan Plot"),ylim=c(0,20))
  dev.off()
}else{
  load("/spin1/users/zhangh24/breast_cancer_data_analysis/whole_genome_age/ICOG/ERPRHER2GRADE_casecase/result/meta_result_shared_1p_filter_1M_Ju.Rdata")
 # load("/spin1/users/zhangh24/breast_cancer_data_analysis/whole_genome_age/ICOG/ERPRHER2GRADE_fixed_baseline/result/extract.list.ld")
  meta_result_shared_1p_filter <- meta_result_shared_1p_filter_Ju
  
  pvalue <- meta_result_shared_1p_filter[,15]
  idx <- which(pvalue<=1.0e-20)
  pvalue[idx] <- 1.0e-20
  gwas_result <- data.frame(SNP=as.character(meta_result_shared_1p_filter$rs_id),CHR=meta_result_shared_1p_filter$CHR,BP=meta_result_shared_1p_filter$position,
                            P=pvalue,stringsAsFactors =F)
 # idx.cut <- which((gwas_result$BP%in%extract.list.ld$position)&(gwas_result$CHR%in%extract.list.ld$CHR))
  #gwas_result <- gwas_result[-idx.cut,]
  
  png(paste0("./result/plot/man_filter_Montse.png"),width = 7.635,height =4.7175,units = "in",res = 600)
  manhattan(gwas_result,suggestiveline = F, cex.axis = 2,main=paste0(" Global Test Association Manhattan Plot"),ylim=c(0,20))
  dev.off()
  
}













