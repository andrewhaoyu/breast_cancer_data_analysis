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
  load("/spin1/users/zhangh24/breast_cancer_data_analysis/whole_genome/ICOG/ERPRHER2_fixed/result/meta_result_shared_1p_filter.Rdata")
  pvalue <- meta_result_shared_1p_filter[,15]
  idx <- which(pvalue<=1.0e-20)
  pvalue[idx] <- 1.0e-20
  gwas_result <- data.frame(SNP=as.character(meta_result_shared_1p_filter$rs_id),CHR=meta_result_shared_1p_filter$CHR,BP=meta_result_shared_1p_filter$position,
                            P=pvalue,stringsAsFactors =F)
  png(paste0("./result/plot/man_filter.png"),width = 7.635,height =4.7175,units = "in",res = 600)
  manhattan(gwas_result,suggestiveline = F, cex.axis = 2,main=paste0(" Global Test Association Manhattan Plot"),ylim=c(0,20))
  dev.off()
}else{
  load("/spin1/users/zhangh24/breast_cancer_data_analysis/whole_genome/ICOG/ERPRHER2_fixed/result/meta_result_shared_1p_filter.Rdata")
  new_filter <- read.csv("/spin1/users/zhangh24/breast_cancer_data_analysis/whole_genome/ICOG/ERPRHER2_fixed/result/Filter_based_on_Montse.csv",header=T,stringsAsFactors = F)
  new_filter[,2] <- as.numeric(gsub(",","",new_filter[,2]))
  
  pvalue <- meta_result_shared_1p_filter[,15]
  idx <- which(pvalue<=1.0e-20)
  pvalue[idx] <- 1.0e-20
  gwas_result <- data.frame(SNP=as.character(meta_result_shared_1p_filter$rs_id),CHR=meta_result_shared_1p_filter$CHR,BP=meta_result_shared_1p_filter$position,
                            P=pvalue,stringsAsFactors =F)
  idx_cut <- NULL
  
  
  for(i in 1:nrow(new_filter)){
    print(i)
    chr_temp <- new_filter[i,3]
    position_temp <- new_filter[i,2]
    position_low <- position_temp-0.5*10^6
    position_high <- position_temp+0.5*10^6
    idx <- which(meta_result_shared_1p_filter$CHR==chr_temp&meta_result_shared_1p_filter$position>position_low&
                   meta_result_shared_1p_filter$position<position_high)
    idx_cut <- c(idx_cut,idx)
  }
  idx_cut <- unique(idx_cut)
  gwas_result <- gwas_result[-idx_cut,]
  
  png(paste0("./result/plot/man_filter_Montse.png"),width = 7.635,height =4.7175,units = "in",res = 600)
  manhattan(gwas_result,suggestiveline = F, cex.axis = 2,main=paste0(" Global Test Association Manhattan Plot"),ylim=c(0,20))
  dev.off()
  
}













