library(qqman)
load("/spin1/users/zhangh24/breast_cancer_data_analysis/poland/result/whole_genome/onco_result_fixed_5p.Rdata")
load("/spin1/users/zhangh24/breast_cancer_data_analysis/poland/result/whole_genome/onco_result_casecase_5p.Rdata")
load("/spin1/users/zhangh24/breast_cancer_data_analysis/poland/result/whole_genome/onco_result_standard_5p.Rdata")
gwas_result1 <- data.frame(SNP=as.character(onco_result_fixed_5p$rs_id),CHR=onco_result_fixed_5p$CHR,BP=onco_result_fixed_5p$position,
                           P=onco_result_fixed_5p$P,stringsAsFactors =F)
gwas_result2 <- data.frame(SNP=as.character(onco_result_casecase_5p$rs_id),CHR=onco_result_casecase_5p$CHR,BP=onco_result_casecase_5p$position,
                           P=onco_result_casecase_5p$P,stringsAsFactors =F)
gwas_result3 <- data.frame(SNP=as.character(onco_result_standard_5p$rs_id),CHR=onco_result_standard_5p$CHR,BP=onco_result_standard_5p$position,
                           P=onco_result_standard_5p$p_value1,stringsAsFactors =F)
gwas_result4 <- data.frame(SNP=as.character(onco_result_standard_5p$rs_id),CHR=onco_result_standard_5p$CHR,BP=onco_result_standard_5p$position,
                           P=p_value2,stringsAsFactors =F)
idx <- which(gwas_result1$P<=5E-08)
temp <- gwas_result1[idx,]
temp[order(temp$P),]

idx <- which(gwas_result2$P<=5E-08)
temp <- gwas_result2[idx,]
temp[order(temp$P),]


idx <- which(gwas_result3$P<=5E-08)
temp <- gwas_result3[idx,]
temp[order(temp$P),]
