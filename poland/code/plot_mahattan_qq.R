setwd("/spin1/users/zhangh24/breast_cancer_data_analysis/poland/result/whole_genome/")
library(qqman)
#load("/spin1/users/zhangh24/breast_cancer_data_analysis/poland/result/whole_genome/onco_result_fixed_5p.Rdata")
#load("/spin1/users/zhangh24/breast_cancer_data_analysis/poland/result/whole_genome/onco_result_casecase_5p.Rdata")
load("/spin1/users/zhangh24/breast_cancer_data_analysis/poland/result/whole_genome/onco_result_standard_5p.Rdata")

#onco_result_fixed_5p
# gwas_result1 <- data.frame(SNP=as.character(onco_result_fixed_5p$rs_id),CHR=onco_result_fixed_5p$CHR,BP=onco_result_fixed_5p$position,
#                           P=onco_result_fixed_5p$P,stringsAsFactors =F)
# 
# png(paste0("./man_fix.png"),width = 7.635,height =4.7175,units = "in",res = 600)
# manhattan(gwas_result1,suggestiveline = F, cex.axis = 2,main=paste0("Global Test Association Manhattan Plot"))
# dev.off()
# png(paste0("./qq_fix.png"),width = 7.635,height =4.7175,units = "in",res = 600)
# qq(gwas_result1$P,main=paste0("Global Test Association Manhattan QQ Plot"))
# dev.off()
# gwas_result2 <- data.frame(SNP=as.character(onco_result_casecase_5p$rs_id),CHR=onco_result_casecase_5p$CHR,BP=onco_result_casecase_5p$position,
#                            P=onco_result_casecase_5p$P,stringsAsFactors =F)
# 
# png(paste0("./man_random.png"),width = 7.635,height =4.7175,units = "in",res = 600)
# manhattan(gwas_result2,suggestiveline = F, cex.axis = 2,main=paste0("Global Test Association Manhattan Plot"))
# dev.off()
# png(paste0("./qq_random.png"),width = 7.635,height =4.7175,units = "in",res = 600)
# qq(gwas_result2$P,main=paste0("Global Test Association Manhattan QQ Plot"))
# dev.off()

gwas_result3 <- data.frame(SNP=as.character(onco_result_standard_5p$rs_id),CHR=onco_result_standard_5p$CHR,BP=onco_result_standard_5p$position,
                           P=onco_result_standard_5p$p_value1,stringsAsFactors =F)
png(paste0("./man_standard.png"),width = 7.635,height =4.7175,units = "in",res = 600)
manhattan(gwas_result3,suggestiveline = F, cex.axis = 2,main=paste0("Global Test Association Manhattan Plot"))
dev.off()
png(paste0("./qq_standard.png"),width = 7.635,height =4.7175,units = "in",res = 600)
qq(gwas_result3$P,main=paste0("Global Test Association Manhattan QQ Plot"))
dev.off()

gwas_result4 <- data.frame(SNP=as.character(onco_result_standard_5p$rs_id),CHR=onco_result_standard_5p$CHR,BP=onco_result_standard_5p$position,
                           P=onco_result_standard_5p$p_value2,stringsAsFactors =F)
png(paste0("./man_poly.png"),width = 7.635,height =4.7175,units = "in",res = 600)
manhattan(gwas_result4,suggestiveline = F, cex.axis = 2,main=paste0("Global Test Association Manhattan Plot"))
dev.off()
png(paste0("./qq_poly.png"),width = 7.635,height =4.7175,units = "in",res = 600)
qq(gwas_result4$P,main=paste0("Global Test Association Manhattan QQ Plot"))
dev.off()










