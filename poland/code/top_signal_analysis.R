load("/spin1/users/zhangh24/breast_cancer_data_analysis/poland/result/whole_genome/onco_result_casecase_5p.Rdata")
idx <-which(onco_result_casecase_5p$P<=5e-08)
sig <- onco_result_casecase_5p[idx,]
