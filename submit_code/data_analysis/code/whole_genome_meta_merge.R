
load("/data/zhangh24/breast_cancer_data_analysis/poland/result/whole_genome/onco_result_fixed_5p.Rdata")
load("/data/zhangh24/breast_cancer_data_analysis/poland/result/whole_genome/onco_result_casecase_5p.Rdata")
n <- nrow(onco_result_fixed_5p)
p_value1 <- rep(0,n)
p_value2 <- rep(0,n)
total <- 0
for(i1 in 1:1000){
  print(i1)
  load(paste0("/data/zhangh24/breast_cancer_data_analysis/poland/result/whole_genome/p_sub",i1,".Rdata"))
temp <- length(p_sub[[1]])
p_value1[total+(1:temp)] <- p_sub[[1]]
p_value2[total+(1:temp)] <- p_sub[[2]]
total = total+temp

}
onco_result_fixed_5p <- cbind(onco_result_fixed_5p,p_value1)
colnames(onco_result_fixed_5p)[42] <- "P"
onco_result_casecase_5p <- cbind(onco_result_casecase_5p,p_value2)
colnames(onco_result_casecase_5p)[24] <- "P"
save(onco_result_fixed_5p,file="/data/zhangh24/breast_cancer_data_analysis/poland/result/whole_genome/onco_result_fixed_5p.Rdata")
save(onco_result_casecase_5p,file="/data/zhangh24/breast_cancer_data_analysis/poland/result/whole_genome/onco_result_casecase_5p.Rdata")
