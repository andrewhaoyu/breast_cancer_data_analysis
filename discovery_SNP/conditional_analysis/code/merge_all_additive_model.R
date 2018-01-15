result.com <- NULL

for(i1 in 1:207){
  print(i1)
  load(paste0("/spin1/users/zhangh24/breast_cancer_data_analysis/discovery_SNP/conditional_analysis/result/result.all",i1,".Rdata"))
result.com <- rbind(result.com,result.all)  
}


write.csv(result.com,file="/spin1/users/zhangh24/breast_cancer_data_analysis/discovery_SNP/conditional_analysis/result/condition.all.csv")
colnames(result.com)[19]="known_flag"
table(result.com$mark)
