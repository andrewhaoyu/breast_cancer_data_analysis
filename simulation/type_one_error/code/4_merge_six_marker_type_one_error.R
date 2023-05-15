setwd("/data/zhangh24/breast_cancer_data_analysis/")
library(data.table)
result_list = list()
for(i1 in 1:1000){
  load(paste0("./simulation/type_one_error/result/topo_result/topo_result_",i1+ 10^6,".rdata"))
  result_list[[i1]] = as.data.frame(result)
}
result = rbindlist(result_list)
cor(result)

