load(paste0("/data/zhangh24/breast_cancer_data_analysis/whole_genome_age/ICOG/Intrinsic_subtypes/result/meta_result_shared_1p_082119.Rdata"))
var_name = meta_result_shared_1p[,"var_name",drop=F]
p_heter_list1 = list()
p_heter_list2 = list()
library(dplyr)
library(data.table)
for(i1 in 1:1000){
  load(paste0("/data/zhangh24/breast_cancer_data_analysis/whole_genome_age/ICOG/Intrinsic_subtypes/result/p_heter_sub_",i1,".rdata"))
  p_heter_list1[[i1]] = as.data.frame(p_list[[1]])
  p_heter_list2[[i1]] = as.data.frame(p_list[[2]])
}
p_heter = rbindlist(p_heter_list1)
p_heter_LA_to_TN = rbindlist(p_heter_list2)
heter_result = cbind(var_name,p_heter,p_heter_LA_to_TN)
snp_list = fread("/data/zhangh24/breast_cancer_data_analysis/whole_genome_age/ICOG/Intrinsic_subtypes/result/SNP_TWAS_GENES.txt",header= F)
colnames(snp_list) = "var_name"
twas_heter_result = left_join(snp_list, heter_result)
colnames(twas_heter_result)[2:3] = c("p_heter_global","p_heter_LAtoTN")
write.table(twas_heter_result, file = "/data/zhangh24/breast_cancer_data_analysis/whole_genome_age/ICOG/Intrinsic_subtypes/result/SNP_TWAS_GENES_heter.txt",
            row.names = F, col.names = T, quote = F)
twas_heter_result %>% filter(var_name == "19_17354586_C_A")
