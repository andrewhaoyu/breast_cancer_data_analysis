rm(list=ls())

load("/data/NC_BW/HZ_SF/cond_results.rdata")
load("/data/NC_BW/HZ_SF/tensnps_info.rdata")
load("/data/NC_BW/HZ_SF/sixsnps_info.rdata")
library(tidyverse)

data = cbind(data_clean_select_1,p.acat1 = cond.res$p.acat)
head(data)
idx <- which(data[,13]<=1E-06)
data[idx,]

data_clean_select_2$p.acat1 = "c"
mydata = rbind(data[idx,],data_clean_select_2)


#create files for https://my.locuszoom.org/
library(dplyr)
library(tidyverse)

load("/data/NC_BW/HZ_SF/data3.rdata")
head(data_select)
head(data_select)

#this way you use CHR position as reference
n.snp = nrow(mydata)
for (i in 1:n.snp){
  CHRi = mydata$CHR[i]
  posi = mydata$position[i]
  
  test_data = data_select %>% 
    filter(CHR==CHRi&(position-posi)<=500000&(position-posi)>=-500000) %>% 
    separate(var_name,into = c("CHR_backup","position_backup","ref.allele","alt.allele"),
             remove = F, sep = "_") %>% 
    select(CHR, position, ref.allele, alt.allele, p.acat)
  
  #   rename(MarkerName = var_name, 
  #          P.value = p.acat)
  idx_order = order(test_data$CHR, test_data$position)
  test_data = test_data[idx_order,]
  colnames(test_data)[5] = "p-value"
  filename = paste0("/data/NC_BW/HZ_SF/test_data_locus_zoom_",i,".txt")
  write.table(test_data, file = filename, row.names = F, col.names = T, quote = F,
              sep = "\t")
  cat("SNP",i,"is done!\n")
}


load("/data/NC_BW/HZ_SF/snp_infor_match.rdata")
head(snp_infor_match)
snp_infor_match_sub = snp_infor_match %>% 
  select(chr.pos,rs_id)
#this way you create back up file using markername as reference
n.snp = nrow(mydata)
for (i in 1:n.snp){
  CHRi = mydata$CHR[i]
  posi = mydata$position[i]
  
  test_data = data_select %>% 
    filter(CHR==CHRi&(position-posi)<=500000&(position-posi)>=-500000) %>% 
    separate(var_name,into = c("CHR_backup","position_backup","ref.allele","alt.allele"),
             remove = F, sep = "_") %>% 
    mutate(chr.pos = paste0(CHR,":",position)) %>% 
    inner_join(y = snp_infor_match_sub, by = "chr.pos") %>% 
    rename(Markername = rs_id) %>% 
    select(Markername,CHR, position, ref.allele, alt.allele, p.acat)
  
  #   rename(MarkerName = var_name, 
  #          P.value = p.acat)
  idx_order = order(test_data$CHR, test_data$position)
  test_data = test_data[idx_order,]
  colnames(test_data)[6] = "p-value"
  filename = paste0("/data/NC_BW/HZ_SF/test_alt_data_locus_zoom_",i,".txt")
  write.table(test_data, file = filename, row.names = F, col.names = T, quote = F,
              sep = "\t")
  cat("SNP",i,"is done!\n")
}






# 
# 
# load("/data/NC_BW/HZ_SF/cond_results.rdata")
# load("/data/NC_BW/HZ_SF/tensnps_info.rdata")
# 
# setwd("")
# 
# 
# 
# data = cbind(data_clean_select_1,cond.res$p.acat)
# head(data)
# idx <- which(data[,13]<=1E-06)
# data[idx,]
# #create files for https://my.locuszoom.org/
# library(dplyr)
# library(tidyverse)
# load("/data/NC_BW/HZ_SF/data3.rdata")
# head(data_select)
# head(data_select)
# test_data = data_select %>% 
#   filter(CHR==1&(position-200342046)<=500000&(position-200342046)>=-500000) %>% 
#   separate(var_name,into = c("CHR_backup","position_backup","ref.allele","alt.allele"),
#            remove = F, sep = "_") %>% 
#   select(CHR, position, ref.allele, alt.allele, p.acat)
# 
# #   rename(MarkerName = var_name, 
# #          P.value = p.acat)
# idx_order = order(test_data$CHR, test_data$position)
# test_data = test_data[idx_order,]
# colnames(test_data)[5] = "p-value"
# write.table(test_data, file = "/data/NC_BW/HZ_SF/test_data_locus_zoom.txt", row.names = F, col.names = T, quote = F,
#             sep = "\t")