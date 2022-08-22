rm(list=ls())

setwd("/data/shengf2/mtop/result/")
library(dplyr)


load("meta_result_shared_1p_final.Rdata")
ftop = meta_result_shared_1p
load("meta_result_shared_1p_pairwise.Rdata")
mtop.p = meta_result_shared_1p
load("meta_result_shared_1p_saturated.Rdata")
mtop.s = meta_result_shared_1p


all0 = inner_join(mtop.p, mtop.s, by = "var_name")
data = inner_join(all0, ftop,by = "var_name")

dim(data)
colnames(data)

all.equal(data$snp_id.x,data$snp_id.y)
all.equal(data$snp_id,data$snp_id.y)



## create a column for ACAT-pvalues
p.values = cbind(data$p.value,data$p.value.x,data$p.value.y)
p.values = t(p.values)
library(ACAT, lib.loc = '/home/shengf2/R/4.1/library/')
data$p.acat = ACAT(Pvals=p.values)


Known_snp = read.csv("210__known_discovery_snp_paper_order.csv",header = TRUE)
#suppose data is your final table with all SNPs information and ACAT p-values.
# idx <- which((data$CHR== Known_snp$CHR) & (data$position>=Known_snp$position-500000) & (data$position<=Known_snp$position+500000))

## i get the index by a for-loop
n.known = nrow(Known_snp)
idx = NULL
for(i in 1:n.known){
  CHRi = Known_snp$CHR[i]
  lb = Known_snp$position[i]-500000
  ub = Known_snp$position[i]+500000
  temp =  which((data$CHR==CHRi) & (data$position>=lb) & (data$position<=ub))
  idx = c(idx,temp)
  cat(i,length(idx),"\n")
}

idx = unique(idx)
length(idx)  # 648754

#data_clean contains the SNPs that are outside of the known SNPs region
data_clean = data[-idx,]

save(data,file = "data.rdata")
save(data_clean,file = "data_clean.rdata")


## significant variants  
id1 = which(data_clean$p.acat<5E-08)
length(id1)  # 505



