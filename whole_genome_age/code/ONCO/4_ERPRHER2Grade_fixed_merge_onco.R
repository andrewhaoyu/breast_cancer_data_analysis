Filesdir <- "/gpfs/gsfs4/users/NC_BW/icogs_onco/genotype/imputed2/onco_imputed"
# Files <- dir(Filesdir,pattern="OncoArray_european_merged_b1_15.",full.names=T)
# Filesex <- dir(Filesdir,pattern="OncoArray_european_merged_b1_15.chr23",full.names=T)
# idx.sex <- Files%in%Filesex
# Files <- Files[!idx.sex]
library(gtools)
Files <- mixedsort(Files)

Files <- gsub("/data/NC_BW/icogs_onco/genotype/imputed2/onco_imputed/OncoArray_european_merged_b1_15.","",Files)
Files <- gsub(".txt.gz","",Files)

Files_sub <- data.frame(chr=rep(1,length(Files)),p1=rep(0,length(Files)),p2=rep(0,length(Files)))

for(i in 1:length(Files)){
  temp <- gsub("/data/NC_BW/icogs_onco/genotype/imputed2/onco_imputed/OncoArray_european_merged_b1_15.","",Files[i])
  temp <- strsplit(temp,"\\.")
  temp <- unlist(temp)
  chr = as.integer(gsub("chr","",temp[1]))
  p_temp <- temp[2]
  p_temp <- strsplit(p_temp,"_")
  p_temp <- unlist(p_temp)
  p1 <- as.integer(p_temp[1])
  p2 <- as.integer(p_temp[2])
  Files_sub[i,] <- c(chr,p1,p2)
}
idx <- order(Files_sub$chr,Files_sub$p1)



setwd("/data/zhangh24/breast_cancer_data_analysis/whole_genome_age/ONCO/ERPRHER2GRADE_fixed_baseline/result/")
library(data.table)
rs_id_list = list()
score_list = list()
infor_list = list()
freq_list = list()

setwd("/data/zhangh24/breast_cancer_data_analysis/whole_genome_age/ICOG/ERPRHER2GRADE_fixed_baseline/result/")

for(i in 1:length(Files)){
  print(i)
  #load the idx[i] file
  load(paste0("ERPRHER2Grade_fixed_onco",idx[i]))
  rs_id_list[[i]] = data.frame(rs_id = result[[1]])
  score_list[[i]] = as.data.frame(result[[2]])
  infor_list[[i]] = as.data.frame(result[[3]])
  freq_list[[i]] = data.frame(freq = result[[4]])
  
}
rs_id = rbindlist(rs_id_list)
score = rbindlist(score_list)
infor = rbindlist(infor_list)
freq.all = rbindlist(freq_list)

load("/data/zhangh24/breast_cancer_data_analysis/whole_genome/ONCO/ERPRHER2GRADE_fixed_baseline/result/onco_info_clean.Rdata")
all.equal(onco_info$rs_id,rs_id)
CHR <- onco_info[,11]
onco_info <- onco_info[,1:10]


onco_result <- data.frame(onco_info,score,infor,CHR)

save(onco_result,file="/data/zhangh24/breast_cancer_data_analysis/whole_genome_age/ONCO/ERPRHER2GRADE_fixed_baseline/result/onco_result.Rdata")
print(1)
# onco_result_baseline <- data.frame(onco_info,score_baseline,infor_baseline,CHR)
# save(onco_result_baseline,file="/data/zhangh24/breast_cancer_data_analysis/wholge_genome_age/ONCO/ERPRHER2GRADE_fixed_baseline/result/onco_result_baseline.Rdata")
# print(2)

