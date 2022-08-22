#Goal: combine the subfiles of iCOGs into one big file
Filesdir <- "/data/NC_BW/icogs_onco/genotype/imputed2/icogs_imputed"
Files <- dir(Filesdir,pattern="icogs_merged_b1_12.",full.names=T)
# Filesex <- dir(Filesdir,pattern="icogs_merged_b1_12.chr23",full.names=T)
# idx.sex <- Files%in%Filesex
# Files <- Files[!idx.sex]
library(gtools)
Files <- mixedsort(Files)
#clean file names to only chromosmes and position
Files <- gsub("/data/NC_BW/icogs_onco/genotype/imputed2/icogs_imputed/icogs_merged_b1_12.","",Files)
Files <- gsub(".txt.gz","",Files)
n_files = length(Files)
#we want to order the files by chromosomes and positions
#mixedsort of files names couldn't achieve the goals
#therefore we manually order the files with the following command
#chr is the chromosme; p1 is the starting position number of a file; p2 is the ending position number of a file
Files_sub <- data.frame(chr=rep(1,length(Files)),p1=rep(0,length(Files)),p2=rep(0,length(Files)))

for(i in 1:length(Files)){
  temp <- gsub("/data/NC_BW/icogs_onco/genotype/imputed2/icogs_imputed/icogs_merged_b1_12.","",Files[i])
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
#idx gives you the order to load the files
idx <- order(Files_sub$chr,Files_sub$p1)
library(data.table)
rs_id_list = list()
score_list = list()
infor_list = list()
freq_list = list()

setwd("/data/zhangh24/breast_cancer_data_analysis/whole_genome_age/ICOG/ERPRHER2GRADE_fixed_baseline/result/")

for(i in 1:length(Files)){
  print(i)
  #load the idx[i] file
  load(paste0("ERPRHER2Grade_fixed_baseline",idx[i]))
  rs_id_list[[i]] = data.frame(rs_id = result[[1]])
  score_list[[i]] = as.data.frame(result[[2]])
  infor_list[[i]] = as.data.frame(result[[3]])
  freq_list[[i]] = data.frame(freq = result[[4]])
  
}
rs_id = rbindlist(rs_id_list)
score = rbindlist(score_list)
infor = rbindlist(infor_list)
freq.all = rbindlist(freq_list)

#load("/data/zhangh24/breast_cancer_data_analysis/whole_genome/ICOG/ERPRHER2GRADE_fixed_baseline/result/icog_info.Rdata")
#load the basic SNP information including ID, CHR, POS, imputation certainty, etc. 
load("/data/zhangh24/breast_cancer_data_analysis/whole_genome/ICOG/ERPRHER2GRADE_fixed_baseline/result/icog_info_clean.Rdata")
#test whether the rs_id in the two files are consistent
#the result of all.equal needs to be true to proceed
all.equal(rs_id$rs_id,icog_info$rs_id)
# icog_info <- cbind(icog_info,CHR)
# save(icog_info,file="/data/zhangh24/breast_cancer_data_analysis/whole_genome/ICOG/ERPRHER2GRADE_fixed_baseline/result/icog_info.Rdata")
CHR <- icog_info[,11]
icog_info <- icog_info[,1:10]

#combine the SNP basic information with score and infor matrix matrix
icog_result <- data.frame(icog_info,score,infor,CHR)


# CHR <- icog_result_baseline[,13]

#icog_result[,41] <- CHR


save(icog_result,file="/data/zhangh24/breast_cancer_data_analysis/whole_genome_age/ICOG/ERPRHER2GRADE_fixed_baseline/result/Icog_result.Rdata")
#icog_result_baseline <- data.frame(icog_info,score_baseline,infor_baseline,CHR)
#save(icog_result_baseline,file="/data/zhangh24/breast_cancer_data_analysis/whole_genome/ICOG/ERPRHER2GRADE_fixed_baseline/result/Icog_result_baseline.Rdata")
print(1)
#sicne sex chromsome result is seperately computed last time, we need to recompute it
# load("/data/zhangh24/breast_cancer_data_analysis/whole_genome_age/ICOG/ERPRHER2GRADE_fixed_baseline/result/Icog_result.Rdata")
# icog_result_auto = icog_result
# load("/data/zhangh24/breast_cancer_data_analysis/whole_genome_age/ICOG/ERPRHER2GRADE_fixed_baseline/result/Icog_result_sex.Rdata")
# icog_result_sex = icog_result
# icog_result = rbind(icog_result_auto,icog_result_sex)

