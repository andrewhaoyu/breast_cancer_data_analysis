Filesdir <- "/gpfs/gsfs4/users/NC_BW/icogs_onco/genotype/imputed2/onco_imputed"
Files <- dir(Filesdir,pattern="OncoArray_european_merged_b1_15.",full.names=T)
Filesex <- dir(Filesdir,pattern="OncoArray_european_merged_b1_15.chr23",full.names=T)
idx.sex <- Files%in%Filesex
Files <- Files[!idx.sex]
library(gtools)
Files <- mixedsort(Files)

Files <- gsub("/gpfs/gsfs4/users/NC_BW/icogs_onco/genotype/imputed2/onco_imputed/OncoArray_european_merged_b1_15.","",Files)
Files <- gsub(".txt.gz","",Files)

Files_sub <- data.frame(chr=rep(1,length(Files)),p1=rep(0,length(Files)),p2=rep(0,length(Files)))

for(i in 1:length(Files)){
  temp <- gsub("/gpfs/gsfs4/users/NC_BW/icogs_onco/genotype/imputed2/onco_imputed/OncoArray_european_merged_b1_15.","",Files[i])
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
File_sub_order <- Files_sub[order(Files_sub$chr,Files_sub$p1),]
result.dir <- "/data/zhangh24/breast_cancer_data_analysis/risk_prediction/FTOP_whole_genome/ONCO/result"
result_Files <- dir(result.dir,pattern="ERPRHER2Grade_fixed_onco")
result_Files <- result_Files[1:567]
result.idx <- rep(0,length(result_Files))
for(i in 1:length(result_Files)){
  result_Files[i] <- gsub("ERPRHER2Grade_fixed_onco","",result_Files[i])
  result.idx.temp <- as.integer(gsub(".Rdata","",result_Files[i]))
  result.idx[i] <- result.idx.temp
}











setwd("/data/zhangh24/breast_cancer_data_analysis/risk_prediction/FTOP_whole_genome/ONCO/result")
num.total <- 0
for(i in 1:length(Files)){
  print(i)

  
  load(paste0("ERPRHER2Grade_fixed_onco",idx[i]))
  temp <- length(result[[1]])
  num.total <- num.total+temp
}

num <- num.total
num.total <- 0
rs_id <- rep("c",num)
number.of.tumor <- 4
score <- matrix(0,nrow=num,ncol = (number.of.tumor+1))
infor <- matrix(0,nrow = num,ncol = (number.of.tumor+1)^2)
freq.all <- rep(0,num)





setwd("/data/zhangh24/breast_cancer_data_analysis/risk_prediction/FTOP_whole_genome/ONCO/result")
num.total <- 0
num.length <- rep(0,length(Files))

###i is in the right order: ordered by chr and position
###idx[i] is the default files order got by the dir function in r

for(i in 1:length(Files)){
  print(i)
  
  load(paste0("ERPRHER2Grade_fixed_onco",idx[i]))
  
  temp <- length(result[[1]])
  #print(paste0("temp:",temp))
  rs_id[num.total+(1:temp)] <- result[[1]]
  score[num.total+(1:temp),] <- result[[2]]
  infor[num.total+(1:temp),] <- result[[3]]
  freq.all[num.total+(1:temp)] <- result[[4]] 
  
  num.length[i] <- length(result[[1]])
  num.total <- temp+num.total
  
  
  
}






load("/data/zhangh24/breast_cancer_data_analysis/risk_prediction/standard_whole_genome/ONCO/result/onco_result.Rdata")
CHR <- onco_result[,13]
onco_info <- onco_result[,1:10]
onco_result <- data.frame(onco_info,score,infor,CHR)

save(onco_result,file="/data/zhangh24/breast_cancer_data_analysis/risk_prediction/FTOP_whole_genome/ONCO/result/onco_result.Rdata")
# onco_result_baseline <- data.frame(onco_info,score_baseline,infor_baseline,CHR)
# save(onco_result_baseline,file="/data/zhangh24/breast_cancer_data_analysis/wholge_genome_age/ONCO/ERPRHER2GRADE_fixed_baseline/result/onco_result_baseline.Rdata")
# print(2)

