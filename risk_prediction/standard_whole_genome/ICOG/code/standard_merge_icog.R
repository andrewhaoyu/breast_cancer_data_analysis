Filesdir <- "/gpfs/gsfs4/users/NC_BW/icogs_onco/genotype/imputed2/icogs_imputed/"
Files <- dir(Filesdir,pattern="icogs_merged_b1_12.",full.names=T)
Filesex <- dir(Filesdir,pattern="icogs_merged_b1_12.chr23",full.names=T)
idx.sex <- Files%in%Filesex
Files <- Files[!idx.sex]
library(gtools)
Files <- mixedsort(Files)

Files <- gsub("/gpfs/gsfs4/users/NC_BW/icogs_onco/genotype/imputed2/icogs_imputed/icogs_merged_b1_12.","",Files)
Files <- gsub(".txt.gz","",Files)


Files_sub <- data.frame(chr=rep(1,length(Files)),p1=rep(0,length(Files)),p2=rep(0,length(Files)))

for(i in 1:length(Files)){
  temp <- gsub("/gpfs/gsfs4/users/NC_BW/icogs_onco/genotype/imputed2/icogs_imputed//icogs_merged_b1_12.","",Files[i])
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
result.dir <- "/spin1/users/zhangh24/breast_cancer_data_analysis/risk_prediction/standard_whole_genome/ICOG/result"

result_Files <- dir(result.dir,pattern="standard")
result_Files <- result_Files[1:564]
result.idx <- rep(0,length(result_Files))
for(i in 1:length(result_Files)){
  result.idx.temp <- as.integer(gsub("standard","",result_Files[i]))
  result.idx[i] <- result.idx.temp
}

setwd("/spin1/users/zhangh24/breast_cancer_data_analysis/risk_prediction/standard_whole_genome/ICOG/result")
num.total <- 0
for(i in 1:length(Files)){
  print(i)
  
  
  load(paste0("standard",idx[i]))
  temp <- length(result[[1]])
  num.total <- num.total+temp
}

num <-  num.total
rs_id <- rep("c",num)
number.of.tumor <- 0
score <- matrix(0,nrow=num,ncol = (number.of.tumor+1))
infor <- matrix(0,nrow = num,ncol = (number.of.tumor+1)^2)
freq.all <- rep(0,num)




setwd("/spin1/users/zhangh24/breast_cancer_data_analysis/risk_prediction/standard_whole_genome/ICOG/result")
num.total <- 0
for(i in 1:length(Files)){
  print(i)
  
  load(paste0("standard",idx[i]))
  temp <- length(result[[1]])
  rs_id[num.total+(1:temp)] <- result[[1]]
  score[num.total+(1:temp),] <- result[[2]]
  infor[num.total+(1:temp),] <- result[[3]]
  freq.all[num.total+(1:temp)] <- result[[4]] 
  num.total <- temp+num.total
  
  
  
}

load("/spin1/users/zhangh24/breast_cancer_data_analysis/whole_genome/ICOG/ERPRHER2GRADE_fixed_baseline/result/icog_info.Rdata")
all.equal(rs_id,icog_info$rs_id)
CHR <- icog_info[,11]
icog_info <- icog_info[,1:10]
icog_result <- data.frame(icog_info,score,infor,CHR)




setwd('/spin1/users/zhangh24/breast_cancer_data_analysis/risk_prediction/standard_whole_genome/ICOG/result')


save(icog_result,file="./Icog_result.Rdata")




