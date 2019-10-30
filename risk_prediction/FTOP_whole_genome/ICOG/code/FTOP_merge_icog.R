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
result.dir <- "/data/zhangh24/breast_cancer_data_analysis/risk_prediction/FTOP_whole_genome/ICOG/result"

result_Files <- dir(result.dir,pattern="ERPRHER2Grade_fixed_baseline")
result_Files <- result_Files[1:564]
result.idx <- rep(0,length(result_Files))
for(i in 1:length(result_Files)){
  result.idx.temp <- as.integer(gsub("ERPRHER2Grade_fixed_baseline","",result_Files[i]))
  result.idx[i] <- result.idx.temp
}

setwd("/data/zhangh24/breast_cancer_data_analysis/risk_prediction/FTOP_whole_genome/ICOG/result")
num.total <- 0
for(i in 1:length(Files)){
  print(i)
  
  load(paste0("ERPRHER2Grade_fixed_baseline",idx[i]))
  temp <- length(result[[1]])
  num.total <- num.total+temp
}

num <-  num.total

rs_id <- rep("c",num)
number.of.tumor <- 4
score <- matrix(0,nrow=num,ncol = (number.of.tumor+1))
infor <- matrix(0,nrow = num,ncol = (number.of.tumor+1)^2)
freq.all <- rep(0,num)





setwd("/data/zhangh24/breast_cancer_data_analysis/whole_genome_age/ICOG/ERPRHER2GRADE_fixed_baseline/result/")
num.total <- 0
for(i in 1:length(Files)){
  print(i)
  
  load(paste0("ERPRHER2Grade_fixed_baseline",idx[i]))
  temp <- length(result[[1]])
  rs_id[num.total+(1:temp)] <- result[[1]]
  score[num.total+(1:temp),] <- result[[2]]
  infor[num.total+(1:temp),] <- result[[3]]
  
  freq.all[num.total+(1:temp)] <- result[[4]] 
  num.total <- temp+num.total
  
  
  
}


load('/data/zhangh24/breast_cancer_data_analysis/risk_prediction/standard_whole_genome/ICOG/result/Icog_result.Rdata')
icog_infor <- icog_result[,1:10]
CHR <- icog_result[,13]


icog_result <- data.frame(icog_infor,score,infor,CHR)



# CHR <- icog_result_baseline[,13]

#icog_result[,41] <- CHR


save(icog_result,file="/data/zhangh24/breast_cancer_data_analysis/risk_prediction/FTOP_whole_genome/ICOG/result/Icog_result.Rdata")
#icog_result_baseline <- data.frame(icog_info,score_baseline,infor_baseline,CHR)
#save(icog_result_baseline,file="/data/zhangh24/breast_cancer_data_analysis/whole_genome/ICOG/ERPRHER2GRADE_fixed_baseline/result/Icog_result_baseline.Rdata")
print(1)



