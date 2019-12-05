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
result.dir <- "/data/zhangh24/breast_cancer_data_analysis/risk_prediction/standard_whole_genome/ONCO/result"








setwd('/data/zhangh24/breast_cancer_data_analysis/risk_prediction/standard_whole_genome/ONCO/result')


num.total <- 0
size <- 6
for(i in 1:length(Files)){
  print(i)
  for(j in 1:size){
    
    
    load(paste0("/data/zhangh24/breast_cancer_data_analysis/risk_prediction/standard_whole_genome/ONCO/result/standard",idx[i],"_",j))
    temp <- length(result[[1]])
    num.total <- num.total+temp  
  }
  
}

num <- num.total
num.total <- 0
rs_id <- rep("c",num)
number.of.tumor <- 0
score <- matrix(0,nrow=num,ncol = (number.of.tumor+1))
infor <- matrix(0,nrow = num,ncol = (number.of.tumor+1)^2)
freq.all <- rep(0,num)
# score_baseline <- rep(0,num)
# infor_baseline <- rep(0,num)




setwd('/data/zhangh24/breast_cancer_data_analysis/risk_prediction/standard_whole_genome/ONCO/result')

num.total <- 0
num.length <- rep(0,length(Files))

###i is in the right order: ordered by chr and position
###idx[i] is the default files order got by the dir function in r

for(i in 1:length(Files)){
  print(i)
  for(j in 1:size){
  
    
    load(paste0("/data/zhangh24/breast_cancer_data_analysis/risk_prediction/standard_whole_genome/ONCO/result/standard",idx[i],"_",j))
    
    temp <- length(result[[1]])
    rs_id[num.total+(1:temp)] <- result[[1]]
    score[num.total+(1:temp),] <- result[[2]]
    infor[num.total+(1:temp),] <- result[[3]]
    freq.all[num.total+(1:temp)] <- result[[4]] 
    num.length[i] <- length(result[[1]])
    num.total <- temp+num.total
    
    
    
  }
  
}








load("/data/zhangh24/breast_cancer_data_analysis/whole_genome/ONCO/ERPRHER2GRADE_fixed_baseline/result/onco_info.Rdata")
all.equal(onco_info$rs_id,rs_id)
CHR <- onco_info[,11]
onco_info <- onco_info[,1:10]
onco_result <- data.frame(onco_info,score,infor,CHR)

save(onco_result,file="./onco_result.Rdata")
