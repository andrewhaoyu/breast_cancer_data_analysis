Filesdir <- "/gpfs/gsfs4/users/NC_BW/icogs_onco/genotype/imputed2/onco_imputed"
Files <- dir(Filesdir,pattern="OncoArray_european_merged_b1_15.",full.names=T)
Filesex <- dir(Filesdir,pattern="OncoArray_european_merged_b1_15.chr23",full.names=T)
idx.sex <- Files%in%Filesex
Files <- Files[!idx.sex]

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
result.dir <- "/spin1/users/zhangh24/breast_cancer_data_analysis/risk_prediction/standard_whole_genome/ONCO/result"
result_Files <- dir(result.dir,pattern="standard")
result_Files <- result_Files[1:567]
result.idx <- rep(0,length(result_Files))
for(i in 1:length(result_Files)){
  result_Files[i] <- gsub("standard","",result_Files[i])
  result.idx.temp <- as.integer(gsub(".Rdata","",result_Files[i]))
  result.idx[i] <- result.idx.temp
}








setwd('/spin1/users/zhangh24/breast_cancer_data_analysis/risk_prediction/standard_whole_genome/ONCO/result')


num.total <- 0
for(i in 1:length(Files)){
  print(i)
 
  load(paste0("standard",idx[i]))
  temp <- length(result[[1]])
  num.total <- num.total+temp
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




setwd('/spin1/users/zhangh24/breast_cancer_data_analysis/risk_prediction/standard_whole_genome/ONCO/result')

num.total <- 0
num.length <- rep(0,length(Files))

###i is in the right order: ordered by chr and position
###idx[i] is the default files order got by the dir function in r

for(i in 1:length(Files)){
  print(i)
  
  load(paste0("standard",idx[i]))
  
  temp <- length(result[[1]])
  print(paste0("temp:",temp))
  rs_id[num.total+(1:temp)] <- result[[1]]
  score[num.total+(1:temp),] <- result[[2]]
  infor[num.total+(1:temp),] <- result[[3]]
  freq.all[num.total+(1:temp)] <- result[[4]] 
  num.length[i] <- length(result[[1]])
  num.total <- temp+num.total
  
  
  
}




num.length.info <- rep(0,length(Files))

onco_info <- data.frame(snp_id = rep("c",num.total),rs_id = rep("c",num.total),
                        position=rep(0,num.total),exp_freq_a1=rep(0,num.total),info=rep(0,num.total),
                        certainty=rep(0,num.total),type=rep(0,num.total),info_type0=rep(0,num.total),
                        concord_type0=rep(0,num.total),r2_type0=rep(0,num.total),stringsAsFactors=F)
CHR <- rep(0,num.total)
num.total <-  0
temp.j <- 0
for(i in 1:22){
  print(i)
  filedir <- paste0("/gpfs/gsfs4/users/NC_BW/icogs_onco/genotype/imputed2/onco_info_files/chr",i)
  files <- dir(filedir,pattern="txt_info",full.names=T)
  files_num <- gsub(paste0(filedir,"/OncoArray_chr",i,"_euro15_phased_"),
                    "",files)
  files_num <- gsub(".txt_info","",files_num)
  files_num <- strsplit(files_num,"_")
  files_num <- as.integer(unlist(files_num)[seq(1,2*length(files_num),2)])
  idx <- order(files_num)
  for(j in 1:length(idx)){
    temp.j <- temp.j +1
    print(temp.j)
    #print(j)
    data <- read.table(files[idx[j]],header=T,stringsAsFactors=F)
    temp <- nrow(data)
    num.length.info[temp.j] <- temp
    onco_info[num.total+(1:temp),1:3] <- data[,1:3]
    onco_info[num.total+(1:temp),4:10] <- data[,6:12]
    CHR[num.total+(1:temp)] <- i
    num.total <- temp+num.total
    
  }
  
}

onco_result <- data.frame(onco_info,score,infor,CHR)

save(onco_result,file="./onco_result.Rdata")
# print(1)
# onco_result_baseline <- data.frame(onco_info,score_baseline,infor_baseline,CHR)
# save(onco_result_baseline,file="/spin1/users/zhangh24/breast_cancer_data_analysis/whole_genome/ONCO/ERPRHER2GRADE_fixed_baseline/result/onco_result_baseline.Rdata")
# print(2)

