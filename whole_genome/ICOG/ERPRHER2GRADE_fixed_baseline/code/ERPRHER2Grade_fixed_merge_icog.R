Filesdir <- "/gpfs/gsfs4/users/NC_BW/icogs_onco/genotype/imputed2/icogs_imputed/"
Files <- dir(Filesdir,pattern="icogs_merged_b1_12.",full.names=T)
Filesex <- dir(Filesdir,pattern="icogs_merged_b1_12.chr23",full.names=T)
idx.sex <- Files%in%Filesex
Files <- Files[!idx.sex]
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
result.dir <- "/spin1/users/zhangh24/breast_cancer_data_analysis/whole_genome/ICOG/ERPRHER2GRADE_fixed_baseline/result/"

result_Files <- dir(result.dir,pattern="ERPRHER2Grade_fixed_baseline")
result_Files <- result_Files[1:564]
result.idx <- rep(0,length(result_Files))
for(i in 1:length(result_Files)){
  result.idx.temp <- as.integer(gsub("ERPRHER2Grade_fixed_baseline","",result_Files[i]))
  result.idx[i] <- result.idx.temp
}


num.total <- 0
for(i in 1:length(Files)){
  print(i)
  setwd("/spin1/users/zhangh24/breast_cancer_data_analysis/whole_genome/ICOG/ERPRHER2GRADE_fixed_baseline/result/")
  
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
score_baseline <- rep(0,num)
infor_baseline <- rep(0,num)



setwd("/spin1/users/zhangh24/breast_cancer_data_analysis/whole_genome/ICOG/ERPRHER2GRADE_fixed_baseline/result/")
num.total <- 0
for(i in 1:length(Files)){
  print(i)
  
  load(paste0("ERPRHER2Grade_fixed_baseline",idx[i]))
  temp <- length(result[[1]])
  rs_id[num.total+(1:temp)] <- result[[1]]
  score[num.total+(1:temp),] <- result[[2]]
  for(j in 1:temp){
    infor_j <- result[[3]][((number.of.tumor+1)*j-number.of.tumor):((number.of.tumor+1)*j),]
    infor[num.total+j,] <- as.vector(infor_j)
  }
  freq.all[num.total+(1:temp)] <- result[[4]] 
  score_baseline[num.total+(1:temp)] <- result[[5]]
  infor_baseline[num.total+(1:temp)] <- result[[6]]
  num.total <- temp+num.total
  
  
  
}












#####to get the total number of SNPs from the information files


icog_info <- data.frame(snp_id = rep("c",num.total),rs_id = rep("c",num.total),
                        position=rep(0,num.total),exp_freq_a1=rep(0,num.total),info=rep(0,num.total),
                        certainty=rep(0,num.total),type=rep(0,num.total),info_type0=rep(0,num.total),
                        concord_type0=rep(0,num.total),r2_type0=rep(0,num.total),stringsAsFactors=F)
CHR <- rep(0,num.total)
num.total <-  0

for(i in 1:22){
  print(i)
  filedir <- paste0("/gpfs/gsfs4/users/NC_BW/icogs_onco/genotype/imputed2/icogs_info_files/chr",i)
  files <- dir(filedir,pattern="txt_info",full.names=T)
  files_num <- gsub(paste0(filedir,"/icogs_euro12_chr",i,"_phased"),
                    "",files)
  files_num <- gsub(".txt_info","",files_num)
  files_num <- strsplit(files_num,"_")
  files_num <- as.integer(unlist(files_num)[seq(1,2*length(files_num),2)])
  idx <- order(files_num)
  for(j in 1:length(idx)){
    #print(j)
    data <- read.table(files[idx[j]],header=T,stringsAsFactors=F)
    temp <- nrow(data)
    icog_info[num.total+(1:temp),] <- data
    CHR[num.total+(1:temp)] <- i
    num.total <- temp+num.total
  }
  
}


icog_result <- data.frame(icog_info,score,infor,CHR)







save(icog_result,file="/spin1/users/zhangh24/breast_cancer_data_analysis/whole_genome/ICOG/ERPRHER2GRADE_fixed_baseline/result/Icog_result.Rdata")
icog_result <- data.frame(icog_info,score_baseline,infor_baseline,CHR)
save(icog_result,file="/spin1/users/zhangh24/breast_cancer_data_analysis/whole_genome/ICOG/ERPRHER2GRADE_fixed_baseline/result/Icog_result_baseline.Rdata")
print(1)



