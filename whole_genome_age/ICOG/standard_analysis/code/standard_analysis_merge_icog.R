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
result.dir <- "/spin1/users/zhangh24/breast_cancer_data_analysis/whole_genome_age/ICOG/standard_analysis/result/"


result_Files <- dir(result.dir,pattern="standard_analysis_")
result_Files <- result_Files[1:564]
# result.idx <- rep(0,length(result_Files))
# for(i in 1:length(result_Files)){
#   result.idx.temp <- as.integer(gsub("standard_analysis_","",result_Files[i]))
#   result.idx[i] <- result.idx.temp
# }


load("/spin1/users/zhangh24/breast_cancer_data_analysis/whole_genome_age/ICOG/ERPRHER2GRADE_fixed_baseline/result/Icog_result.Rdata")
rs_id <- icog_result$rs_id
num <- nrow(icog_result)

# num.total <- 0
# for(i in 1:564){
#   print(i)
#  
# }

#rs_id <- rep("c",num)
number.of.tumor <- 2
score <- matrix(0,nrow=num,ncol = (number.of.tumor-1))
infor <- matrix(0,nrow = num,ncol = (number.of.tumor-1)^2)
freq.all <- rep(0,num)

setwd("/spin1/users/zhangh24/breast_cancer_data_analysis/")


num.total <- 0
for(i in 1:length(Files)){
  print(i)
  for(k in 1:3){
    load(paste0("./whole_genome_age/ICOG/standard_analysis/result/standard_analysis_",idx[i],"_",k))
    temp <- nrow(result[[2]])
    score[num.total+(1:temp),] <- result[[2]]
    infor[num.total+(1:temp),] <- result[[3]]
    # for(j in 1:temp){
    #   infor_j <- result[[3]][(number.of.tumor*j-(number.of.tumor-1)):((number.of.tumor)*j),]
    #   infor[num.total+j,] <- as.vector(infor_j)
    # }
    # if(num.total< 12327300&(num.total+temp)> 12327300){
    #   print(c(i,k))
    # }
    num.total <- temp+num.total
    
    
  }
  
  
  
  
}




load("/spin1/users/zhangh24/breast_cancer_data_analysis/whole_genome/ICOG/ERPRHER2GRADE_fixed_baseline/result/icog_info.Rdata")
# icog_info <- cbind(icog_info,CHR)
# save(icog_info,file="/spin1/users/zhangh24/breast_cancer_data_analysis/whole_genome/ICOG/ERPRHER2GRADE_fixed_baseline/result/icog_info.Rdata")
all.equal(icog_info$rs_id,rs_id)
CHR <- icog_info[,11]
icog_info <- icog_info[,1:10]

icog_result_casecase <- data.frame(icog_info,score,infor,CHR)






save(icog_result_casecase,file="./whole_genome_age/ICOG/standard_analysis/result/icog_result_standard.Rdata")
# icog_result_baseline <- data.frame(icog_info,score_baseline,infor_baseline,CHR)
# save(icog_result_baseline,file="/spin1/users/zhangh24/breast_cancer_data_analysis/whole_genome/ICOG/ERPRHER2GRADE_fixed_baseline/result/Icog_result_baseline.Rdata")
# print(1)



