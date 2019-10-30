setwd("/data/zhangh24/breast_cancer_data_analysis/")
filedir <- './risk_prediction/intrinsic_subtypes_whole_genome/ICOG/result/'
#files <- dir(filedir,pattern="intrinsic_subytpe_icog_resubmit")
files <- dir(filedir,pattern="intrinsic_subytpe_icog")
result_files <- dir(filedir,pattern="intrinsic_subytpe_icog")





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
result.dir <- './risk_prediction/intrinsic_subtypes_whole_genome/ICOG/result'



load("/data/zhangh24/breast_cancer_data_analysis/whole_genome_age/ICOG/ERPRHER2GRADE_fixed_baseline/result/Icog_result.Rdata")
#rs_id <- icog_result$rs_id
num <- nrow(icog_result)

# num.total <- 0
# for(i in 1:564){
#   print(i)
#  
# }

#rs_id <- rep("c",num)
number.of.tumor <- 4
score <- matrix(0,nrow=num,ncol = (number.of.tumor+1))
infor <- matrix(0,nrow = num,ncol = (number.of.tumor+1)^2)
freq.all <- rep(0,num)
rs_id <- rep("c",num)


resubimt_resubmimt_id <- c(231,281,486)

#resubmit_id <- matrix(0,100,2)
#resubmit_temp <- 0
num.total <- 0
for(i in 1:length(Files)){
  print(i)
  file_load = paste0("intrinsic_subytpe_icog_resubmit",idx[i],"_",1)
  if(idx[i]%in%resubimt_resubmimt_id){
    for(k in 1:750){
      load(paste0("./risk_prediction/intrinsic_subtypes_whole_genome/ICOG/result/intrinsic_subytpe_icog_resubmit",idx[i],"_",k))
      temp <- nrow(result[[2]])
      rs_id[num.total+(1:temp)] <- result[[1]]
      score[num.total+(1:temp),] <- result[[2]]
      infor[num.total+(1:temp),] <- result[[3]]
      num.total <- temp+num.total
      if(sum(result[[1]]=="c")!=0){
        resubmit_temp <- resubmit_temp+1
        resubmit_id[resubmit_temp,1] <- idx[i]
        resubmit_id[resubmit_temp,2] <- k
      }
    }
  }else{
    for(k in 1:2){
      load(paste0("./risk_prediction/intrinsic_subtypes_whole_genome/ICOG/result/intrinsic_subytpe_icog",idx[i],"_",k))
      temp <- nrow(result[[2]])
      rs_id[num.total+(1:temp)] <- result[[1]]
      score[num.total+(1:temp),] <- result[[2]]
      infor[num.total+(1:temp),] <- result[[3]]
      num.total <- temp+num.total
    }
  } 
  }

load("/data/zhangh24/breast_cancer_data_analysis/whole_genome/ICOG/ERPRHER2GRADE_fixed_baseline/result/icog_info.Rdata")
# icog_info <- cbind(icog_info,CHR)
# save(icog_info,file="/data/zhangh24/breast_cancer_data_analysis/whole_genome/ICOG/ERPRHER2GRADE_fixed_baseline/result/icog_info.Rdata")
all.equal(icog_info$rs_id,rs_id)
#idx.diff <- which(icog_info$rs_id!=rs_id)
CHR <- icog_info[,11]
icog_info <- icog_info[,1:10]

icog_result_casecase <- data.frame(icog_info,score,infor,CHR)






save(icog_result_casecase,file="/data/zhangh24/breast_cancer_data_analysis/risk_prediction/intrinsic_subtypes_whole_genome/ICOG/result/Icog_result_intrinsic_subtype.Rdata")
# icog_result_baseline <- data.frame(icog_info,score_baseline,infor_baseline,CHR)
# save(icog_result_baseline,file="/data/zhangh24/breast_cancer_data_analysis/whole_genome/ICOG/ERPRHER2GRADE_fixed_baseline/result/Icog_result_baseline.Rdata")
# print(1)



