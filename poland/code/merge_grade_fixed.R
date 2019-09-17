#merge the score test results for er and grade as fixed
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
# result.dir <- "/spin1/users/zhangh24/breast_cancer_data_analysis/whole_genome_age/ONCO/ERPRHER2GRADE_fixed_baseline/result/"
# result_Files <- dir(result.dir,pattern="ERPRHER2Grade_fixed_onco")
# result_Files <- result_Files[1:567]
# result.idx <- rep(0,length(result_Files))
# for(i in 1:length(result_Files)){
#   result_Files[i] <- gsub("ERPRHER2Grade_fixed_onco","",result_Files[i])
#   result.idx.temp <- as.integer(gsub(".Rdata","",result_Files[i]))
#   result.idx[i] <- result.idx.temp
# }


#load("/spin1/users/zhangh24/breast_cancer_data_analysis/whole_genome/ONCO/ERPRHER2GRADE_fixed_baseline/result/onco_result.Rdata")
load("/spin1/users/zhangh24/breast_cancer_data_analysis/whole_genome_age/ONCO/ERPRHER2GRADE_fixed_baseline/result/onco_result.Rdata")
num.total <- nrow(onco_result)



num <- num.total
#num.total <- 0
rs_id <- rep("c",num.total)
number.of.tumor <- 4
score <- matrix(0,nrow=num,ncol = (number.of.tumor+1))
infor <- matrix(0,nrow = num,ncol = (number.of.tumor+1)^2)
score2 <- matrix(0,nrow=num,ncol = (number.of.tumor-2))
infor2 <- matrix(0,nrow = num,ncol = (number.of.tumor-2)^2)

freq.all <- rep(0,num)
#score_baseline <- rep(0,num)
#infor_baseline <- rep(0,num)


job.sub.length <- rep(0,567)

setwd("/spin1/users/zhangh24/breast_cancer_data_analysis/")
num.total <- 0
for(i in 1:567){
  print(i)
  for (k in 1:5) {
    #print(k)
   
    
    load(
      paste0("./poland/result/whole_genome/grade_fixed_",idx[i],"_",k)
    )
    temp <- length(result[[2]])
    rs_id[num.total+(1:temp)] <- result[[1]]
    
    freq.all[num.total+(1:temp)] <- result[[2]] 
    score2[num.total+(1:temp),] <- result[[3]]
    infor2[num.total+(1:temp),] <- result[[4]]
    # 
    num.total <- temp+num.total
    job.sub.length[i] <- job.sub.length[i]+temp
  }  
  
  
}




# 
# 
# num.length.info <- rep(0,length(Files))
# 
# onco_info <- data.frame(snp_id = rep("c",num.total),rs_id = rep("c",num.total),
#                         position=rep(0,num.total),exp_freq_a1=rep(0,num.total),info=rep(0,num.total),
#                         certainty=rep(0,num.total),type=rep(0,num.total),info_type0=rep(0,num.total),
#                         concord_type0=rep(0,num.total),r2_type0=rep(0,num.total),stringsAsFactors=F)
# CHR <- rep(0,num.total)
# num.total <-  0
# temp.j <- 0
# for(i in 1:22){
#   print(i)
#   filedir <- paste0("/gpfs/gsfs4/users/NC_BW/icogs_onco/genotype/imputed2/onco_info_files/chr",i)
#   files <- dir(filedir,pattern="txt_info",full.names=T)
#   files_num <- gsub(paste0(filedir,"/OncoArray_chr",i,"_euro15_phased_"),
#                     "",files)
#   files_num <- gsub(".txt_info","",files_num)
#   files_num <- strsplit(files_num,"_")
#   files_num <- as.integer(unlist(files_num)[seq(1,2*length(files_num),2)])
#   idx <- order(files_num)
#   for(j in 1:length(idx)){
#     temp.j <- temp.j +1
#     print(temp.j)
#     #print(j)
#     data <- read.table(files[idx[j]],header=T,stringsAsFactors=F)
#     temp <- nrow(data)
#     num.length.info[temp.j] <- temp
#     onco_info[num.total+(1:temp),1:3] <- data[,1:3]
#     onco_info[num.total+(1:temp),4:10] <- data[,6:12]
#     CHR[num.total+(1:temp)] <- i
#     num.total <- temp+num.total
# 
#   }
# 
#  }
# 
#  load("/spin1/users/zhangh24/breast_cancer_data_analysis/whole_genome/ONCO/ERPRHER2GRADE_fixed_baseline/result/onco_result.Rdata")
# 
# 
#  onco_info <- onco_result[,1:10]
#  CHR <- onco_result[,41]
#  onco_info <- cbind(onco_info,CHR)
#  save(onco_info,file="/spin1/users/zhangh24/breast_cancer_data_analysis/whole_genome/ONCO/ERPRHER2GRADE_fixed_baseline/result/onco_info.Rdata")




load("/spin1/users/zhangh24/breast_cancer_data_analysis/whole_genome/ONCO/ERPRHER2GRADE_fixed_baseline/result/onco_info.Rdata")

all.equal(onco_info$rs_id,rs_id)

CHR <- onco_info[,11]
onco_info <- onco_info[,1:10]


onco_result_fixed <- data.frame(onco_info,score,infor,CHR)
onco_result_casecase <- data.frame(onco_info,score2,infor2,CHR)


#save(onco_result_fixed,file="/spin1/users/zhangh24/breast_cancer_data_analysis/poland/result/whole_genome/onco_result_fixed.Rdata")
save(onco_result_casecase_grade,file="/spin1/users/zhangh24/breast_cancer_data_analysis/poland/result/whole_genome/onco_result_casecase_grade.Rdata")

idx <- which(onco_result_fixed$exp_freq_a1>=0.05&
               onco_result_fixed$exp_freq_a1<=0.95)
onco_result_fixed_5p <- onco_result_fixed[idx,]
onco_result_casecase_5p <- onco_result_casecase[idx,]
save(onco_result_fixed_5p,file="/spin1/users/zhangh24/breast_cancer_data_analysis/poland/result/whole_genome/onco_result_fixed_5p.Rdata")
save(onco_result_casecase_5p,file="/spin1/users/zhangh24/breast_cancer_data_analysis/poland/result/whole_genome/onco_result_casecase_5p.Rdata")

print(1)
