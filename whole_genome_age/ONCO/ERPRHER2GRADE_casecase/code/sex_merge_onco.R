Filesdir <- "/gpfs/gsfs4/users/NC_BW/icogs_onco/genotype/imputed2/onco_imputed"
Files <- dir(Filesdir,pattern="OncoArray_european_merged_b1_15.",full.names=T)
Filesex <- dir(Filesdir,pattern="OncoArray_european_merged_b1_15.chr23",full.names=T)
idx.sex <- Files%in%Filesex
Files <- Files[idx.sex]
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
# result.dir <- "/data/zhangh24/breast_cancer_data_analysis/whole_genome_age/ONCO/ERPRHER2GRADE_fixed_baseline/result/"
# result_Files <- dir(result.dir,pattern="ERPRHER2Grade_fixed_onco")
# result_Files <- result_Files[1:567]
# result.idx <- rep(0,length(result_Files))
# for(i in 1:length(result_Files)){
#   result_Files[i] <- gsub("ERPRHER2Grade_fixed_onco","",result_Files[i])
#   result.idx.temp <- as.integer(gsub(".Rdata","",result_Files[i]))
#   result.idx[i] <- result.idx.temp
# }


#load("/data/zhangh24/breast_cancer_data_analysis/whole_genome/ONCO/ERPRHER2GRADE_fixed_baseline/result/onco_result.Rdata")
load("/data/zhangh24/breast_cancer_data_analysis/whole_genome_age/ONCO/ERPRHER2GRADE_fixed_baseline/result/onco_result_sex.Rdata")
num.total <- nrow(onco_result)



num <- num.total
#num.total <- 0
rs_id <- rep("c",num.total)
number.of.tumor <- 4
score <- matrix(0,nrow=num,ncol = (number.of.tumor-1))
infor <- matrix(0,nrow = num,ncol = (number.of.tumor-1)^2)
freq.all <- rep(0,num)
#score_baseline <- rep(0,num)
#infor_baseline <- rep(0,num)





num.total <- 0
for(i in 1:32){
  
    print(i)
    for (k in 1:30) {
      #print(k)
      load(
        paste0("/data/zhangh24/breast_cancer_data_analysis/whole_genome_age/ONCO/ERPRHER2GRADE_casecase/result/ERPRHER2Grade_casecase_onco_sex",idx[i],"_",k)
      )
      temp <- nrow(result[[2]])
      rs_id[num.total+(1:temp)] <- result[[1]]
      score[num.total+(1:temp),] <- result[[2]]
      infor[num.total+(1:temp),] <- result[[3]]
      freq.all[num.total+(1:temp)] <- result[[4]] 
      num.total <- temp+num.total
      
    }  
  
  
  
  
  
}







load("/data/zhangh24/breast_cancer_data_analysis/whole_genome_age/ONCO/ERPRHER2GRADE_fixed_baseline/result/onco_infor_sex.Rdata")

all.equal(onco_info$rs_id,rs_id)

CHR <- rep(23,nrow(onco_info))
onco_info <- onco_info[,1:10]


onco_result_casecase <- data.frame(onco_info,score,infor,CHR)

save(onco_result_casecase,file="/data/zhangh24/breast_cancer_data_analysis/whole_genome_age/ONCO/ERPRHER2GRADE_casecase/result/onco_result_casecase_sex.Rdata")
print(1)
