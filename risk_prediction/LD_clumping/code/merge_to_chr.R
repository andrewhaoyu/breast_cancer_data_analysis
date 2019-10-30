#-------------------------------------------------------------------
# Update Date: 11/14/2018
# Create Date: 11/13/2018
# Goal: merge the imputed gen format data of BCAC by chr
# Author: Haoyu Zhang
#-------------------------------------------------------------------
Filesdir <- "/gpfs/gsfs4/users/NC_BW/icogs_onco/genotype/imputed2/onco_imputed"
Files <- dir(Filesdir,pattern="OncoArray_european_merged_b1_15.",full.names=T)
Filesex <- dir(Filesdir,pattern="OncoArray_european_merged_b1_15.chr23",full.names=T)
idx.sex <- Files%in%Filesex
library(gtools)
Files <- Files[!idx.sex]
Files <- mixedsort(Files)
library(dplyr)
merge.to.chr <- rep("c",22)
merge.to.chr <- data.frame(merge.to.chr,stringsAsFactors=F)



for(i in 1:22){
  #take out all the files with chri with grep command
  #fixed option = T only take out the exact match
  idx <- grep(
    paste0(".chr",i,"."),
    Files,
    fixed=T)
  File.chr <- Files[idx]
  #cat multiple files > one file
  merge.code <- paste0("cat ")
  for(j in 1:length(idx)){
    merge.code <- paste0(" ",
                          merge.code,
                          File.chr[j],
                          " ")  
  }
  merge.code <- paste0(merge.code,
                       "> /data/zhangh24/BCAC/impute_onco/chr",
                       i,
                       ".gz")
  merge.to.chr[i,1] <- merge.code
   
}
#write out the command and submit use cluster
write.table(merge.to.chr,file = paste0("/data/zhangh24/breast_cancer_data_analysis/risk_prediction/LD_clumping/code/merge.to.chr.sh"),col.names = F,row.names = F,quote=F)



