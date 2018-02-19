
load("/data/zhangh20/onco/gwas_info/onco_info.Rdata")

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


# result.dir <- "/data/zhangh20/breast_cancer/standard_gwas/case_control/onco/"
# result_Files <- dir(result.dir,pattern="result")
# result.idx <- rep(0,length(result_Files))
# for(i in 1:length(result_Files)){
#   result.idx.temp <- as.integer(gsub("result","",result_Files[i]))
#   result.idx[i] <- result.idx.temp
# }

num.total <- 0


setwd("/data/zhangh20/breast_cancer/standard_gwas/ER+_control/onco/")
for(i in 1:length(Files)){
print(i)

  load(paste0("result",idx[i]))
   temp <- length(result[[1]])
   num.total <- temp+num.total
 }

    id = rep(0,num.total)
   logodds <- rep(0,num.total)
   sd <- rep(0,num.total)

   setwd("/data/zhangh20/breast_cancer/standard_gwas/ER+_control/onco/")
   temp <- 0
   num.total <- 0
   for(i in 1:length(Files)){
     print(i)
     load(paste0("result",idx[i]))
     temp <- length(result[[1]])
     id[(num.total+1):(num.total+temp)] = result[[1]]
     logodds[(num.total+1):(num.total+temp)] <- result[[2]]
     sd[(num.total+1):(num.total+temp)] <- result[[3]]
     num.total = temp+num.total
   }

    onco_result <- cbind(onco_info,id,logodds,sd)

   save(onco_result,file=paste0("/data/zhangh20/breast_cancer/standard_gwas/ER+_control/onco/onco_result_odds_sd.Rdata"))
