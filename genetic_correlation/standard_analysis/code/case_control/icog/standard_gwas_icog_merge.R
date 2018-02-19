
load("/data/zhangh20/Icog/gwas_info/icog_info.Rdata")


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
###order the SNP files to the right order
idx <- order(Files_sub$chr,Files_sub$p1)
File_sub_order <- Files_sub[order(Files_sub$chr,Files_sub$p1),]



# result.dir <- "/data/zhangh20/breast_cancer/standard_gwas/case_control/icog/"
# result_Files <- dir(result.dir,pattern="result")
# result.idx <- rep(0,length(result_Files))
# for(i in 1:length(result_Files)){
#   result.idx.temp <- as.integer(gsub("result","",result_Files[i]))
#   result.idx[i] <- result.idx.temp
# }

num.total <- 0




setwd("/data/zhangh20/breast_cancer/standard_gwas/case_control/icog/")
for(i in 1:length(Files)){
print(i)

  load(paste0("result",idx[i]))
   temp <- length(result[[1]])
   num.total <- temp+num.total
 }

    id = rep(0,num.total)
   logodds <- rep(0,num.total)
   sd <- rep(0,num.total)



setwd("/data/zhangh20/breast_cancer/standard_gwas/case_control/icog/")
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

 icog_result <- cbind(icog_info,id,logodds,sd)

save(icog_result,file=paste0("/data/zhangh20/breast_cancer/standard_gwas/case_control/icog/icog_result_odds_sd.Rdata"))
