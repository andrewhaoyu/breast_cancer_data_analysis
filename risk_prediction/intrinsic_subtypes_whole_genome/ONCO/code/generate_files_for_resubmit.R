#goal: load the results files and generate the resubmission jobs id
setwd("/data/zhangh24/breast_cancer_data_analysis/")

filedir <- './risk_prediction/intrinsic_subtypes_whole_genome/ONCO/result'

files <- dir(filedir,pattern="intrinsic_subytpe_onco_112519")

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
  #print(i)
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
load("/data/zhangh24/breast_cancer_data_analysis/whole_genome_age/ONCO/ERPRHER2GRADE_fixed_baseline/result/onco_result.Rdata")
num.total <- nrow(onco_result)

resubmit_id1 <- NULL
resubmit_id2 <- NULL

#for(i in 1:428){
  for(i in 1:567){ 
 
    for (k in 1:7) {
      #print(k)
      file_temp = 
        paste0("intrinsic_subytpe_onco_112519_",idx[i],"_",k)
      if(file_temp%in%files==F){
        resubmit_id1 <- c(resubmit_id1,idx[i])
        resubmit_id2 <- c(resubmit_id2,k)
      }
      
    }
  }    


code <- rep("c",length(unique(resubmit_id1))*37)
temp <- 1
for(i in 1:length(unique(resubmit_id1))){
  for(j in 1:37){
    code[temp] <- paste0("Rscript /gpfs/gsfs11/users/zhangh24/breast_cancer_data_analysis/risk_prediction/intrinsic_subtypes_whole_genome/ONCO/code/intrinsic_subtypes_onco_resubmit_tr.R ",unique(resubmit_id1)[i]," ",j)
    temp <- temp+1
  }
}
write.table(code,file ="/gpfs/gsfs11/users/zhangh24/breast_cancer_data_analysis/risk_prediction/intrinsic_subtypes_whole_genome/ONCO/code/intrinsic_subtypes_onco_resubmit_tr.sh",row.names = F,col.names = F,quote = F)
resubmit_id <- unique(resubmit_id1)
save(resubmit_id,file = "./risk_prediction/intrinsic_subtypes_whole_genome/ONCO/result/resubmit_id.rdata")
