setwd("/data/zhangh24/breast_cancer_data_analysis/")

filedir <- './risk_prediction/intrinsic_subtypes_whole_genome/ONCO/result'
#result_files[1:2000]
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



num <- num.total
#num.total <- 0
rs_id <- rep("c",num.total)
number.of.tumor <- 4
score <- matrix(0,nrow=num,ncol = (number.of.tumor+1))
infor <- matrix(0,nrow = num,ncol = (number.of.tumor+1)^2)
freq.all <- rep(0,num)
#score_baseline <- rep(0,num)
#infor_baseline <- rep(0,num)


load("./risk_prediction/intrinsic_subtypes_whole_genome/ONCO/result/resubmit_id.rdata")
resubmit_resubmit_id <- c(232)
#resubmit_temp <- 0

num.total <- 0
for(i in 1:567){
  print(i)
if(idx[i]%in%resubmit_resubmit_id){
      for (k in 1:1000) {
        #print(k)
        load(
          paste0("./risk_prediction/intrinsic_subtypes_whole_genome/ONCO/result/intrinsic_subytpe_onco_112519_resubmit_resubmit_",idx[i],"_",k)
        )
        temp <- nrow(result[[2]])
        rs_id[num.total+(1:temp)] <- result[[1]]
        score[num.total+(1:temp),] <- result[[2]]
        infor[num.total+(1:temp),] <- result[[3]]
        freq.all[num.total+(1:temp)] <- result[[4]]
        num.total <- temp+num.total
      }
    }else if(idx[i]%in%resubmit_id){
      for (k in 1:37) {
        #print(k)
        load(
          paste0("./risk_prediction/intrinsic_subtypes_whole_genome/ONCO/result/intrinsic_subytpe_onco_112519_resubmit_",idx[i],"_",k)
        )
        temp <- nrow(result[[2]])
        rs_id[num.total+(1:temp)] <- result[[1]]
        score[num.total+(1:temp),] <- result[[2]]
        infor[num.total+(1:temp),] <- result[[3]]
        freq.all[num.total+(1:temp)] <- result[[4]]
        num.total <- temp+num.total
      }
    }
  else{
      for (k in 1:7) {
        #print(k)
        load(
          paste0("./risk_prediction/intrinsic_subtypes_whole_genome/ONCO/result/intrinsic_subytpe_onco_112519_",idx[i],"_",k)
        )
        temp <- nrow(result[[2]])
        rs_id[num.total+(1:temp)] <- result[[1]]
        score[num.total+(1:temp),] <- result[[2]]
        infor[num.total+(1:temp),] <- result[[3]]
        freq.all[num.total+(1:temp)] <- result[[4]] 
        num.total <- temp+num.total
        
      }
    }    
    }
  

load("/data/zhangh24/breast_cancer_data_analysis/whole_genome/ONCO/ERPRHER2GRADE_fixed_baseline/result/onco_info.Rdata")

all.equal(onco_info$rs_id,rs_id)
# idx.diff <- which(onco_info$rs_id!=rs_id)
# length(idx.diff)
# length(which(rs_id=="c"))
CHR <- onco_info[,11]
onco_info <- onco_info[,1:10]


onco_result_casecase <- data.frame(onco_info,score,infor,CHR)

save(onco_result_casecase,file="./risk_prediction/intrinsic_subtypes_whole_genome/ONCO/result/onco_result_intrinsic_subtype.Rdata")
print(1)


idx <- which(onco_result_casecase$rs_id=="chr1_100880328_A_T")
onco_result_casecase[idx,]
