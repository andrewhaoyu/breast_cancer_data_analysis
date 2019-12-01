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

for(i in 1:428){
 
 
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


code <- rep("c",length(unique(resubmit_id1))*155)
temp <- 1
for(i in 1:length(unique(resubmit_id1))){
  for(j in 1:40){
    code[temp] <- paste0("Rscript /gpfs/gsfs11/users/zhangh24/breast_cancer_data_analysis/risk_prediction/intrinsic_subtypes_whole_genome/ONCO/code/intrinsic_subtypes_onco_resubmit_tr.R ",i," ",j)
    temp <- temp+1
  }
}
write.table(code,file ="/gpfs/gsfs11/users/zhangh24/breast_cancer_data_analysis/risk_prediction/intrinsic_subtypes_whole_genome/ONCO/code/intrinsic_subtypes_onco_resubmit_tr.sh",row.names = F,col.names = F,quote = F)



#   if(idx[i]%in%resubimt_resubmimt_id){
#     for (k in 1:70) {
#       #print(k)
#       load(
#         paste0("./whole_genome_age/ONCO/intrinsic_subtypes/result/intrinsic_subytpe_onco_resubmit_resubmit_resubmit",idx[i],"_",k)
#       )
#       temp <- nrow(result[[2]])
#       rs_id[num.total+(1:temp)] <- result[[1]]
#       score[num.total+(1:temp),] <- result[[2]]
#       infor[num.total+(1:temp),] <- result[[3]]
#       freq.all[num.total+(1:temp)] <- result[[4]] 
#       num.total <- temp+num.total
#       if(sum(result[[1]]=="c")!=0){
#         resubmit_temp <- resubmit_temp+1
#         resubmit_id[resubmit_temp,1] <- idx[i]
#         resubmit_id[resubmit_temp,2] <- k
#       }
#     }
#   }else if(idx[i]==327){
#     for (k in 1:1000) {
#       #print(k)
#       load(
#         paste0("./whole_genome_age/ONCO/intrinsic_subtypes/result/intrinsic_subytpe_onco_resubmit_resubmit",idx[i],"_",k)
#       )
#       temp <- nrow(result[[2]])
#       rs_id[num.total+(1:temp)] <- result[[1]]
#       score[num.total+(1:temp),] <- result[[2]]
#       infor[num.total+(1:temp),] <- result[[3]]
#       freq.all[num.total+(1:temp)] <- result[[4]] 
#       num.total <- temp+num.total
#       if(sum(result[[1]]=="c")!=0){
#         resubmit_temp <- resubmit_temp+1
#         resubmit_id[resubmit_temp,1] <- idx[i]
#         resubmit_id[resubmit_temp,2] <- k
#       }
#     }
#   }else if(file_load%in%result_files){
#     for (k in 1:15) {
#       #print(k)
#       load(
#         paste0("./whole_genome_age/ONCO/intrinsic_subtypes/result/intrinsic_subytpe_onco_resubmit",idx[i],"_",k)
#       )
#       temp <- nrow(result[[2]])
#       rs_id[num.total+(1:temp)] <- result[[1]]
#       score[num.total+(1:temp),] <- result[[2]]
#       infor[num.total+(1:temp),] <- result[[3]]
#       freq.all[num.total+(1:temp)] <- result[[4]] 
#       num.total <- temp+num.total
#       if(sum(result[[1]]=="c")!=0){
#         resubmit_temp <- resubmit_temp+1
#         resubmit_id[resubmit_temp,1] <- idx[i]
#         resubmit_id[resubmit_temp,2] <- k
#       }
#     }  
#   }else{
#     for(k in 1:5){
#       load(
#         paste0("./whole_genome_age/ONCO/intrinsic_subtypes/result/intrinsic_subytpe_onco",idx[i],"_",k)
#       )
#       temp <- nrow(result[[2]])
#       rs_id[num.total+(1:temp)] <- result[[1]]
#       score[num.total+(1:temp),] <- result[[2]]
#       infor[num.total+(1:temp),] <- result[[3]]
#       freq.all[num.total+(1:temp)] <- result[[4]]
#       num.total <- temp+num.total
#       if(sum(result[[1]]=="c")!=0){
#         resubmit_temp <- resubmit_temp+1
#         resubmit_id[resubmit_temp,1] <- idx[i]
#         resubmit_id[resubmit_temp,2] <- k
#       }
#     }
#     
#   }
#   
#   
#   
#   
# }

# resubmit_id <- resubmit_id[1:resubmit_temp,]
# unique(resubmit_id[,1])
# length(unique(resubmit_id[,1]))


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
#  load("/data/zhangh24/breast_cancer_data_analysis/whole_genome/ONCO/ERPRHER2GRADE_fixed_baseline/result/onco_result.Rdata")
# 
# 
#  onco_info <- onco_result[,1:10]
#  CHR <- onco_result[,41]
#  onco_info <- cbind(onco_info,CHR)
#  save(onco_info,file="/data/zhangh24/breast_cancer_data_analysis/whole_genome/ONCO/ERPRHER2GRADE_fixed_baseline/result/onco_info.Rdata")




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


idx <- which(onco_result_casecase$rs_id=="rs1416885:100010765:C:G")
onco_result_casecase[idx,]
