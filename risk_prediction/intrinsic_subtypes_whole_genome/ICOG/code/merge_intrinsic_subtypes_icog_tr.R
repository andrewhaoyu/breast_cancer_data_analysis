# setwd("/spin1/users/zhangh24/breast_cancer_data_analysis/")
# filedir <- './whole_genome_age/ICOG/Intrinsic_subtypes/result/'
# files <- dir(filedir,pattern="intrinsic_subytpe_icog")
# total <- 564*5
# missingid <- matrix(0,total,2)
# temp <- 0
# for(i1 in 1:564){
#   print(i1)
#   for(i2 in 1:5){
#     text <- paste0("intrinsic_subytpe_icog",i1,"_",i2)
#     if((text%in%files)==F){
#       temp <- temp+1
#       missingid[temp,] <- c(i1,i2)
#     }
#   }
# }
# missingid <- missingid[1:temp,]
# icog.unique.resubmit <- unique(missingid[,1])
# save(icog.unique.resubmit,file="./whole_genome_age/ICOG/Intrinsic_subtypes/result/icog.unique.resubmit.Rdata")
# submit <- rep("c",length(icog.unique.resubmit)*15)
# temp <- 1
# for(i in 1:length(icog.unique.resubmit)){
#   for(j in 1:15){
#     submit[temp] <- paste0("Rscript /spin1/users/zhangh24/breast_cancer_data_analysis/whole_genome_age/ICOG/Intrinsic_subtypes/code/intrinsic_subtype_icog.R ",icog.unique.resubmit[i]," ",j)
#     temp <- temp+1
#   }
#   
# }
# write.table(submit,file="/spin1/users/zhangh24/breast_cancer_data_analysis/whole_genome_age/ICOG/Intrinsic_subtypes/code/icog_resubmit.sh",
#             row.names=F,quote=F,col.names=F)







setwd("/spin1/users/zhangh24/breast_cancer_data_analysis/")
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



load("/spin1/users/zhangh24/breast_cancer_data_analysis/whole_genome_age/ICOG/ERPRHER2GRADE_fixed_baseline/result/Icog_result.Rdata")
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
  # if(idx[i]%in%resubimt_resubmimt_id){
  #   for(k in 1:2){
  #     load(paste0("./whole_genome_age/ICOG/Intrinsic_subtypes/result/intrinsic_subytpe_icog_resubmit_resubmit_resubmit",idx[i],"_",k))
  #     temp <- nrow(result[[2]])
  #     rs_id[num.total+(1:temp)] <- result[[1]]
  #     score[num.total+(1:temp),] <- result[[2]]
  #     infor[num.total+(1:temp),] <- result[[3]]
  #     num.total <- temp+num.total
  #     if(sum(result[[1]]=="c")!=0){
  #       resubmit_temp <- resubmit_temp+1
  #       resubmit_id[resubmit_temp,1] <- idx[i]
  #       resubmit_id[resubmit_temp,2] <- k
  #     }
  #   }
  # }else if(idx[i]==413){
  #   for(k in 1:1000){
  #     load(paste0("./whole_genome_age/ICOG/Intrinsic_subtypes/result/intrinsic_subytpe_icog_resubmit_resubmit",idx[i],"_",k))
  #     temp <- nrow(result[[2]])
  #     rs_id[num.total+(1:temp)] <- result[[1]]
  #     score[num.total+(1:temp),] <- result[[2]]
  #     infor[num.total+(1:temp),] <- result[[3]]
  #     num.total <- temp+num.total
  #     if(sum(result[[1]]=="c")!=0){
  #       resubmit_temp <- resubmit_temp+1
  #       resubmit_id[resubmit_temp,1] <- idx[i]
  #       resubmit_id[resubmit_temp,2] <- k
  #     }
  #   }
  # }else if(file_load%in%result_files){
  #   for(k in 1:15){
  #     load(paste0("./whole_genome_age/ICOG/Intrinsic_subtypes/result/intrinsic_subytpe_icog_resubmit",idx[i],"_",k))
  #     temp <- nrow(result[[2]])
  #     rs_id[num.total+(1:temp)] <- result[[1]]
  #     score[num.total+(1:temp),] <- result[[2]]
  #     infor[num.total+(1:temp),] <- result[[3]]
  #     num.total <- temp+num.total
  #     if(sum(result[[1]]=="c")!=0){
  #       resubmit_temp <- resubmit_temp+1
  #       resubmit_id[resubmit_temp,1] <- idx[i]
  #       resubmit_id[resubmit_temp,2] <- k
  #     }
  #   }
  # }else{
  #   for(k in 1:5){
  #     load(paste0("./whole_genome_age/ICOG/Intrinsic_subtypes/result/intrinsic_subytpe_icog",idx[i],"_",k))
  #     temp <- nrow(result[[2]])
  #     rs_id[num.total+(1:temp)] <- result[[1]]
  #     score[num.total+(1:temp),] <- result[[2]]
  #     infor[num.total+(1:temp),] <- result[[3]]
  #     # for(j in 1:temp){
  #     #   infor_j <- result[[3]][(number.of.tumor*j-(number.of.tumor-1)):((number.of.tumor)*j),]
  #     #   infor[num.total+j,] <- as.vector(infor_j)
  #     # }
  #     # if(num.total< 12327300&(num.total+temp)> 12327300){
  #     #   print(c(i,k))
  #     # }
  #     num.total <- temp+num.total
  #     if(sum(result[[1]]=="c")!=0){
  #       resubmit_temp <- resubmit_temp+1
  #       resubmit_id[resubmit_temp,1] <- idx[i]
  #       resubmit_id[resubmit_temp,2] <- k
  #     }
  #   }
  #   
    
  
  
  
  
  
#}
# resubmit_id <- resubmit_id[1:resubmit_temp,]
# unique(resubmit_id[,1])



# k <- 1
# load(paste0("./whole_genome_age/ICOG/Intrinsic_subtypes/result/intrinsic_subytpe_icog",idx[i],"_",k))
# idx.try <- which(result[[1]]=="c")
# print(length(idx.try))
#try <- merge(icog_info,rs_id,by.x=rs_id,by.y=rs_id,all=T)



#####to get the total number of SNPs from the information files


# icog_info <- data.frame(snp_id = rep("c",num.total),rs_id = rep("c",num.total),
#                         position=rep(0,num.total),exp_freq_a1=rep(0,num.total),info=rep(0,num.total),
#                         certainty=rep(0,num.total),type=rep(0,num.total),info_type0=rep(0,num.total),
#                         concord_type0=rep(0,num.total),r2_type0=rep(0,num.total),stringsAsFactors=F)
# CHR <- rep(0,num.total)
# num.total <-  0
# library(data.table)
# for(i in 1:22){
#   print(i)
#   filedir <- paste0("/gpfs/gsfs4/users/NC_BW/icogs_onco/genotype/imputed2/icogs_info_files/chr",i)
#   files <- dir(filedir,pattern="txt_info",full.names=T)
#   files_num <- gsub(paste0(filedir,"/icogs_euro12_chr",i,"_phased"),
#                     "",files)
#   files_num <- gsub(".txt_info","",files_num)
#   files_num <- strsplit(files_num,"_")
#   files_num <- as.integer(unlist(files_num)[seq(1,2*length(files_num),2)])
#   idx <- order(files_num)
#   for(j in 1:length(idx)){
#     #print(j)
#     data <- as.data.frame(fread(files[idx[j]],header=T,stringsAsFactors=F))
#     
#     temp <- nrow(data)
#     icog_info[num.total+(1:temp),] <- data
#     CHR[num.total+(1:temp)] <- i
#     num.total <- temp+num.total
#   }
#   
# }
# 
load("/spin1/users/zhangh24/breast_cancer_data_analysis/whole_genome/ICOG/ERPRHER2GRADE_fixed_baseline/result/icog_info.Rdata")
# icog_info <- cbind(icog_info,CHR)
# save(icog_info,file="/spin1/users/zhangh24/breast_cancer_data_analysis/whole_genome/ICOG/ERPRHER2GRADE_fixed_baseline/result/icog_info.Rdata")
all.equal(icog_info$rs_id,rs_id)
#idx.diff <- which(icog_info$rs_id!=rs_id)
CHR <- icog_info[,11]
icog_info <- icog_info[,1:10]

icog_result_casecase <- data.frame(icog_info,score,infor,CHR)






save(icog_result_casecase,file="/spin1/users/zhangh24/breast_cancer_data_analysis/risk_prediction/intrinsic_subtypes_whole_genome/ICOG/result/Icog_result_intrinsic_subtype.Rdata")
# icog_result_baseline <- data.frame(icog_info,score_baseline,infor_baseline,CHR)
# save(icog_result_baseline,file="/spin1/users/zhangh24/breast_cancer_data_analysis/whole_genome/ICOG/ERPRHER2GRADE_fixed_baseline/result/Icog_result_baseline.Rdata")
# print(1)



