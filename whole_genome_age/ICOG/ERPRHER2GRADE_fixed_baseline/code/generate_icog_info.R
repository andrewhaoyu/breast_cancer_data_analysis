#####to get the total number of SNPs from the information files

num.total <- 20667573
icog_info <- data.frame(snp_id = rep("c",num.total),rs_id = rep("c",num.total),
                        position=rep(0,num.total),exp_freq_a1=rep(0,num.total),info=rep(0,num.total),
                        certainty=rep(0,num.total),type=rep(0,num.total),info_type0=rep(0,num.total),
                        concord_type0=rep(0,num.total),r2_type0=rep(0,num.total),stringsAsFactors=F)
CHR <- rep(0,num.total)
num.total <-  0

for(i in 1:22){
  print(i)
  filedir <- paste0("/gpfs/gsfs4/users/NC_BW/icogs_onco/genotype/imputed2/icogs_info_files/chr",i)
  files <- dir(filedir,pattern="txt_info",full.names=T)
  files_num <- gsub(paste0(filedir,"/icogs_euro12_chr",i,"_phased"),
                    "",files)
  files_num <- gsub(".txt_info","",files_num)
  files_num <- strsplit(files_num,"_")
  files_num <- as.integer(unlist(files_num)[seq(1,2*length(files_num),2)])
  idx <- order(files_num)
  for(j in 1:length(idx)){
    print(j)
    #print(j)
    data <- as.data.frame(fread(files[idx[j]],header=T,stringsAsFactors=F))
    temp <- nrow(data)
    icog_info[num.total+(1:temp),] <- data
    CHR[num.total+(1:temp)] <- i
    num.total <- temp+num.total
  }

}

######
try <- strsplit(data$rs_id,":")
A1 <- rep("c",nrow(data))
A2 <- rep("c",nrow(data))
temp <- strsplit(data$rs_id,":")
for(i in 1:nrow(data)){
  print(i)
  
  A1[i] <- temp[[i]][3]
  A2[i] <- temp[[i]][4]
}


icog_info <- icog_result[,1:10]

save(icog_info,file="/spin1/users/zhangh24/breast_cancer_data_analysis/whole_genome/ICOG/ERPRHER2GRADE_fixed_baseline/result/icog_info.Rdata")






#####to get the total number of SNPs from the information files for sex chromosome

num.total <- 728693
icog_info <- data.frame(snp_id = rep("c",num.total),rs_id = rep("c",num.total),
                        position=rep(0,num.total),
                        a0 = rep("c",num.total),
                        a1 = rep("c",num.total),
                        exp_freq_a1=rep(0,num.total),
                        info=rep(0,num.total),
                        certainty=rep(0,num.total),type=rep(0,num.total),info_type0=rep(0,num.total),
                        concord_type0=rep(0,num.total),r2_type0=rep(0,num.total),stringsAsFactors=F)
CHR <- rep(0,num.total)
num.total <-  0
library(data.table)
for(i in 23){
  print(i)
  filedir <- paste0("/gpfs/gsfs4/users/NC_BW/icogs_onco/genotype/imputed2/icogs_info_files/chr",i)
  files <- dir(filedir,pattern="txt_info",full.names=T)
  files_num <- gsub(paste0(filedir,"/icogs_chr",i,"_euro12_phased_"),
                    "",files)
  files_num <- gsub("_new_phase3.txt_info","",files_num)
  files_num <- strsplit(files_num,"_")
  files_num <- as.integer(unlist(files_num)[seq(1,2*length(files_num),2)])
  idx <- order(files_num)
  for(j in 1:length(idx)){
    print(j)
    #print(j)
    data <- as.data.frame(fread(files[idx[j]],header=T,stringsAsFactors=F))
    temp <- nrow(data)
    icog_info[num.total+(1:temp),] <- data
    CHR[num.total+(1:temp)] <- i
    num.total <- temp+num.total
  }
  
}



save(icog_info,file="/spin1/users/zhangh24/breast_cancer_data_analysis/whole_genome/ICOG/ERPRHER2GRADE_fixed_baseline/result/icog_info_sex.Rdata")
