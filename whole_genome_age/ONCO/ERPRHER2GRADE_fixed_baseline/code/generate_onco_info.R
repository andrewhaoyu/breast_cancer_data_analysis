#goal: generate the information matrix for onco
#generate the information matrix for 22 autochromosomes
onco_info <- data.frame(snp_id = rep("c",num.total),rs_id = rep("c",num.total),
                        position=rep(0,num.total),exp_freq_a1=rep(0,num.total),info=rep(0,num.total),
                        certainty=rep(0,num.total),type=rep(0,num.total),info_type0=rep(0,num.total),
                        concord_type0=rep(0,num.total),r2_type0=rep(0,num.total),stringsAsFactors=F)
CHR <- rep(0,num.total)
num.total <-  0
for(i in 1:22){
  print(i)
  filedir <- paste0("/gpfs/gsfs4/users/NC_BW/icogs_onco/genotype/imputed2/onco_info_files/chr",i)
  files <- dir(filedir,pattern="txt_info",full.names=T)
  files_num <- gsub(paste0(filedir,"/OncoArray_chr",i,"_euro15_phased_"),
                    "",files)
  files_num <- gsub(".txt_info","",files_num)
  files_num <- strsplit(files_num,"_")
  files_num <- as.integer(unlist(files_num)[seq(1,2*length(files_num),2)])
  idx <- order(files_num)
  for(j in 1:length(idx)){
    #print(j)
    data <- read.table(files[idx[j]],header=T,stringsAsFactors=F)
    temp <- nrow(data)
    onco_info[num.total+(1:temp),1:3] <- data[,1:3]
    onco_info[num.total+(1:temp),4:10] <- data[,6:12]
    CHR[num.total+(1:temp)] <- i
    num.total <- temp+num.total
  }
  
}

onco_result <- data.frame(onco_info,score,infor,CHR)

save(onco_result,file="/spin1/users/zhangh24/breast_cancer_data_analysis/whole_genome/ONCO/ERPRHER2_fixed/result/onco_result.Rdata")
print(1)



##generate information matrix for sex chromosome
num.total <- 730733 
onco_info <- data.frame(snp_id = rep("c",num.total),rs_id = rep("c",num.total),
                        position=rep(0,num.total),
                        A0 = rep("c",num.total),
                        A1 = rep("c",num.total),
                        exp_freq_a1=rep(0,num.total),info=rep(0,num.total),
                        certainty=rep(0,num.total),type=rep(0,num.total),info_type0=rep(0,num.total),
                        concord_type0=rep(0,num.total),r2_type0=rep(0,num.total),stringsAsFactors=F)
CHR <- rep(0,num.total)
num.total <-  0
for(i in 23){
  print(i)
  filedir <- paste0("/gpfs/gsfs4/users/NC_BW/icogs_onco/genotype/imputed2/onco_info_files/chr",i)
  files <- dir(filedir,pattern="txt_info",full.names=T)
  files_num <- gsub(paste0(filedir,"/OncoArray_chr",i,"_euro15_phased_"),
                    "",files)
  files_num <- gsub(".txt_info","",files_num)
  files_num <- strsplit(files_num,"_")
  files_num <- as.integer(unlist(files_num)[seq(1,2*length(files_num),2)])
  idx <- order(files_num)
  for(j in 1:length(idx)){
    #print(j)
    data <- read.table(files[idx[j]],header=T,stringsAsFactors=F)
    temp <- nrow(data)
    onco_info[num.total+(1:temp),] <- data
    CHR[num.total+(1:temp)] <- i
    num.total <- temp+num.total
  }
  
}



save(onco_info,file="/spin1/users/zhangh24/breast_cancer_data_analysis/whole_genome_age/ONCO/ERPRHER2GRADE_fixed_baseline/result/onco_infor_sex.Rdata")
print(1)

