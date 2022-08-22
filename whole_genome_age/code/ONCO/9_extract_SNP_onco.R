
setwd("/data/zhangh24/breast_cancer_data_analysis/")

subject.file <- "/data/NC_BW/icogs_onco/genotype/imputed2/onco_order.txt.gz"
onco.order <- read.table(gzfile(subject.file))

library(data.table)
data2 <- fread("./data/Onco_euro_v10_10232017.csv",header=T)
data2 <- as.data.frame(data2)

#Icog.order.match <- Icog.order[idx.fil,1][idx.match]
#library(bc2)
#load("./whole_genome/ONCO/ERPRHER2_fixed/result/score.test.support.onco.ERPRHER2.Rdata")


######### you need to modify this section#############################
#########load data_clean_select_sig in the investigate results code
load("/data/zhangh24/breast_cancer_data_analysis/whole_genome_age/ICOG/ERPRHER2GRADE_fixed_baseline/result/extract_list.Rdata")
#########snp.onco.extract.id is the column SNP.ONCO
#find the corresponding column for SNP.ICOG in data_clean_select_sig
#your outcome format is different with the old version, you can just do snp.onco.extract.id = data_clean_select_sig$SNP.ONCO
snp.onco.extract.id <- extract.list[,(ncol(extract.list)-1),drop=F]
########change the directory to your local directory
write.table(snp.onco.extract.id,file = paste0("./whole_genome_age/ONCO/ERPRHER2GRADE_fixed_baseline/result/extract_id_onco.txt"),quote = F,row.names=F)
######################################################################

#we don't need CHR 23 during extrac steps since there is no significant SNPs
Filesdir <- "/data/NC_BW/icogs_onco/genotype/imputed2/onco_imputed"
Files <- dir(Filesdir,pattern="OncoArray_european_merged_b1_15.",full.names=T)
Filesex <- dir(Filesdir,pattern="OncoArray_european_merged_b1_15.chr23",full.names=T)
idx.sex <- Files%in%Filesex
Files <- Files[!idx.sex]
library(gtools)
Files <- mixedsort(Files)

#for onco array, the 1-567 files are for auto chromosomes
qctool.command <- rep("c",567)
qctool.command <- data.frame(qctool.command,stringsAsFactors=F)


for(i in 1:567){
  geno.file <- Files[i]
  temp <- paste0("/data/zhangh24/qctool_v1.4-linux-x86_64/qctool -g ",Files[i]," -incl-rsids /data/zhangh24/breast_cancer_data_analysis/whole_genome_age/ONCO/ERPRHER2GRADE_fixed_baseline/result/extract_id_onco.txt -og /data/zhangh24/breast_cancer_data_analysis/whole_genome_age/ONCO/ERPRHER2GRADE_fixed_baseline/result/ERPRHER2_extract",i,".txt")
  qctool.command[i,1] <- temp
  
}

#change the directory to your local directory
write.table(qctool.command,file = paste0("/data/zhangh24/breast_cancer_data_analysis/whole_genome_age/ONCO/ERPRHER2GRADE_fixed_baseline/code/qc_extract_snp.sh"),col.names = F,row.names = F,quote=F)
#write qctool.command to a .sh file to submit to biowulf for running
#you can do swarm -f qc_extract_snp.sh -g 20 --time 4:00:00















# tryCatch(
#   {
#     num <- as.integer(system(paste0("zcat ",geno.file,"| wc -l"),intern=T))
#   },
#   error=function(cond){
#     num <- countLines(geno.file)[1]
#   }
# )
# #num = 22349
# #num <- countLines(geno.file)[1];
# #num <- as.integer(system(paste0("zcat ",geno.file,"| wc -l"),intern=T))
# num.of.tumor <- ncol(y.pheno.mis2)-1
# 
# 
# 
# extract.max <- nrow(extract.list)
# 
# 
# snpid_result <- rep("c",extract.max)
# snpvalue_result <- matrix(0,nrow(y.pheno.mis2),extract.max)
# temp <- 0
# con <- gzfile(geno.file)
# open(con)
# for(i in 1:num){
#   if(i%%500==0){
#     print(i)
#   }
#   oneLine <- readLines(con,n=1)
#   myVector <- strsplit(oneLine," ")
#   snpid <- as.character(myVector[[1]][2])
#   
#   
#   
#   if(snpid%in%extract.list$SNP.ONCO==T){
#     temp <- temp+1
#     snpid_result[temp] <- snpid
#     snpvalue <- rep(0,n)
#     snppro <- as.numeric(unlist(myVector)[6:length(myVector[[1]])])
#     snpvalue <- convert(snppro,n)
#     snpvalue <- snpvalue[idx.fil][idx.match]
#     snpvalue_result[,temp] <- snpvalue
#     
#   }
#   
#   
#   #print(paste0("freq",freq))
#   
#   
# }
# close(con)
# 
# if(temp!=0){
#   snpid_result <- snpid_result[1:temp]
#   snpvalue_result <- snpvalue_result[,1:temp]
# }else{
#   snpid_result <- NULL
#   snpvalue_result <- NULL
# }
# 
# result <- list(snpid_reuslt=snpid_result,snpvalue_result)
# save(result,file=paste0("./whole_genome/ONCO/ERPRHER2_fixed/result/ERPRHER2_fixed_onco_extract",i1,".Rdata"))
