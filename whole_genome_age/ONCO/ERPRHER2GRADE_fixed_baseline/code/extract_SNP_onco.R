
setwd("/spin1/users/zhangh24/breast_cancer_data_analysis/")

subject.file <- "/gpfs/gsfs4/users/NC_BW/icogs_onco/genotype/imputed2/onco_order.txt.gz"
onco.order <- read.table(gzfile(subject.file))
pheno.file <- "./data/pheno.onco"
load(pheno.file)
n.sub = nrow(pheno)
y.pheno.mis2 <- cbind(pheno$Behaviour1,pheno$PR_status1,pheno$ER_status1,pheno$HER2_status1)
colnames(y.pheno.mis2) = c("Behaviour1","PR_status1",
                           "ER_status1","HER2_status1")


idx.fil <- onco.order[,1]%in%pheno$Onc_ID
n <- length(idx.fil)
snpvalue <- rep(0,n)
idx.match <- match(pheno$Onc_ID,onco.order[idx.fil,1])
#Icog.order.match <- Icog.order[idx.fil,1][idx.match]
library(bc2)
#load("./whole_genome/ONCO/ERPRHER2_fixed/result/score.test.support.onco.ERPRHER2.Rdata")
load("/spin1/users/zhangh24/breast_cancer_data_analysis/whole_genome/ICOG/ERPRHER2GRADE_fixed_baseline/result/extract_list.Rdata")

#extract.list <- extract.list[-c(1721,1722,1724,1726,1727,1728),]
#save(extract.list,file="/spin1/users/zhangh24/breast_cancer_data_analysis/whole_genome/ICOG/ERPRHER2GRADE_fixed_baseline/result/extract_list.Rdata")

snp.onco.extract.id <- extract.list[,(ncol(extract.list)-1),drop=F]
write.table(snp.onco.extract.id,file = paste0("./whole_genome/ONCO/ERPRHER2GRADE_fixed_baseline/result/extract_id_onco.txt"),quote = F,row.names=F)



Filesdir <- "/gpfs/gsfs4/users/NC_BW/icogs_onco/genotype/imputed2/onco_imputed"
Files <- dir(Filesdir,pattern="OncoArray_european_merged_b1_15.",full.names=T)
Filesex <- dir(Filesdir,pattern="OncoArray_european_merged_b1_15.chr23",full.names=T)
idx.sex <- Files%in%Filesex
Files <- Files[!idx.sex]


qctool.command <- rep("c",567)
qctool.command <- data.frame(qctool.command,stringsAsFactors=F)


for(i in 1:567){
  geno.file <- Files[i]
  temp <- paste0("/spin1/users/zhangh24/qctool_v1.4-linux-x86_64/qctool -g ",Files[i]," -incl-rsids /spin1/users/zhangh24/breast_cancer_data_analysis/whole_genome/ONCO/ERPRHER2GRADE_fixed_baseline/result/extract_id_onco.txt -og /spin1/users/zhangh24/breast_cancer_data_analysis/whole_genome/ONCO/ERPRHER2GRADE_fixed_baseline/result/ERPRHER2_extract",i,".txt")
  qctool.command[i,1] <- temp
  
}


write.table(qctool.command,file = paste0("/spin1/users/zhangh24/breast_cancer_data_analysis/whole_genome/ONCO/ERPRHER2GRADE_fixed_baseline/code/qc_extract_snp.sh"),col.names = F,row.names = F,quote=F)















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
