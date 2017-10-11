rm(list=ls())
#args=(commandArgs(TRUE))
#for(p in 1:length(args)){
#       eval(parse(text=args[[p]]))
#  }
#print(i)
#i1 <- i
#commandarg <- commandArgs(trailingOnly=F)
#myarg <- commandarg[length(commandarg)]
#myarg <- sub("-","",myarg)
#i <- as.numeric(myarg)
#print(i)
#pheno is ICOGS,data2 is onco_array
arg <- commandArgs(trailingOnly=T)
i <- as.numeric(arg[[1]])
i1 <- i
print(i)
library(R.utils)
setwd("/spin1/users/zhangh24/breast_cancer_data_analysis/")

n <- 109713
snpvalue <- rep(0,n)
subject.file <- "/gpfs/gsfs4/users/NC_BW/icogs_onco/genotype/imputed2/icogs_order.txt.gz"
Icog.order <- read.table(gzfile(subject.file))
pheno.file <- "./data/pheno.Icog"
load(pheno.file)
n.sub = nrow(pheno)
y.pheno.mis1 <- cbind(pheno$Behaviour1,pheno$PR_status1,pheno$ER_status1,pheno$HER2_status1)
colnames(y.pheno.mis1) = c("Behaviour1","PR_status1",
                           "ER_status1","HER2_status1")
idx.fil <- Icog.order[,1]%in%pheno$SG_ID
idx.match <- match(pheno$SG_ID,Icog.order[idx.fil,1])
#Icog.order.match <- Icog.order[idx.fil,1][idx.match]
library(bc2)
#load("./whole_genome/ICOG/ERPRHER2_fixed/result/score.test.support.icog.ERPRHER2.Rdata")
load("/spin1/users/zhangh24/breast_cancer_data_analysis/whole_genome/ICOG/ERPRHER2_fixed/result/extract_list.Rdata")


Filesdir <- "/gpfs/gsfs4/users/NC_BW/icogs_onco/genotype/imputed2/icogs_imputed/"
Files <- dir(Filesdir,pattern="icogs_merged_b1_12.",full.names=T)
Filesex <- dir(Filesdir,pattern="icogs_merged_b1_12.chr23",full.names=T)
idx.sex <- Files%in%Filesex
Files <- Files[!idx.sex]
geno.file <- Files[i1]
tryCatch(
  {
    num <- as.integer(system(paste0("zcat ",geno.file,"| wc -l"),intern=T))
  },
  error=function(cond){
    num <- countLines(geno.file)[1]
  }
)
#num = 22349
#num <- countLines(geno.file)[1];
#num <- as.integer(system(paste0("zcat ",geno.file,"| wc -l"),intern=T))
num.of.tumor <- ncol(y.pheno.mis1)-1

extract.max <- nrow(extract.list)


snpid_result <- rep("c",extract.max)
snpvalue_result <- matrix(0,nrow(y.pheno.mis1),extract.max)
temp <- 0

con <- gzfile(geno.file)
open(con)
for(i in 1:num){
  if(i%%500==0){
    print(i)
  }
  oneLine <- readLines(con,n=1)
  myVector <- strsplit(oneLine," ")
  snpid <- as.character(myVector[[1]][2])
  
  if(snpid%in%extract.list$SNP.ICOGS){
    temp <- temp+1
    snpid_result[temp] <- snpid
    snpvalue <- rep(0,n)
    
    
    snppro <- as.numeric(unlist(myVector)[6:length(myVector[[1]])])
    
    snpvalue <- convert(snppro,n)
    snpvalue <- snpvalue[idx.fil][idx.match]
    snpvalue_result[,temp] <- snpvalue
    
  }
}
  close(con)
  if(temp!=0){
    snpid_result <- snpid_result[1:temp]
    snpvalue_result <- snpvalue_result[,1:temp]
  }else{
    snpid_result <- NULL
    snpvalue_result <- NULL
  }
result <- list(snpid_reuslt=snpid_result,snpvalue_result=snpvalue_result)
save(result,file=paste0("./whole_genome/ICOG/ERPRHER2_fixed/result/ERPRHER2_fixed_extracted.Rdata",i1))
