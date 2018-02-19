rm(list=ls())
arg <- commandArgs(trailingOnly=T)
i <- as.numeric(arg[[1]])
i1 <- i
print(i)
#args=(commandArgs(TRUE))
#for(p in 1:length(args)){
 #       eval(parse(text=args[[p]]))
  #  }
#print(i)
#i1 <- i
library(R.utils)
setwd("/home/zhangh20/breast_cancer")

sof <- "/home/zhangh20/breast_cancer/V10/convert.so"
dyn.load(sof)

subject.file <- "/gpfs/gsfs4/users/NC_BW/icogs_onco/genotype/imputed2/onco_order.txt.gz"
onco.order <- read.table(gzfile(subject.file))

pheno  <- read.csv("/data/zhangh20/breast_cancer/known_SNPs_anlysis/Onco_euro_v10_05242017.csv",1)
idx.ER = which(pheno$ER_status1==0|pheno$Behaviour1==0)
pheno = pheno[idx.ER,]
idx.fil <- onco.order[,1]%in%pheno$Onc_ID
n <- length(idx.fil)
snpvalue <- rep(0,n)
idx.match <- match(pheno$Onc_ID,onco.order[idx.fil,1])
#Icog.order.match <- Icog.order[idx.fil,1][idx.match]
x.covar = pheno[,c(4,5:14)]
country = x.covar[,1]
for(q in 1:10){
  eval(parse(text=paste0("pc",q,"= x.covar[,",q+1,"]")))
}


Filesdir <- "/gpfs/gsfs4/users/NC_BW/icogs_onco/genotype/imputed2/onco_imputed"
Files <- dir(Filesdir,pattern="OncoArray_european_merged_b1_15.",full.names=T)
Filesex <- dir(Filesdir,pattern="OncoArray_european_merged_b1_15.chr23",full.names=T)
idx.sex <- Files%in%Filesex
Files <- Files[!idx.sex]
geno.file <- Files[i1]
#if(geno.file%in%Filesex){
 # result <- NULL
 # save(result,file=paste0("/gpfs/gsfs4/users/NC_BW/icogs_onco/genotype/imputed2/onco_gwas/result_sex",i1))
#}else{
#geno.file <- Files[i1]
#num <- countLines(geno.file)[1];

#tryCatch(
#{
num <- as.integer(system(paste0("zcat ",geno.file,"| wc -l"),intern=T))
#},
#error=function(cond){
#num <- countLines(geno.file)[1]
#}
#)

#num <- as.integer(system(paste0("zcat ",geno.file,"| wc -l"),intern=T))


pvalue_result <- rep(0,num)
logodds <- rep(0,num)
sd_odds <- rep(0,num)
snpid_result <- rep("c",num)
con <- gzfile(geno.file)
open(con)

for(i in 1:num){
if(i%%500==0){
print(i)
}
  oneLine <- readLines(con,n=1)
  myVector <- strsplit(oneLine," ")
  snpid <- as.character(myVector[[1]][2])
  snpid_result[i] <- snpid
  snpvalue <- rep(0,n)

  snppro <- as.numeric(unlist(myVector)[6:length(myVector[[1]])])
  if(length(snppro)!=(3*n)){
    break
  }
  temp <- .C("convert",as.integer(n),as.numeric(snppro),as.numeric(snpvalue))
  snpvalue <- temp[[3]]
  snpvalue <- snpvalue[idx.fil][idx.match]

  tryCatch(
        {

                      if(sum(snpvalue)==0){
            pvalue_result[i] <- 1
            logodds[i] <- 0
            sd_odds[i] <- 0
     }else{
       model = glm(pheno$Behaviour1~snpvalue+country+pc1+pc2+pc3+pc4+pc5+pc6+pc7+pc8+pc9+pc10
         ,family=binomial())
         model.sum = summary(model)
         logodds[i] = model.sum$coefficients[2,1]
         sd_odds[i] = model.sum$coefficients[2,2]
         pvalue_result[i] = model.sum$coefficients[2,4]
     }

   },
   error=function(cond) {
     pvalue_result[i] <- 1
     logodds[i] <- 0
     sd_odds[i] <- 0
   })

}
close(con)
if(i !=num){
snpid_result <- snpid_result[1:i]
pvalue_result <- pvalue_result[1:i]
logodds = logodds[1:i]
sd_odds = sd_odds[1:i]
}
result <- list(snpid_reuslt=snpid_result,logodds=logodds,sd_odds=sd_odds,
    pvalue_result=pvalue_result)
  save(result,file=paste0("/data/zhangh20/breast_cancer/standard_gwas/ER-_control/onco/result",i1))
