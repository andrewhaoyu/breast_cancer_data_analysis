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
setwd("/home/zhangh20/breast_cancer")
#snp_names <- read.csv("iCOGS_Oncoarray_SNPs_names.csv",header = T)
sof <- "/home/zhangh20/breast_cancer/V10/convert.so"
dyn.load(sof)
n <- 109713
snpvalue <- rep(0,n)
subject.file <- "/gpfs/gsfs4/users/NC_BW/icogs_onco/genotype/imputed2/icogs_order.txt.gz"
Icog.order <- read.table(gzfile(subject.file))
pheno.file <- "/home/zhangh20/breast_cancer/V10/pheno_Icog.Rdata"
load(pheno.file)
pheno <- read.csv("/data/zhangh20/breast_cancer/known_SNPs_anlysis/iCOGS_euro_v10_05242017.csv",1)

idx.fil <- Icog.order[,1]%in%pheno$SG_ID
idx.match <- match(pheno$SG_ID,Icog.order[idx.fil,1])
x.covar = pheno[,c(3,5:14)]
country = x.covar[,1]
for(q in 1:10){
  eval(parse(text=paste0("pc",q,"= x.covar[,",q+1,"]")))
}

#Icog.order.match <- Icog.order[idx.fil,1][idx.match]



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
  pvalue_result <- rep(0,num)
  logodds <- rep(0,num)
  sd_odds <- rep(0,num)
  snpid_result <- rep("c",num)
  con <- gzfile(geno.file)
  open(con)
#  for(i in 1:num){
for(i in 1:10){
#if(i%%500==0){
print(i)
#}
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
    save(result,file=paste0("/data/zhangh20/breast_cancer/standard_gwas/case_control/icog/result",i1))
