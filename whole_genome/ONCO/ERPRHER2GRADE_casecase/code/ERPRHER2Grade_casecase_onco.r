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
i1 <- as.numeric(arg[[1]])
i2 <- as.numeric(arg[[2]])
print(i1)
print(i2)

library(R.utils)
setwd("/spin1/users/zhangh24/breast_cancer_data_analysis/")

subject.file <- "/gpfs/gsfs4/users/NC_BW/icogs_onco/genotype/imputed2/onco_order.txt.gz"
onco.order <- read.table(gzfile(subject.file))
pheno.file <- "./data/pheno.onco"
load(pheno.file)
n.sub = nrow(pheno)
y.pheno.mis2 <- cbind(pheno$Behaviour1,pheno$PR_status1,pheno$ER_status1,pheno$HER2_status1,pheno$Grade)
colnames(y.pheno.mis2) = c("Behaviour1","PR_status1",
                           "ER_status1","HER2_status1",
                           "Grade")


idx.fil <- onco.order[,1]%in%pheno$Onc_ID
n <- length(idx.fil)
snpvalue <- rep(0,n)
idx.match <- match(pheno$Onc_ID,onco.order[idx.fil,1])
#Icog.order.match <- Icog.order[idx.fil,1][idx.match]
library(bc2)
#load("./whole_genome/ONCO/ERPRHER2GRADE_fixed_baseline/result/score.test.support.onco.ERPRHER2Grade.Rdata")
load("./whole_genome/ONCO/ERPRHER2GRADE_fixed_baseline/result/delta0.onco.Rdata")
load("./whole_genome/ONCO/ERPRHER2GRADE_fixed_baseline/result/z.standard.Rdata")
x.all.covar <- pheno[,2:11]






Filesdir <- "/gpfs/gsfs4/users/NC_BW/icogs_onco/genotype/imputed2/onco_imputed"
Files <- dir(Filesdir,pattern="OncoArray_european_merged_b1_15.",full.names=T)
Filesex <- dir(Filesdir,pattern="OncoArray_european_merged_b1_15.chr23",full.names=T)
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
rm(pheno)
num.of.tumor <- ncol(y.pheno.mis2)-1
size = 5
start.end <- startend(num,size,i2)
start <- start.end[1]
end <- start.end[2]
file.num <- end-start+1


score_result <- matrix(0.1,file.num,num.of.tumor)
infor_result <- matrix(0.1,(num.of.tumor)*file.num,num.of.tumor)
snpid_result <- rep("c",file.num)

freq.all <- rep(0,file.num)
temp <- 0



con <- gzfile(geno.file)
open(con)
for(i in 1:num){
  if(i%%500==0){
    print(i)
  }

  oneLine <- readLines(con,n=1)
  if(i>=start){
    if(temp%%100==0){
      print(paste0("temp",temp))
    }
    temp = temp+1
    if(temp%%100==0){
      print(paste0("temp",temp))
    }
    myVector <- strsplit(oneLine," ")
    snpid <- as.character(myVector[[1]][2])
    snpid_result[i] <- snpid
    snpvalue <- rep(0,n)
    
    
    snppro <- as.numeric(unlist(myVector)[6:length(myVector[[1]])])
    if(length(snppro)!=(3*n)){
      break
    }
    
    snpvalue <- convert(snppro,n)
    snpvalue <- snpvalue[idx.fil][idx.match]
    freq <- sum(snpvalue)/(2*n.sub)
    freq.all[i] <- freq
    #print(paste0("freq",freq))
    
    tryCatch(
      {
        
        if(freq<0.005|freq>0.995){
          
          score_result[temp,] <- 0
          infor_result[((num.of.tumor)*temp-(num.of.tumor-1)):((num.of.tumor)*temp),] <- 0
        }else{
          
          score.test.support.onco.casecase <- ScoreTestSupportMixedModel(y=y.pheno.mis2,
                                                                         baselineonly = snpvalue,
                                                                         additive=x.all.covar,
                                                                         missingTumorIndicator = 888,
                                                                         delta0=delta0.onco)
          score.test.onco.casecase<- ScoreTestMixedModel(y=y.pheno.mis2,
                                                         x=snpvalue,
                                                         z.design = z.standard, 
                                                         score.test.support= score.test.support.onco.casecase,
                                                         missingTumorIndicator=888)
          
          score_result[temp,]  <- score.test.onco.casecase[[1]]
          infor_result[((num.of.tumor)*temp-(num.of.tumor-1)):((num.of.tumor)*temp),] <- score.test.onco.casecase[[2]]
          rm(score.test.support.onco.casecase)
          rm(score.test.onco.casecase)
          gc()
        }
        
      },
      error=function(cond) {
        
        score_result[temp,] <- 0
        infor_result[((num.of.tumor)*temp-(num.of.tumor-1)):((num.of.tumor)*temp),] <- 0
        
      })
    
  }
  if(i==end){
    break
  }
 
}
close(con)

result <- list(snpid_reuslt=snpid_result,score_result=score_result,infor_result=infor_result,freq.all=freq.all)
save(result,file=paste0("./whole_genome/ONCO/ERPRHER2GRADE_casecase/result/ERPRHER2Grade_casecase_onco",i1,"_",i2))
