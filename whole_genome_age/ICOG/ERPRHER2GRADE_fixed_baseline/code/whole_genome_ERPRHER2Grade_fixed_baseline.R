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
y.pheno.mis1 <- cbind(pheno$Behaviour1,pheno$PR_status1,pheno$ER_status1,pheno$HER2_status1,pheno$Grade1)
colnames(y.pheno.mis1) = c("Behaviour1","PR_status1",
                           "ER_status1","HER2_status1",
                          "Grade")
idx.fil <- Icog.order[,1]%in%pheno$SG_ID
idx.match <- match(pheno$SG_ID,Icog.order[idx.fil,1])
#Icog.order.match <- Icog.order[idx.fil,1][idx.match]
library(bc2)
load("./whole_genome/ICOG/ERPRHER2GRADE_fixed_baseline/result/score.test.support.icog.ERPRHER2.Rdata")


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



score_result <- matrix(0,num,num.of.tumor+1)
infor_result <- matrix(0,(num.of.tumor+1)*num,num.of.tumor+1)
snpid_result <- rep("c",num)
score_result_baseline <- rep(0,num)
infor_result_baseline <- rep(0,num)

freq.all <- rep(0,num)
con <- gzfile(geno.file)
open(con)
for(i in 1:num){
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
  
  snpvalue <- convert(snppro,n)
  snpvalue <- snpvalue[idx.fil][idx.match]
  freq <- sum(snpvalue)/(2*n.sub)
  freq.all[i] <- freq
  #print(paste0("freq",freq))
  
  tryCatch(
    {
      
      if(freq<0.005|freq>0.995){
        
        score_result[i,] <- 0
        infor_result[((num.of.tumor+1)*i-(num.of.tumor)):((num.of.tumor+1)*i),] <- 0
      }else{
        
        score.test.icog<- ScoreTest(y=y.pheno.mis1,
                                    x=snpvalue,
                                    second.stage.structure="additive",
                                    score.test.support=score.test.support.icog.ERPRHER2Grade,
                                    missingTumorIndicator=888)
        
        score_result[i,]  <- score.test.icog[[1]]
        infor_result[((num.of.tumor+1)*i-(num.of.tumor)):((num.of.tumor+1)*i),] <- score.test.icog[[2]]
        
        score.test.icog.baseline<- ScoreTest(y=y.pheno.mis1,
                                    x=snpvalue,
                                    second.stage.structure="baselineonly",
                                    score.test.support=score.test.support.icog.ERPRHER2Grade,
                                    missingTumorIndicator=888)
        score_result_baseline[i]  <- score.test.icog.baseline[[1]]
        infor_result_baseline[i] <- score.test.icog.baseline[[2]]
      }
      
    },
    error=function(cond) {
      
      score_result[i,] <- 0
      infor_result[((num.of.tumor+1)*i-(num.of.tumor)):((num.of.tumor+1)*i),] <- 0
      score_result_baseline[i]  <- 0
      infor_result_baseline[i] <- 0
      
    })
  
}
close(con)
if(i !=num){
  snpid_result <- snpid_result[1:i]
  score_result <- score_result[1:i,]
  infor_result <- infor_result[1:((num.of.tumor+1)*i),]
  score_result_baseline <- score_result_baseline[1:i]
  infor_result_baseline <- infor_result_baseline[1:i]
  
}
result <- list(snpid_reuslt=snpid_result,score_result=score_result,infor_result=infor_result,freq.all=freq.all,score_result_baseline = score_result_baseline,
               infor_result_baseline=infor_result_baseline)
save(result,file=paste0("./whole_genome/ICOG/ERPRHER2GRADE_fixed_baseline/result/ERPRHER2Grade_fixed_baseline",i1))

