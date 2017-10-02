

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
# arg <- commandArgs(trailingOnly=T)
# i1 <- as.numeric(arg[[1]])
# i2 <- as.numeric(arg[[2]])
# print(i1)
# print(i2)
i1= 1

library(R.utils)
setwd("/spin1/users/zhangh24/breast_cancer_data_analysis/")
load("./whole_genome/ICOG/ERPRHER2GRADE_casecase/result/ERPRHER2Grade_casecase_test")

n <- 109713
snpvalue <- rep(0,n)
subject.file <- "/gpfs/gsfs4/users/NC_BW/icogs_onco/genotype/imputed2/icogs_order.txt.gz"
Icog.order <- read.table(gzfile(subject.file))
pheno.file <- "./data/pheno.Icog"
load(pheno.file)
n.sub = nrow(pheno)

y.pheno.mis1 <- cbind(pheno$Behaviour1,pheno$PR_status1,pheno$ER_status1,pheno$HER2_status1,pheno$Grade1)
x.all.covar <- pheno[,2:11]
colnames(y.pheno.mis1) = c("Behaviour1","PR_status1",
                           "ER_status1","HER2_status1",
                           "Grade")
idx.fil <- Icog.order[,1]%in%pheno$SG_ID
idx.match <- match(pheno$SG_ID,Icog.order[idx.fil,1])
#Icog.order.match <- Icog.order[idx.fil,1][idx.match]
library(bc2)
load("./whole_genome/ICOG/ERPRHER2GRADE_fixed_baseline/result/delta0.icog.Rdata")
load("./whole_genome/ICOG/ERPRHER2GRADE_fixed_baseline/result/z.standard.Rdata")


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
# size = 5
# start.end <- startend(num,size,i2)
# start <- start.end[1]
# end <- start.end[2]
# file.num <- end-start+1
#num = 22349
#num <- countLines(geno.file)[1];
#num <- as.integer(system(paste0("zcat ",geno.file,"| wc -l"),intern=T))
rm(pheno)
num.of.tumor <- ncol(y.pheno.mis1)-1



score_result <- matrix(0.1,num,num.of.tumor)
infor_result <- matrix(0.1,(num.of.tumor)*num,num.of.tumor)
snpid_result <- rep("c",num)


freq.all <- rep(0,file.num)
#temp <- 0
con <- gzfile(geno.file)
open(con)
for(i in 1:num){
  if(i%%500==0){
    print(i)
  }
  oneLine <- readLines(con,n=1)
  
 # if(i>=start){
    #if(temp%%100==0){
     # print(paste0("temp",temp))
    #}
    #print(i)
 #   temp = temp+1
    #print(i)
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
          # score.test.support.icog.casecase <- ScoreTestSupportMixedModel(y=y.pheno.mis1,
          #                                                                baselineonly = snpvalue,
          #                                                                additive=x.all.covar,
          #                                                                missingTumorIndicator = 888,
          #                                                                delta0=delta0.icog)

          score.test.icog.casecase<- ScoreTestMixedModel(y=y.pheno.mis1,
                                                         x=snpvalue,
                                                         z.design = z.standard,

                                                         score.test.support= score.test.support.icog.casecase,
                                                         missingTumorIndicator=888)

          score_result[temp,]  <- score.test.icog.casecase[[1]]
          infor_result[((num.of.tumor)*temp-(num.of.tumor-1)):((num.of.tumor)*i),] <- score.test.icog.casecase[[2]]
          # print(mem_used())
          # rm(score.test.support.icog.casecase)
          rm(score.test.icog.casecase)
           gc()
          # print(memory.profile())
          # print(lsos())
          # print(mem_used())
          
        }
        
      },
      error=function(cond) {
        
        score_result[temp,] <- 0
        infor_result[((num.of.tumor)*temp-(num.of.tumor-1)):((num.of.tumor)*temp),] <- 0
        
        
      })
    
   # if(i==end){
    #  break
    #}
    
  #}
  
  
  
  
  
  
}
close(con)
# if(i !=num){
#   snpid_result <- snpid_result[1:i]
#   score_result <- score_result[1:i,]
#   infor_result <- infor_result[1:(num.of.tumor*i),]
#   
#   
# }
result <- list(snpid_reuslt=snpid_result,score_result=score_result,infor_result=infor_result,freq.all=freq.all)
save(result,file=paste0("./whole_genome/ICOG/ERPRHER2GRADE_casecase/result/ERPRHER2Grade_casecase_debug"))

