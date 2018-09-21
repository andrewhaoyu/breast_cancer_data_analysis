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

print(i1)
library(R.utils)
library(data.table)
setwd("/spin1/users/zhangh24/breast_cancer_data_analysis/")

n <- 109713
snpvalue <- rep(0,n)
subject.file <- "/gpfs/gsfs4/users/NC_BW/icogs_onco/genotype/imputed2/icogs_order.txt.gz"
Icog.order <- read.table(gzfile(subject.file))

setwd("/spin1/users/zhangh24/breast_cancer_data_analysis/")
######load in the data and take out the training data
data1 <- as.data.frame(fread("/spin1/users/zhangh24/breast_cancer_data_analysis/data/sig_snp_icog_prs.csv",header=T))
data1 <- as.data.frame(data1[,-1])
icog.train <- which(data1[,1]%in%icog.train.id)
data1 <- data1[icog.train,]
y.pheno.mis1 <- cbind(data1$Behavior,data1$ER,data1$PR,data1$HER2,data1$Grade)
##########clean phenotype file
idx <- which(y.pheno.mis1[,1]==888|y.pheno.mis1[,1]==2)
y.pheno.mis1[idx,1] <- 1
colnames(y.pheno.mis1) = c("Behavior","ER","PR","HER2","Grade")
#x.test.all.mis1 <- data1[,c(27:206)]
SG_ID <- data1$ID
x.covar.mis1 <- data1[,c(7:16)]

table(y.pheno.mis1[,1])
colnames(y.pheno.mis1) = c("Behavior","ER","PR","HER2","Grade")
#x.test.all.mis1 <- data1[,c(27:206)]
SG_ID <- data1$ID
x.covar.mis1 <- data1[,c(7:16)]




idx.fil <- Icog.order[,1]%in%SG_ID
idx.match <- match(SG_ID,Icog.order[idx.fil,1])
#Icog.order.match <- Icog.order[idx.fil,1][idx.match]
library(bc2)
load("./risk_prediction/FTOP_whole_genome/ICOG/result/score.test.support.icog.ERPRHER2Grade.Rdata")


Filesdir <- "/gpfs/gsfs4/users/NC_BW/icogs_onco/genotype/imputed2/icogs_imputed/"
Files <- dir(Filesdir,pattern="icogs_merged_b1_12.",full.names=T)
Filesex <- dir(Filesdir,pattern="icogs_merged_b1_12.chr23",full.names=T)
idx.sex <- Files%in%Filesex
Files <- Files[!idx.sex]
library(gtools)
Files <- mixedsort(Files)
geno.file <- Files[i1]


num <- as.integer(system(paste0("zcat ",geno.file,"| wc -l"),intern=T))


num.of.tumor <- ncol(y.pheno.mis1)-1
n.sub <- nrow(y.pheno.mis1)
idx.control <- which(y.pheno.mis1[,1]==0)
n.control <- length(idx.control)

library(doParallel)
library(foreach)

no.cores <- 2
inner.size <- 2
registerDoParallel(no.cores)
result.list <- foreach(job.i = 1:inner.size)%dopar%{
  print(job.i)
  start.end <- startend(num,inner.size,job.i)
  start <- start.end[1]
  end <- start.end[2]
  inner.num <- end-start+1
  score_result <- matrix(0,inner.num,num.of.tumor+1)
  infor_result <- matrix(0,inner.num,(num.of.tumor+1)^2)
  snpid_result <- rep("c",inner.num)
  freq.all <- rep(0,inner.num)
  temp <- 0
  con <- gzfile(geno.file)
  open(con)
  for(i in 1:num){
   #if(i%%500==0){
    print(i)
    #}
    oneLine <- readLines(con,n=1)
    
    if(i>=start){
      temp <- temp+1
      myVector <- strsplit(oneLine," ")
      snpid <- as.character(myVector[[1]][2])
      snpid_result[temp] <- snpid
      snpvalue <- rep(0,n)
      
      
      snppro <- as.numeric(unlist(myVector)[6:length(myVector[[1]])])
      if(length(snppro)!=(3*n)){
        break
      }
      
      snpvalue <- convert(snppro,n)
      snpvalue <- snpvalue[idx.fil][idx.match]
      snpvalue.control <- snpvalue[idx.control]
      freq <- sum(snpvalue.control)/(2*n.control)
      freq.all[temp] <- freq
      #print(paste0("freq",freq))
      

          
          if(freq<0.006|freq>0.994){
            
            score_result[temp,] <- 0
            infor_result[temp,] <- 0.1
          }else{
            
            score.test.icog<- ScoreTest(y=y.pheno.mis1,
                                        x=snpvalue,
                                        second.stage.structure="additive",
                                        score.test.support=score.test.support.icog.ERPRHER2Grade,
                                        missingTumorIndicator=888)
            
            score_result[temp,]  <- score.test.icog[[1]]
            infor_result[temp,] <- as.vector(score.test.icog[[2]])
            
            
          }
          
   
      
    }
   
    if(i==end){
      break
    }
    
  }
  close(con)
  result <- list(snpid_result,score_result,infor_result,freq.all)
  return(result)
  }

stopImplicitCluster()





score_result <- matrix(0,num,num.of.tumor+1)
infor_result <- matrix(0,num,(num.of.tumor+1)^2)
snpid_result <- rep("c",num)
freq.all <- rep(0,num)

total <- 0
for(i in 1:inner.size){
  result.temp <- result.list[[i]]
  temp <- length(result.temp[[1]])
  snpid_result[total+(1:temp)] <- result.temp[[1]]
  score_result[total+(1:temp),] <- result.temp[[2]]
  infor_result[total+(1:temp),] <- result.temp[[3]]
  freq.all[total+(1:temp)] <- result.temp[[4]]
   total <- total+temp
}




result <- list(snpid_reuslt=snpid_result,score_result=score_result,infor_result=infor_result,freq.all=freq.all)
save(result,file=paste0("./risk_prediction/FTOP_whole_genome/ICOG/result/ERPRHER2Grade_fixed_baseline",i1))

