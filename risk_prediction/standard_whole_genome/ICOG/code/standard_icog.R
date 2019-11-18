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
library(R.utils)
library(data.table)
setwd("/data/zhangh24/breast_cancer_data_analysis/")

n <- 109713
snpvalue <- rep(0,n)
subject.file <- "/gpfs/gsfs4/users/NC_BW/icogs_onco/genotype/imputed2/icogs_order.txt.gz"
Icog.order <- read.table(gzfile(subject.file))

setwd("/data/zhangh24/breast_cancer_data_analysis/")
load(paste0("./risk_prediction/result/split.id.rdata"))
#icog.test.id <- Generatetestid(subtypes.icog)
icog.train.id <- split.id[[1]]
#onco.train.id <- split.id[[2]]
#onco.test.id <- split.id[[3]]
#icog.cohort.id <- split.id[[4]]
#onco.cohort.id <- split.id[[5]]
#Icog.order <- read.table(gzfile(subject.file))
######load in the data and take out the training data
data1 <- as.data.frame(fread("/data/zhangh24/breast_cancer_data_analysis/data/sig_snp_icog_prs.csv",header=T))
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
#x.covar.mis1 <- data1[,c(7:16)]

table(y.pheno.mis1[,1])
colnames(y.pheno.mis1) = c("Behavior","ER","PR","HER2","Grade")
#x.test.all.mis1 <- data1[,c(27:206)]
SG_ID <- data1$ID
x.covar.mis1 <- as.matrix(data1[,c(7:16)])




idx.fil <- Icog.order[,1]%in%SG_ID
idx.match <- match(SG_ID,Icog.order[idx.fil,1])
#Icog.order.match <- Icog.order[idx.fil,1][idx.match]
library(bc2, lib.loc ="/home/zhangh24/R/x86_64-pc-linux-gnu-library/3.6/")

#load("./risk_prediction/FTOP_whole_genome/ICOG/result/score.test.support.icog.ERPRHER2Grade.Rdata")


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

size = 5
start.end <- startend(num,size,i2)
start <- start.end[1]
end <- start.end[2]
file.num <- end-start+1
num.of.tumor = 0

no.cores <- 2
inner.size <- 2
registerDoParallel(no.cores)
result.list <- foreach(job.i = 1:inner.size)%dopar%{
  inner.start.end <- startend(file.num,inner.size,job.i)
  inner.start <- inner.start.end[1]
  inner.end <- inner.start.end[2]
  inner.file.num <- inner.end-inner.start+1
  true.start <- start+inner.start-1
  true.end <- start+inner.end-1
  score_result <- matrix(0,inner.file.num,num.of.tumor+1)
  infor_result <- matrix(0,inner.file.num,(num.of.tumor+1)^2)
  snpid_result <- rep("c",inner.file.num)
  freq.all <- rep(0,inner.file.num)
  temp <- 0
  con <- gzfile(geno.file)
  open(con)
  for(i in 1:num){
   #if(i%%500==0){
    print(i)
    #}
    oneLine <- readLines(con,n=1)
    
    if(i>=true.start){
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
      

          
          if(freq<0.008|freq>0.992){
            
            score_result[temp,] <- 0
            infor_result[temp,] <- 100
          }else{
          
            model.icog<- glm(y.pheno.mis1[,1]~snpvalue+x.covar.mis1,family = binomial())
            summmary.icog <- summary(model.icog)
             
            
            score_result[temp]  <- summmary.icog$coefficients[2,1]
            infor_result[temp] <- summmary.icog$coefficients[2,2]^2
            
            
          }
          
   
      
    }
   
    if(i==true.end){
      break
    }
    
  }
  close(con)
  result <- list(snpid_result,score_result,infor_result,freq.all)
  return(result)
  }

stopImplicitCluster()

num.of.tumor <- 0

score_result <- matrix(0.1,file.num,num.of.tumor+1)
infor_result <- matrix(0.1,file.num,(num.of.tumor+1)^2)
snpid_result <- rep("c",file.num)

freq.all <- rep(0,file.num)



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
save(result,file=paste0("./risk_prediction/standard_whole_genome/ICOG/result/standard",i1,"_",i2))

