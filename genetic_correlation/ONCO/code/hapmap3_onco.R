rm(list=ls())#args=(commandArgs(TRUE))
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
#i1 <- i
print(i1)
library(R.utils)
library(data.table)
setwd("/spin1/users/zhangh24/breast_cancer_data_analysis/")

subject.file <- "/gpfs/gsfs4/users/NC_BW/icogs_onco/genotype/imputed2/onco_order.txt.gz"
onco.order <- read.table(gzfile(subject.file))
setwd("/spin1/users/zhangh24/breast_cancer_data_analysis/")
data2 <- fread("./data/Onco_euro_v10_10232017.csv",header=T)
data2 <- as.data.frame(data2)
y.pheno.mis2 <- cbind(data2$Behaviour1,data2$PR_status1,data2$ER_status1,data2$HER2_status1,data2$Grade1)
#y.pheno.mis2 <- cbind(data2$Behaviour1,data2$PR_status1,data2$ER_status1,data2$HER2_status1)
colnames(y.pheno.mis2) = c("Behaviour","PR",
                           "ER","HER2","Grade")

x.covar.mis2 <- data2[,c(5:14,204)]
ages <- data2[,204]
idx.complete <- which(ages!=888)

y.pheno.mis2 <- y.pheno.mis2[idx.complete,]
x.covar.mis2 <- x.covar.mis2[idx.complete,]
Onc_ID <- data2$Onc_ID
Onc_ID <- Onc_ID[idx.complete]
rm(data2)
gc()



idx.fil <- onco.order[,1]%in%Onc_ID
n <- length(idx.fil)
snpvalue <- rep(0,n)
idx.match <- match(Onc_ID,onco.order[idx.fil,1])
#Icog.order.match <- Icog.order[idx.fil,1][idx.match]
library(bc2)


geno.file <- paste0("./genetic_correlation/ONCO/result/hapmap_onco",i1,".txt")
z.design <- matrix(c(
  c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1),
  c(0,0,0,0,0,1,1,0,0,0,0,0,1,1,1,0,0,0,0,0,1,1,1),
  c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0),
  c(0,0,0,0,1,0,0,0,0,0,0,1,0,0,0,0,0,0,0,1,0,0,0),
  c(1,0,0,0,0,0,0,1,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0)
),ncol=5)
output <- system(paste0("wc -l ",geno.file),intern=T)
num <- as.integer(gsub(geno.file,"",output))

if(num!=0){
  num.of.tumor <- ncol(y.pheno.mis2)-1
  n.sub <- nrow(y.pheno.mis2)
  idx.control <- which(y.pheno.mis2[,1]==0)
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
      if(i%%500==0){
        print(i)
      }
      oneLine <- readLines(con,n=1)
      
      if(i>=start){
        temp <- temp+1
        myVector <- strsplit(oneLine," ")
        snpid <- as.character(myVector[[1]][3])
        snpid_result[temp] <- snpid
        snpvalue <- rep(0,n)
        
        
        snppro <- as.numeric(unlist(myVector)[7:length(myVector[[1]])])
        if(length(snppro)!=(3*n)){
          break
        }
        
        snpvalue <- convert(snppro,n)
        snpvalue <- snpvalue[idx.fil][idx.match]
        snpvalue.control <- snpvalue[idx.control]
        freq <- sum(snpvalue.control)/(2*n.control)
        freq.all[temp] <- freq
        #print(paste0("freq",freq))
        
        
        
        
        
        Heter.result.Onco = EMmvpolySelfDesign(y.pheno.mis2,x.self.design = snpvalue,z.design = z.design,baselineonly = NULL,additive = x.covar.mis2,pairwise.interaction = NULL,saturated = NULL,missingTumorIndicator = 888)
        z.standard <- Heter.result.Onco[[12]]
        M <- nrow(z.standard)
        number.of.tumor <- ncol(z.standard)
        log.odds.onco <- Heter.result.Onco[[1]][(M+1):(M+1+number.of.tumor)]
        nparm <- length(Heter.result.Onco[[1]])
        sigma.log.odds.onco <- Heter.result.Onco[[2]][(M+1):(M+1+number.of.tumor),(M+1):(M+1+number.of.tumor)]
        
        score_result[temp,]  <- log.odds.onco
        infor_result[temp,] <- as.vector(sigma.log.odds.onco)
        
        
        
        
        
        
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
  save(result,file=paste0("./genetic_correlation/ONCO/result/intrinsic_i1",i1))
  
}else{
result <- NULL
save(result,file=paste0("./genetic_correlation/ONCO/result/intrinsic_i1",i1))
}


