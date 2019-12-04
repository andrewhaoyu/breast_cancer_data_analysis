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
#i1 <- i
print(i1)
library(R.utils)
library(data.table)
setwd("/data/zhangh24/breast_cancer_data_analysis/")
load(paste0("./risk_prediction/result/split.id.rdata"))
#icog.test.id <- Generatetestid(subtypes.icog)
#icog.train.id <- split.id[[1]]
onco.train.id <- split.id[[2]]
#onco.test.id <- split.id[[3]]
#icog.cohort.id <- split.id[[4]]
#onco.cohort.id <- split.id[[5]]
subject.file <- "/gpfs/gsfs4/users/NC_BW/icogs_onco/genotype/imputed2/onco_order.txt.gz"
onco.order <- read.table(gzfile(subject.file))
library(bc2, lib.loc ="/home/zhangh24/R/x86_64-pc-linux-gnu-library/3.6/")

setwd("/data/zhangh24/breast_cancer_data_analysis/")
data2 <- as.data.frame(fread("/data/zhangh24/breast_cancer_data_analysis/data/sig_snp_onco_prs.csv",header=T))
data2 <- data2[,-1]
#onco.train <- which(data2[,1]%in%onco.train.id)
onco.train <- which(data2[,1]%in%onco.train.id)
data2 <- data2[onco.train,]
y.pheno.mis2 <- cbind(data2$Behavior,data2$ER,data2$PR,data2$HER2,data2$Grade)
#######clean y.phneo.mis2 
idx <- which(y.pheno.mis2[,1]==888|y.pheno.mis2[,1]==2)
y.pheno.mis2[idx,1] <- 1
#y.pheno.mis2 <- cbind(data2$Behaviour1,data2$PR_status1,data2$ER_status1,data2$HER2_status1)
colnames(y.pheno.mis2) = c("Behaviour","ER",
                           "PR","HER2","Grade")
z.additive <- cbind(1,GenerateZstandard(y.pheno.mis2))
x.covar.mis2 <- data2[,c(7:16)]
#ages <- data2[,230]
#idx.complete <- which(ages!=888)

#y.pheno.mis2 <- y.pheno.mis2[idx.complete,]
#x.covar.mis2 <- x.covar.mis2[idx.complete,]
Onc_ID <- data2$ID
#Onc_ID <- Onc_ID[idx.complete]
rm(data2)
gc()



idx.fil <- onco.order[,1]%in%Onc_ID
n <- length(idx.fil)
snpvalue <- rep(0,n)
idx.match <- match(Onc_ID,onco.order[idx.fil,1])
#Icog.order.match <- Icog.order[idx.fil,1][idx.match]
load("./risk_prediction/FTOP_whole_genome/ONCO/result/score.test.support.onco.ERPRHER2Grade.Rdata")





Filesdir <- "/gpfs/gsfs4/users/NC_BW/icogs_onco/genotype/imputed2/onco_imputed"
Files <- dir(Filesdir,pattern="OncoArray_european_merged_b1_15.",full.names=T)
Filesex <- dir(Filesdir,pattern="OncoArray_european_merged_b1_15.chr23",full.names=T)
idx.sex <- Files%in%Filesex
library(gtools)
Files <- Files[!idx.sex]
Files <- mixedsort(Files)
geno.file <- Files[i1]

    num <- as.integer(system(paste0("zcat ",geno.file,"| wc -l"),intern=T))




    num.of.tumor <- ncol(y.pheno.mis2)-1
    n.sub <- nrow(y.pheno.mis2)
    idx.control <- which(y.pheno.mis2[,1]==0)
    n.control <- length(idx.control)
    
    library(doParallel)
    library(foreach)
    
    no.cores <- 2
    inner.size <- 2
    size = 6
    start.end <- startend(num,size,i2)
    start <- start.end[1]
    end <- start.end[2]
    file.num <- end-start+1
    
    
    registerDoParallel(no.cores)
    result.list <- foreach(job.i = 1:inner.size)%dopar%{
      print(job.i)
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
      temp <- 0
      con <- gzfile(geno.file)
      open(con)
      for(i in 1:num){
      #  if(i%%500==0){
        print(i)
       # }
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
          
          
          
          if(freq<0.006|freq>0.994){
            
            score_result[temp,] <- 0
            infor_result[temp,] <- 0.1
          }else{
            
            score.test.onco<- ScoreTestMixedModel(y=y.pheno.mis2,
x=snpvalue,
                                        z.design = z.additive,
                                        score.test.support=score.test.support.onco.ERPRHER2Grade,
                                        missingTumorIndicator=888)
            
            score_result[temp,]  <- score.test.onco[[1]]
            infor_result[temp,] <- as.vector(score.test.onco[[2]])
            
            
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
    
    
    
    
    
    score_result <- matrix(0,file.num,num.of.tumor+1)
    infor_result <- matrix(0,file.num,(num.of.tumor+1)^2)
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
save(result,file=paste0("./risk_prediction/FTOP_whole_genome/ONCO/result/ERPRHER2Grade_fixed_onco_120419_",i1,"_",i2))
