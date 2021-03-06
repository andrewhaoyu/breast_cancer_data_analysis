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

print(i1)
library(R.utils)
library(data.table)
setwd("/data/zhangh24/breast_cancer_data_analysis/")

n <- 109713
snpvalue <- rep(0,n)
subject.file <- "/gpfs/gsfs4/users/NC_BW/icogs_onco/genotype/imputed2/icogs_order.txt.gz"
Icog.order <- read.table(gzfile(subject.file))

setwd("/data/zhangh24/breast_cancer_data_analysis/")
data1 <- fread("./data/iCOGS_euro_v10_10232017.csv",header=T)
data1 <- as.data.frame(data1)
y.pheno.mis1 <- cbind(data1$Behaviour1,data1$PR_status1,data1$ER_status1,data1$HER2_status1,data1$Grade1)
colnames(y.pheno.mis1) = c("Behavior","PR","ER","HER2","Grade")
#x.test.all.mis1 <- data1[,c(27:206)]
SG_ID <- data1$SG_ID
x.covar.mis1 <- data1[,c(5:14,204)]
age <- data1[,204]
idx.complete <- which(age!=888)
y.pheno.mis1 <- y.pheno.mis1[idx.complete,]
x.covar.mis1 <- x.covar.mis1[idx.complete,]
SG_ID <- SG_ID[idx.complete]




idx.fil <- Icog.order[,1]%in%SG_ID
idx.match <- match(SG_ID,Icog.order[idx.fil,1])
#Icog.order.match <- Icog.order[idx.fil,1][idx.match]
library(bc2)
# load("./whole_genome_age/ICOG/ERPRHER2GRADE_fixed_baseline/result/score.test.support.icog.ERPRHER2Grade.Rdata")
# 
# 
# Filesdir <- "/gpfs/gsfs4/users/NC_BW/icogs_onco/genotype/imputed2/icogs_imputed/"
# Files <- dir(Filesdir,pattern="icogs_merged_b1_12.",full.names=T)
# Filesex <- dir(Filesdir,pattern="icogs_merged_b1_12.chr23",full.names=T)
# idx.sex <- Files%in%Filesex
# Files <- Files[!idx.sex]
# library(gtools)
# Files <- mixedsort(Files)
geno.file <- paste0("./genetic_correlation/ICOG/result/hapmap_icog",i1,".txt")

z.design <- matrix(c(
  c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1),
  c(0,0,0,0,0,1,1,0,0,0,0,0,1,1,1,0,0,0,0,0,1,1,1),
  c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0),
  c(0,0,0,0,1,0,0,0,0,0,0,1,0,0,0,0,0,0,0,1,0,0,0),
  c(1,0,0,0,0,0,0,1,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0)
),ncol=5)


colnames(z.design) <- c("Luminial A","Luminal B",
                        "Luminal B HER2Neg - Luminal A",
                        "HER2 Enriched - Luminal A",
                        "Triple Negative - Luminal A")

output <- system(paste0("wc -l ",geno.file),intern=T)
num <- as.integer(gsub(geno.file,"",output))

size = 1000
start.end <- startend(num,size,i2)
start <- start.end[1]
end <- start.end[2]
file.num <- end-start+1

if(num!=0){
  num.of.tumor <- ncol(y.pheno.mis1)-1
  n.sub <- nrow(y.pheno.mis1)
  idx.control <- which(y.pheno.mis1[,1]==0)
  n.control <- length(idx.control)
  
  #library(doParallel)
#  library(foreach)
  
  #no.cores <- 2
  #inner.size <- 2
  #registerDoParallel(no.cores)
 #result.list <- foreach(job.i = 1:inner.size)%dopar%{
    
    #inner.start.end <- startend(file.num,inner.size,job.i)
    #inner.start <- inner.start.end[1]
    #inner.end <- inner.start.end[2]
    inner.file.num <- end-start+1
    true.start <- start
    true.end <- end
    score_result <- matrix(0,inner.file.num,num.of.tumor+1)
    infor_result <- matrix(0,inner.file.num,(num.of.tumor+1)^2)
    snpid_result <- rep("c",inner.file.num)
    freq.all <- rep(0,inner.file.num)
    temp <- 0
    con <- gzfile(geno.file)
    open(con)
    # for(i in 1:399){
    #   print(i)
    #   oneLine <- readLines(con,n=1)
    #   if(i==399){
    #     myVector <- strsplit(oneLine," ")
    #     snpid <- as.character(myVector[[1]][3])
    #     print(snpid)
    #   }
    # }
    # for(i in 1:num){
      #print(i)
      if(i%%500==0){
        print(i)
      }
      oneLine <- readLines(con,n=1)
      
      if(i>=true.start){
        if(i==399){
          temp = temp+1
          score_result[temp,]  <- rep(NA,5)
          infor_result[temp,] <- rep(NA,25)
          
          
        }else{
          if(temp%%100==0){
            print(paste0("temp",temp))
          }
          #print(i)
          temp = temp+1
          #print(i)
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
          
          # tryCatch(
          #   {
          
          
          Heter.result.Icog = EMmvpolySelfDesign(y.pheno.mis1,x.self.design = snpvalue,z.design=z.design,baselineonly = NULL,additive = x.covar.mis1,pairwise.interaction = NULL,saturated = NULL,missingTumorIndicator = 888)
          z.standard <- Heter.result.Icog[[12]]
          M <- nrow(z.standard)
          number.of.tumor <- ncol(z.standard)
          log.odds.icog <- Heter.result.Icog[[1]][(M+1):(M+1+number.of.tumor)]
          nparm <- length(Heter.result.Icog[[1]])  
          sigma.log.odds.icog <- Heter.result.Icog[[2]][(M+1):(M+1+number.of.tumor),(M+1):(M+1+number.of.tumor)]
          
          score_result[temp,]  <- log.odds.icog
          infor_result[temp,] <- as.vector(sigma.log.odds.icog)
          
          
          
          
          
          
          
          
        }
        
        
        if(i==true.end){
          break
        }
        
        
        
        
        
        
        
        
        
        
        }
       
    
    
    
       
    }
    close(con)
    result <- list(snpid_result,score_result,infor_result,freq.all)
    #return(result)
  }
  
  #stopImplicitCluster()
  
  
  
  
  
  # score_result <- matrix(0,num,num.of.tumor+1)
  # infor_result <- matrix(0,num,(num.of.tumor+1)^2)
  # snpid_result <- rep("c",num)
  # freq.all <- rep(0,num)
  # 
  # total <- 0
  # for(i in 1:inner.size){
  #   result.temp <- result.list[[i]]
  #   temp <- length(result.temp[[1]])
  #   snpid_result[total+(1:temp)] <- result.temp[[1]]
  #   score_result[total+(1:temp),] <- result.temp[[2]]
  #   infor_result[total+(1:temp),] <- result.temp[[3]]
  #   freq.all[total+(1:temp)] <- result.temp[[4]]
  #   total <- total+temp
  # }
  # 
  # 
  # 
  
  
 # result <- list(snpid_reuslt=snpid_result,score_result=score_result,infor_result=infor_result,freq.all=freq.all)
  save(result,file=paste0("./genetic_correlation/ICOG/result/intrinsic_i1",i1,"_",i2))
  
  
