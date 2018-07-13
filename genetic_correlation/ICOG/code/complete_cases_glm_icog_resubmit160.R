#############ICOG job 160 get some problem with line 399
#############resubmit the iCOG jobs 160 with 1000 subjobs
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
library(data.table)
setwd("/spin1/users/zhangh24/breast_cancer_data_analysis/")

n <- 109713
snpvalue <- rep(0,n)
subject.file <- "/gpfs/gsfs4/users/NC_BW/icogs_onco/genotype/imputed2/icogs_order.txt.gz"
Icog.order <- read.table(gzfile(subject.file))

setwd("/spin1/users/zhangh24/breast_cancer_data_analysis/")
data1 <- fread("./data/iCOGS_euro_v10_10232017.csv",header=T)
data1 <- as.data.frame(data1)

enriched.study <- c("2SISTER","ABCS-F","BBCS","BCFR-NY", "BCFR-PA", "BCFR-UT", "BOCS", "CNIO-BCS", "FHRISK","GC-HBOC", "HEBCS", "HEBON", "HKBCS", "IPOBCS", "KARBAC", "kConFab/AOCS", "KOHBRA", "MBCSG", "MSKCC","MYBRCA", "NBCS", "NC-BCFR", "OFBCR", "RBCS", "SUCCESSB", "SUCCESSC")
idx.enriched <- which(data1$study%in%enriched.study)
data1 <- data1[-idx.enriched,]
y.pheno.mis1 <- cbind(data1$Behaviour1,data1$ER_status1,data1$PR_status1,data1$HER2_status1,data1$Grade1)
idx.mis <- which(data1$ER_status1==888|data1$PR_status1==888|data1$HER2_status1==888|data1$Grade1==888)
data1 <- data1[-idx.mis,]
y.pheno.mis1 <- cbind(data1$Behaviour1,data1$ER_status1,data1$PR_status1,data1$HER2_status1,data1$Grade1)
colnames(y.pheno.mis1) = c("Behavior","ER","PR","HER2","Grade")
#x.test.all.mis1 <- data1[,c(27:206)]
SG_ID <- data1$SG_ID
x.covar.mis1 <- data1[,c(5:14,204)]
age <- data1[,204]
idx.complete <- which(age!=888)
y.pheno.mis1 <- y.pheno.mis1[idx.complete,]
x.covar.mis1 <- x.covar.mis1[idx.complete,]
SG_ID <- SG_ID[idx.complete]

GenerateIntrinsic <- function(ER,PR,HER2,Grade){
  n <- length(ER)
  idx.LA <- which(HER2==0&(ER==1|PR==1)&Grade!=3)
  idx.LB <- which(HER2==1&(ER==1|PR==1))
  idx.LUBHER2 <- which(HER2==0&(ER==1|PR==1)&Grade==3)
  idx.HER2 <- which(HER2==1&ER==0&PR==0)
  idx.Tp <- which(HER2==0&ER==0&PR==0)
  subtypes <- rep("control",n)
  subtypes[idx.LA] <- "Luminal_A"
  subtypes[idx.LB] <- "Luminal_B"
  subtypes[idx.LUBHER2] <- "Luminal_B_HER2Enriched"
  subtypes[idx.HER2] <- "HER2Enriched"
  subtypes[idx.Tp] <- "TripleNeg"
  subtypes <- factor(subtypes,levels=c("control",
                                       "Luminal_A",
                                       "Luminal_B",
                                       "Luminal_B_HER2Enriched",
                                       "HER2Enriched",
                                       "TripleNeg"))
  return(subtypes)
}

subtypes <- GenerateIntrinsic(y.pheno.mis1[,2],y.pheno.mis1[,3],y.pheno.mis1[,4],y.pheno.mis1[,5])

library(nnet)


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



output <- system(paste0("wc -l ",geno.file),intern=T)
num <- as.integer(gsub(geno.file,"",output))
size = 1000
start.end <- startend(num,size,i2)
start <- start.end[1]
end <- start.end[2]
file.num <- end-start+1


#if(num!=0){
  num.of.tumor <- ncol(y.pheno.mis1)-1
  n.sub <- nrow(y.pheno.mis1)
  idx.control <- which(y.pheno.mis1[,1]==0)
  n.control <- length(idx.control)
  
  # library(doParallel)
  # library(foreach)
  # 
  # no.cores <- 2
  # inner.size <- 2
  # registerDoParallel(no.cores)
  # result.list <- foreach(job.i = 1:inner.size)%dopar%{
  #   print(job.i)
    # start.end <- startend(num,inner.size,job.i)
    # start <- start.end[1]
    # end <- start.end[2]
    inner.num <- end-start+1
    true.start <- start
    true.end <- end
    score_result <- matrix(0,inner.num,num.of.tumor+1)
    infor_result <- matrix(0,inner.num,(num.of.tumor+1)^2)
    snpid_result <- rep("c",inner.num)
    freq.all <- rep(0,inner.num)
    temp <- 0
    con <- file(geno.file)
    open(con)
    for(i in 1:num){
      #if(i%%500==0){
      print(i)
      #}
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
          model1 <- multinom(subtypes~snpvalue+as.matrix(x.covar.mis1))
          score_result[temp,]  <- as.vector(coef(model1)[,2])
          infor_result[temp,] <- as.vector(vcov(model1)[2+13*(0:4),2+13*(0:4)])
          
          
          
          
          
          
        }
        
        if(i==true.end){
          break
        }
        
      
      }
    }
    close(con)
    result <- list(snpid_result,score_result,infor_result,freq.all)
    
  
  save(result,file=paste0("./genetic_correlation/ICOG/result/complete_glm",i1,"_",i2))
  
  
#}