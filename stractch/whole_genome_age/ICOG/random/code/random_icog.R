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
###########n is the total number of people in imputed data file
n <- 109713
snpvalue <- rep(0,n)
######this is the SNP order in the imputed data file
subject.file <- "/gpfs/gsfs4/users/NC_BW/icogs_onco/genotype/imputed2/icogs_order.txt.gz"
Icog.order <- read.table(gzfile(subject.file))

setwd("/spin1/users/zhangh24/breast_cancer_data_analysis/")
data1 <- fread("./data/iCOGS_euro_v10_10232017.csv",header=T)
data1 <- as.data.frame(data1)
y.pheno.mis1 <- cbind(data1$Behaviour1,data1$ER_status1,data1$PR_status1,data1$HER2_status1,data1$Grade1)
colnames(y.pheno.mis1) = c("Behavior","ER","PR","HER2","Grade")
#x.test.all.mis1 <- data1[,c(27:206)]
SG_ID <- data1$SG_ID
########We need to incorporate age in the final model
########Some of the people don't have age, so I exclude them
x.covar.mis1 <- data1[,c(5:14,204)]
age <- data1[,204]
idx.complete <- which(age!=888)
y.pheno.mis1 <- y.pheno.mis1[idx.complete,]
x.covar.mis1 <- x.covar.mis1[idx.complete,]
SG_ID <- SG_ID[idx.complete]



###########we need to use idx.fil and idx.match to organize the SNP order as the phenotype data file
idx.fil <- Icog.order[,1]%in%SG_ID
idx.match <- match(SG_ID,Icog.order[idx.fil,1])
#Icog.order.match <- Icog.order[idx.fil,1][idx.match]
library(bc2)

#In real data analysis, I run scoresupport function in the a separate job. Then I load it here when I run whole genome score test. So it could save some time
#load("./whole_genome_age/ICOG/ERPRHER2GRADE_fixed_baseline/result/score.test.support.icog.ERPRHER2Grade.Rdata")
z.standard <- GenerateZstandard(y.pheno.mis1,missingTumorIndicator=888)
z.triple <- c(1,0,0,0,0,0,0,1,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0)
z.design <- cbind(1,z.standard,z.triple)
z.support <- z.design[,c(1,2,6)]
z.test <- z.design[,c(3:5)]

Filesdir <- "/gpfs/gsfs4/users/NC_BW/icogs_onco/genotype/imputed2/icogs_imputed/"
Files <- dir(Filesdir,pattern="icogs_merged_b1_12.",full.names=T)
Filesex <- dir(Filesdir,pattern="icogs_merged_b1_12.chr23",full.names=T)
idx.sex <- Files%in%Filesex
Files <- Files[!idx.sex]
library(gtools)
Files <- mixedsort(Files)
geno.file <- Files[i1]

####get the total number of length in imputed data
num <- as.integer(system(paste0("zcat ",geno.file,"| wc -l"),intern=T))


num.of.tumor <- ncol(y.pheno.mis1)-1
n.sub <- nrow(y.pheno.mis1)
idx.control <- which(y.pheno.mis1[,1]==0)
n.control <- length(idx.control)

library(doParallel)
library(foreach)
#######Here I parallel the job on two different node
no.cores <- 2
inner.size <- 2
registerDoParallel(no.cores)
result.list <- foreach(job.i = 1:inner.size)%dopar%{
  print(job.i)
  #startend function will automatically split the total lines into several equal parts based on how many jobs you want to parallel
  start.end <- startend(num,inner.size,job.i)
  ##for job i, you will start at start.end[1] and end at start.end[2]
  start <- start.end[1]
  end <- start.end[2]
  inner.num <- end-start+1
  ####creat matrix for score,infor, snpid and frequency
  score_result <- matrix(0,inner.num,num.of.tumor-1)
  infor_result <- matrix(0,inner.num,(num.of.tumor-1)^2)
  snpid_result <- rep("c",inner.num)
  freq.all <- rep(0,inner.num)
  temp <- 0
  #read the code line by line
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
      ######convert the imputed SNP into additive value
      snpvalue <- convert(snppro,n)
      snpvalue <- snpvalue[idx.fil][idx.match]
      snpvalue.control <- snpvalue[idx.control]
      #####calculate alleles frequency
      freq <- sum(snpvalue.control)/(2*n.control)
      freq.all[temp] <- freq
      #print(paste0("freq",freq))
      
      
      
      if(freq<0.006|freq>0.994){
        
        score_result[temp,] <- 0
        infor_result[temp,] <- 0.1
      }else{
        score.test.support.random.icog <- 
          score.test.random.icog<- ScoreTestSupportMixedModelSelfDesign            (y.pheno.mis1,
            x.self.design  = snpvalue,
            z.design = z.support,
            additive = x.covar.mis1,
            pairwise.interaction = NULL,
            saturated = NULL,
            missingTumorIndicator = 888
          )
        score.test.icog<- ScoreTestMixedModel(y=y.pheno.mis1,
                                              x=snpvalue,
                                              z.design=z.test,
                                              score.test.support=score.test.support.random.icog,
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





score_result <- matrix(0,num,num.of.tumor-1)
infor_result <- matrix(0,num,(num.of.tumor-1)^2)
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
save(result,file=paste0("./whole_genome_age/ICOG/ERPRHER2GRADE_fixed_baseline/result/ERPRHER2Grade_random_icog",i1))

