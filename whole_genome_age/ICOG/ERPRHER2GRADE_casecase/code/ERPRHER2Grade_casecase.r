rm(list=ls())
#install_github("andrewhaoyu/bc2", ref = "master",args = c('--library="/home/zhangh24/R/x86_64-pc-linux-gnu-library/3.4"'))
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

n <- 109713
snpvalue <- rep(0,n)
subject.file <- "/gpfs/gsfs4/users/NC_BW/icogs_onco/genotype/imputed2/icogs_order.txt.gz"

Icog.order <- read.table(gzfile(subject.file))
library(data.table)
setwd("/spin1/users/zhangh24/breast_cancer_data_analysis/")
data1 <- fread("./data/iCOGS_euro_v10_10232017.csv",header=T)
data1 <- as.data.frame(data1)
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
rm(data1)
gc()




idx.fil <- Icog.order[,1]%in%SG_ID
idx.match <- match(SG_ID,Icog.order[idx.fil,1])
#Icog.order.match <- Icog.order[idx.fil,1][idx.match]
library(bc2)
load("./whole_genome_age/ICOG/ERPRHER2GRADE_fixed_baseline/result/delta0.icog.Rdata")
load("./whole_genome_age/ICOG/ERPRHER2GRADE_fixed_baseline/result/z.standard.Rdata")
z.design.support <- cbind(1,z.standard[,1])

Filesdir <- "/gpfs/gsfs4/users/NC_BW/icogs_onco/genotype/imputed2/icogs_imputed/"
Files <- dir(Filesdir,pattern="icogs_merged_b1_12.",full.names=T)
Filesex <- dir(Filesdir,pattern="icogs_merged_b1_12.chr23",full.names=T)
idx.sex <- Files%in%Filesex
Files <- Files[idx.sex]
library(gtools)
Files <- mixedsort(Files)
geno.file <- Files[i1]

# tryCatch(
#   {
    num <- as.integer(system(paste0("zcat ",geno.file,"| wc -l"),intern=T))
#   },
#   error=function(cond){
#     num <- countLines(geno.file)[1]
#   }
# )
size = 5
start.end <- startend(num,size,i2)
start <- start.end[1]
end <- start.end[2]
file.num <- end-start+1
#num = 22349
#num <- countLines(geno.file)[1];
#num <- as.integer(system(paste0("zcat ",geno.file,"| wc -l"),intern=T))
#rm(pheno)
num.of.tumor <- ncol(y.pheno.mis1)-1
n.sub <- nrow(y.pheno.mis1)
idx.control <- which(y.pheno.mis1[,1]==0)
n.control <- length(idx.control)


z.design.test <- z.standard[,2:4]


no.cores <- 2
library(foreach)
library(doParallel)

inner.size <- 2

registerDoParallel(no.cores)

result.list <- foreach(job.i = 1:2)%dopar%{
  inner.start.end <- startend(file.num,inner.size,job.i)
  inner.start <- inner.start.end[1]
  inner.end <- inner.start.end[2]
  inner.file.num <- inner.end-inner.start+1
  true.start <- start+inner.start-1
  true.end <- start+inner.end-1
  score_result <- matrix(0,inner.file.num,num.of.tumor-1)
  infor_result <- matrix(0,inner.file.num,(num.of.tumor-1)^2)
  snpid_result <- rep("c",inner.file.num)
  freq.all <- rep(0,inner.file.num)
  temp <- 0
  con <- gzfile(geno.file)
  open(con)
  for(i in 1:num){
    if(i%%500==0){
      print(i)
    }
    oneLine <- readLines(con,n=1)
    
    if(i>=true.start){
      if(temp%%100==0){
        print(paste0("temp",temp))
      }
      #print(i)
      temp = temp+1
      #print(i)
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
      
      # tryCatch(
      #   {
          
          if(freq<0.006|freq>0.994){
            
            score_result[temp,] <- 0
            infor_result[temp,] <- 0
          }else{
            score.test.support.icog.casecase <- ScoreTestSupportMixedModelSelfDesign(y=y.pheno.mis1,
                                                                                     x.self.design = snpvalue,
                                                                                     z.design = z.design.support,
                                                                           additive=x.covar.mis1,
                                                                           missingTumorIndicator = 888,
                                                                           delta0=delta0)
            
            score.test.icog.casecase<- ScoreTestMixedModel(y=y.pheno.mis1,
                                                           x=snpvalue,
                                                           z.design = z.design.test,
                                                           
                                                           score.test.support= score.test.support.icog.casecase,
                                                           missingTumorIndicator=888)
            
            score_result[temp,]  <- score.test.icog.casecase[[1]]
            infor_result[temp,] <- as.vector(score.test.icog.casecase[[2]])
           
            
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

score_result <- matrix(0.1,file.num,num.of.tumor-1)
infor_result <- matrix(0.1,file.num,(num.of.tumor-1)^2)
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
  total <- temp+total
}

result <- list(snpid_reuslt=snpid_result,score_result=score_result,infor_result=infor_result,freq.all=freq.all)

save(result,file=paste0("./whole_genome_age/ICOG/ERPRHER2GRADE_casecase/result/ERPRHER2Grade_casecase",i1,"_",i2))

