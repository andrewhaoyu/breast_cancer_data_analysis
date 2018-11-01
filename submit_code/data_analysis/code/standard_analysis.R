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

Generatesubtypes<- function(ER,PR,HER2,Grade){
  n <- length(ER)
  subtypes <- rep("control",n)
  temp = 1
  for(i in 0:1){
    for(j in 0:1){
      for(k in 0:1){
        for(l in 1:3){
          
          idx <- which(ER==i&PR==j&HER2==k&Grade==l)
          if(length(idx)!=0){
            subtypes[idx] <- temp
            temp = temp+1  
          }
          
        }
      }
    }
  }
  
  subtypes <- factor(subtypes,levels=c("control",
                                       c(1:temp)))
  sum <- table(subtypes)
  idx.cat <- which(sum<=10)
  idx.remove <- which((subtypes%in%(unique(idx.cat)-1))==T)
  subtypes <- subtypes[-idx.remove]
  return(list(subtypes,idx.remove))
}


library(R.utils)
library(data.table)
library(nnet)
setwd("/spin1/users/zhangh24/breast_cancer_data_analysis/")

subject.file <- "/gpfs/gsfs4/users/NC_BW/icogs_onco/genotype/imputed2/onco_order.txt.gz"
onco.order <- read.table(gzfile(subject.file))
setwd("/spin1/users/zhangh24/breast_cancer_data_analysis/")
data2 <- fread("./data/Onco_euro_v10_10232017.csv",header=T)
data2 <- as.data.frame(data2)
idx <- which(data2$StudyCountry=="Poland")
data2 <- data2[idx,]

y.pheno.mis2 <- cbind(data2$Behaviour1,data2$ER_status1,data2$PR_status1,data2$HER2_status1,data2$Grade1)
#y.pheno.mis2 <- cbind(data2$Behaviour1,data2$PR_status1,data2$ER_status1,data2$HER2_status1)
colnames(y.pheno.mis2) = c("Behaviour","ER",
                           "PR","HER2","Grade")

x.covar.mis2 <- data2[,c(5:8)]

Onc_ID <- data2$Onc_ID
# 
# load("./poland/result/whole_genome/delta0.onco.Rdata")
# load("./poland/result/whole_genome/z.standard.Rdata")
# z.design.support <- cbind(1,z.standard[,1])
# z.design.test <- z.standard[,2:4]
#Onc_ID <- data2$Onc_ID



idx.mis <- which(data2$ER_status1==888|data2$PR_status1==888|
               data2$HER2_status1==888|data2$Grade1==888)
data2.com <- data2[-idx.mis,]
x.covar.com <- data2[-idx.mis,c(5:8)]
Onc_ID.com <-  data2$Onc_ID[-idx.mis]
y.pheno.com <- cbind(data2.com$Behaviour1,data2.com$ER_status1,data2.com$PR_status1,data2.com$HER2_status1,data2.com$Grade1)

#y.pheno.mis2 <- cbind(data2$Behaviour1,data2$PR_status1,data2$ER_status1,data2$HER2_status1)
colnames(y.pheno.com) = c("Behaviour","ER",
                           "PR","HER2","Grade")



temp <-  Generatesubtypes(data2.com$ER_status1,data2.com$PR_status1,data2.com$HER2_status1,data2.com$Grade1)
subtypes <- temp[[1]]
idx.remove <- temp[[2]]
Onc_ID2 <- Onc_ID.com[-idx.remove]
x.covar.complete <- x.covar.com[-idx.remove,]









rm(data2)
 gc()
# load("./poland/result/whole_genome/score.test.support.onco.ERPRHER2Grade.Rdata")


idx.fil <- onco.order[,1]%in%Onc_ID
idx.fil2 <- onco.order[,1]%in%Onc_ID2
n <- length(idx.fil)
snpvalue <- rep(0,n)
idx.match <- match(Onc_ID,onco.order[idx.fil,1])
idx.match2 <- match(Onc_ID2,onco.order[idx.fil2,1])
#Icog.order.match <- Icog.order[idx.fil,1][idx.match]
library(bc2)
#load("./whole_genome/ONCO/ERPRHER2GRADE_fixed_baseline/result/score.test.support.onco.ERPRHER2Grade.Rdata")








Filesdir <- "/gpfs/gsfs4/users/NC_BW/icogs_onco/genotype/imputed2/onco_imputed"
Files <- dir(Filesdir,pattern="OncoArray_european_merged_b1_15.",full.names=T)
Filesex <- dir(Filesdir,pattern="OncoArray_european_merged_b1_15.chr23",full.names=T)
idx.sex <- Files%in%Filesex
Files <- Files[!idx.sex]
library(gtools)
Files <- mixedsort(Files)
geno.file <- Files[i1]

num <- as.integer(system(paste0("zcat ",geno.file,"| wc -l"),intern=T))


num.of.tumor <- ncol(y.pheno.mis2)-1

n.sub <- nrow(y.pheno.mis2)
idx.control <- which(y.pheno.mis2[,1]==0)
n.control <- length(idx.control)


size = 5
start.end <- startend(num,size,i2)
start <- start.end[1]
end <- start.end[2]
file.num <- end-start+1


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
  p_value1 <- rep(0,inner.file.num)
  p_value2 <- rep(0,inner.file.num)
  # score_result <- matrix(0,inner.file.num,num.of.tumor+1)
  # infor_result <- matrix(0,inner.file.num,(num.of.tumor+1)^2)
  # score_result2 <- matrix(0,inner.file.num,num.of.tumor-1)
  # infor_result2 <- matrix(0,inner.file.num,(num.of.tumor-1)^2)
  # 
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
      
      rawsnpvalue <- convert(snppro,n)
      snpvalue <- rawsnpvalue[idx.fil][idx.match]
      snpvalue2 <- rawsnpvalue[idx.fil2][idx.match2]
      snpvalue.control <- snpvalue[idx.control]
      freq <- sum(snpvalue.control)/(2*n.control)
      freq.all[temp] <- freq
      #print(paste0("freq",freq))
      
      # tryCatch(
      #   {
      
      if(freq<0.04|freq>0.96){
       p_value1[temp] <- 1
       p_value2[temp] <- 1
        # score_result[temp,] <- 0
        # infor_result[temp,] <- 0
        # score_result2[temp,] <- 0
        # infor_result2[temp,] <- 0
      }else{
        
        x1 <- cbind(snpvalue,as.matrix(x.covar.mis2)) 
        x2 <- cbind(snpvalue2,as.matrix(x.covar.complete))
        model.standard <- glm(y.pheno.mis2[,1]~x1,family = binomial(link='logit'))
        
        p_value1[temp] <- summary(model.standard)$coefficients[2,4]
        
        poly.model <- multinom(subtypes~x2)
        if(poly.model$convergence==0){
          tryCatch({
            poly.model.coef <- coef(poly.model)
            M <- nrow(poly.model.coef)
            p.covariate <- ncol(poly.model.coef)
            snp.cov <- vcov(poly.model)[2+p.covariate*(0:(M-1)),2+p.covariate*(0:(M-1))]
            snp.coef <- poly.model.coef[,2]
            p_value2[temp] <- DisplaySecondStageTestResult(snp.coef,snp.cov)[35]  
          },
          error = function(e){
            p_value2[temp] <- 1
          }
            
          )
          
          
        }else{
          p_value2[temp] = 1
        }
        
        
        

        
      }
      
      
      
      
    }
    
    
    if(i==true.end){
      break
    }
    
    
    
    
  }
  close(con)
  result <- list(snpid_result,p_value1,p_value2,freq.all)
  
  return(result)
}
stopImplicitCluster()

# score_result <- matrix(0.1,file.num,num.of.tumor+1)
# infor_result <- matrix(0.1,file.num,(num.of.tumor+1)^2)
# 
# 
# score_result2 <- matrix(0.1,file.num,num.of.tumor-1)
# infor_result2 <- matrix(0.1,file.num,(num.of.tumor-1)^2)
snpid_result <- rep("c",file.num)
p_value1 <- rep(1,file.num)
p_value2 <- rep(1,file.num)
freq.all <- rep(0,file.num)

total <- 0
for(i in 1:inner.size){
  result.temp <- result.list[[i]]
  temp <- length(result.temp[[1]])
  snpid_result[total+(1:temp)] <- result.temp[[1]]
  p_value1[total+(1:temp)] <- result.temp[[2]]
  p_value2[total+(1:temp)] <- as.numeric(result.temp[[3]])
  freq.all[total+(1:temp)] <- result.temp[[4]]
  # score_result2[total+(1:temp),] <- result.temp[[5]]
  # infor_result2[total+(1:temp),] <- result.temp[[6]]
  # 
  total <- temp+total
}

result <- list(snpid_reuslt=snpid_result,p_value1=p_value1,p_value2=p_value2,freq.all=freq.all)

save(result,file=paste0("./poland/result/whole_genome/standard_analysis",i1,"_",i2))
