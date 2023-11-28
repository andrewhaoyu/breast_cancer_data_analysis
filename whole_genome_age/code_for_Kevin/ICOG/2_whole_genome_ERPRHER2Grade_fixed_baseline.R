rm(list=ls())
#i1 represent the index of genotype file
#i1 ranges from 1 to 596
arg <- commandArgs(trailingOnly=T)
i1 <- as.numeric(arg[[1]])

print(i1)
library(R.utils)
library(data.table)
library(devtools)
library(withr)
library(gtools)
library(doParallel)
library(foreach)
#install R package
#bc2 is a development version of TOP package
#I used bc2 in my previous analyses
#the function of bc2 and TOP are almost the same
#TOP has more documentation
#to install bc2 or TOP, one needs to use install_github function
#you can specify the directory to your local directory
#with_libpaths(new = "/home/zhangh24/R/x86_64-pc-linux-gnu-library/4.2/", install_github('andrewhaoyu/bc2'))
library(bc2, 
        lib.loc ="/home/zhangh24/R/x86_64-pc-linux-gnu-library/4.2/")

setwd("/data/zhangh24/breast_cancer_data_analysis/")

#imputation file subject order
if(i1<=564){
  subject.file <- "/data/NC_BW/icogs_onco/genotype/imputed2/icogs_order.txt.gz"
  Icog.order <- read.table(gzfile(subject.file))
  
}else{
  subject.file <- "/data/NC_BW/icogs_onco/genotype/imputed2/icogs_order_23.txt.gz"
  Icog.order <- read.table(gzfile(subject.file))
  
}

setwd("/data/zhangh24/breast_cancer_data_analysis/")
#load the phenotypes data
data1 <- fread("./data/iCOGS_euro_v10_10232017.csv",header=T)
data1 <- as.data.frame(data1)
y.pheno.mis1 <- cbind(data1$Behaviour1,data1$ER_status1,data1$PR_status1,data1$HER2_status1,data1$Grade1)
colnames(y.pheno.mis1) = c("Behavior","ER","PR","HER2","Grade")
#x.test.all.mis1 <- data1[,c(27:206)]
SG_ID <- data1$SG_ID
#load the covariates for the model: PC1-10, age
x.covar.mis1 <- data1[,c(5:14,204)]
age <- data1[,204]
#find the people with missing ages
idx.incomplete <- which(age==888)
table(y.pheno.mis1[idx.incomplete,1])
idx.complete <- which(age!=888)
#remove people with missing age
y.pheno.mis1 <- y.pheno.mis1[idx.complete,]
x.covar.mis1 <- x.covar.mis1[idx.complete,]
SG_ID <- SG_ID[idx.complete]


#number of subject in the genotype file is n
n <- length(Icog.order[,1])
#creat a intial value for snpvalue
snpvalue <- rep(0,n)
#the phenotype data is a subset a the genotype data
#find the correponding subset
idx.fil <- Icog.order[,1]%in%SG_ID
#match the phenotype data with genotype data
idx.match <- match(SG_ID,Icog.order[idx.fil,1])
#idx.fil and idx.match will be used in later step for matching phenotype and genotype

#load the null hypothesis results for other covariates
#this component will be needed in later ScoreTest
load("./whole_genome_age/ICOG/ERPRHER2GRADE_fixed_baseline/result/score.test.support.icog.ERPRHER2Grade.Rdata")

#load all the imputed files
Filesdir <- "/data/NC_BW/icogs_onco/genotype/imputed2/icogs_imputed/"
Files <- dir(Filesdir,pattern="icogs_merged_b1_12.",full.names=T)
#order the imputed files
Files <- mixedsort(Files)
#specific one filegeno.file
geno.file <- Files[i1]

#count the number of variants in the file
num <- as.integer(system(paste0("zcat ",geno.file,"| wc -l"),intern=T))
#num = 10
#number of tumor characteristis is four
num.of.tumor <- ncol(y.pheno.mis1)-1
#number of subject in the phenotype files
n.sub <- nrow(y.pheno.mis1)
idx.control <- which(y.pheno.mis1[,1]==0)
#count the number of control in the data
n.control <- length(idx.control)

#get the three different z design matrix
z.design.list = GenerateZDesignCombination(y.pheno.mis1)
z.additive = z.design.list[[1]]
z.interaction = z.design.list[[2]]
z.saturated = z.design.list[[3]]
#number of second stage parameters
#if use additive model
n.second = ncol(z.additive)
#if use pair-wise interaction model
#n.second = ncol(z.interaction)
#if use saturated model
#n.second = ncol(z.saturated)



#parallel computing with foreach function
#the default of biowulf job allocation is two cores
#without parallel, we are only using 50% of the computing resources
#the job is running on two cores simultaneously
#parallel computing is faster, but sometimes also hard to debug
#it's also okay to just use a single for loop
#single for loop is easier to debug
#here I am splitting the jobs onto two cores

no.cores <- 2
inner.size <- 2
registerDoParallel(no.cores)
result.list <- foreach(job.i = 1:inner.size)%dopar%{
  print(job.i)
  #startend is a function in bc2 package
  #specific the total loop number, the number of inner jobs
  #startend will equally split the total loop
  #startend will return with the start and the end of the job line
  #job.i is the index of the inner jobs
  #for example, if num = 10, inner.size =2, job.i = 1, then start = 1, end = 5
  #for example, if num = 10, inner.size =2, job.i = 2, then start = 6, end = 10
  start.end <- startend(num,inner.size,job.i)
  start <- start.end[1]
  end <- start.end[2]
  inner.num <- end-start+1
  #score_matrix, each row is the score vector for a genetic marker
  score_result <- matrix(0,inner.num,n.second)
  #information matrix, each row is the as.vector(information matrix) for a genetic marker
  infor_result <- matrix(0,inner.num,(n.second)^2)
  #snpid information
  snpid_result <- rep("c",inner.num)
  #frequencies of the genetic marker
  freq.all <- rep(0,inner.num)
  temp <- 0
  #open the file
  con <- gzfile(geno.file)
  open(con)
  for(i in 1:num){
   #print the index every 500 SNPs
    #if(i%%500==0){
    print(i)
    #}
    #read one line of genetic file
    oneLine <- readLines(con,n=1)
    
    #the total number of SNPs are split into two sub-jobs
    #only start run the test after the start location
    if(i>=start){
      temp <- temp+1
      #the readLine result is a vector
      myVector <- strsplit(oneLine," ")
      #load the SNP ID
      snpid <- as.character(myVector[[1]][2])
      snpid_result[temp] <- snpid
      snpvalue <- rep(0,n)
      
      #load the imputed score for the genetic marker
      #3 * number of subjects length
      #every three columns are the probality for aa, Aa, AA for one subject
      snppro <- as.numeric(unlist(myVector)[6:length(myVector[[1]])])
      if(length(snppro)!=(3*n)){
        break
      }
      
      #calculate the expected genotype score of the subject. Value between 0 to 2.
      snpvalue <- convert(snppro,n)
      #match the genotype to the phenotype data
      snpvalue <- snpvalue[idx.fil][idx.match]
      #calculate the allele frequencies only use controls
      snpvalue.control <- snpvalue[idx.control]
      freq <- sum(snpvalue.control)/(2*n.control)
      freq.all[temp] <- freq
      #print(paste0("freq",freq))
      

      #only keep SNPs with allele frequency between 0.006 to 0.994    
      if(freq<0.006|freq>0.994){
            
        #if the SNP is too rare, just keep as score 0.
            score_result[temp,] <- 0
            infor_result[temp,] <- 0.1
          }else{
            
            #fit the ScoreTest
            #change second.stage.structure to second.stage.structure = pairwise.interaction for interaction model
            #change second.stage.structure to second.stage.structure = saturated for saturated
            score.test.icog<- ScoreTest(y=y.pheno.mis1,
                                        x=snpvalue,
                                        second.stage.structure="additive",
                                        score.test.support=score.test.support.icog.ERPRHER2Grade,
                                        missingTumorIndicator=888)
            
            #the first element is score
            score_result[temp,]  <- score.test.icog[[1]]
            #the second element is the efficient information matrix
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


#the output of foreach is saved as two list
#the attached code combine the two list as one


score_result <- matrix(0,num,n.second)
infor_result <- matrix(0,num,(n.second)^2)
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
#change the directory to your local directory
save(result,file=paste0("./whole_genome_age/ICOG/ERPRHER2GRADE_fixed_baseline/result/ERPRHER2Grade_fixed_baseline",i1))

