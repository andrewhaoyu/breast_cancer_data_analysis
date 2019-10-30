
##########extract list was generated after we filter out all the SNPs with 1M around the known SNPs region###
##########all the SNPs with pvalue <= 5E-06 was token out
new_filter <- read.table("./discovery_SNP/CIMBA_conditional_analysis/result/extract_id_icog_CIMBA.txt",header=T,stringsAsFactors = F)
#new_filter[,2] <- as.numeric(gsub(",","",new_filter[,2]))


setwd("/data/zhangh24/breast_cancer_data_analysis/")

n.raw <- 109713
snpvalue <- rep(0,n.raw)
subject.file <- "/gpfs/gsfs4/users/NC_BW/icogs_onco/genotype/imputed2/icogs_order.txt.gz"
library(data.table)
Icog.order <- read.table(gzfile(subject.file))
setwd("/data/zhangh24/breast_cancer_data_analysis/")
data1 <- fread("./data/iCOGS_euro_v10_10232017.csv",header=T)
data1 <- as.data.frame(data1)
y.pheno.mis1 <- cbind(data1$Behaviour1,data1$ER_status1,data1$PR_status1,data1$HER2_status1,data1$Grade1)
colnames(y.pheno.mis1) = c("Behavior","ER","PR","HER2","Grade")
#x.test.all.mis1 <- data1[,c(27:206)]
SG_ID <- data1$SG_ID
x.covar.mis1 <- data1[,c(5:14,204)]
idx.fil <- Icog.order[,1]%in%SG_ID
idx.match <- match(SG_ID,Icog.order[idx.fil,1])
#Icog.order.match <- Icog.order[idx.fil,1][idx.match]
library(bc2)


extract.num <- nrow(new_filter)
snpid.result <- rep("c",extract.num)
n.sub <- 72411
snpvalue.result <- matrix(0,n.sub,extract.num)



total <- 0

for(i in 394:394){
  
  print(i)  
  geno.file <- paste0("/data/zhangh24/breast_cancer_data_analysis/discovery_SNP/CIMBA_conditional_analysis/result/ICOG/CIMBA_Icog",i,".txt"
  )
  temp.out <- system(paste0('wc -l ',geno.file),intern=T)
  temp.out <- gsub("/data/zhangh24/breast_cancer_data_analysis/discovery_SNP/CIMBA_conditional_analysis/result/ICOG/CIMBA_Icog394.txt","",temp.out)
  num <- as.numeric(temp.out)
  
  if(num!=0){
    
    
    con <- file(geno.file)
    temp <- 0
    open(con)
    for(i in 1:num){
      
      oneLine <- readLines(con,n=1)
      myVector <- strsplit(oneLine," ")
      snpid <- as.character(myVector[[1]][3])
      
      
      temp <- temp+1
      print(temp)
      snpid.result[temp+total] <- snpid
      #snpvalue <- rep(0,n)
      
      
      snppro <- as.numeric(unlist(myVector)[7:length(myVector[[1]])])
      
      snpvalue <- convert(snppro,n.raw)
      snpvalue <- snpvalue[idx.fil][idx.match]
      snpvalue.result[,temp+total] <- snpvalue
      
      
    }
    close(con)
    
    total <- total+num
    
    
  }
  
  # if(is.null(result[[1]])==0){
  #   temp <- length(result[[1]])
  #   snpid.result[total+(1:temp)] <- result[[1]]
  #   snpvalue.result[,total+(1:temp)] <- result[[2]]
  #   total <- temp+total
  # }
}
snpid.result <- snpid.result[1:total]
snpvalue.result <- snpvalue.result[,1:total]
#load("./whole_genome_age/ICOG/ERPRHER2GRADE_fixed_baseline/result/Julie_snp_name_match.Rdata")



idx.match <- match(CIMBA_snp$SNP.ICOGS,snpid.result)
snpid.result <- snpid.result[idx.match]
all.equal(snpid.result,CIMBA_snp$SNP.ICOGS)
snpvalue.result <- snpvalue.result[,idx.match]
extract.result <- list(snpid.result,snpvalue.result)
colnames(snpvalue.result) <- snpid.result
snpvalue.result.icog <- snpvalue.result
save(snpvalue.result.icog,file="./discovery_SNP/CIMBA_conditional_analysis/result/ICOG/CIMBA_icog_snpvalue.Rdata")




########merge CIMBA onco
n.raw <- 142273
snpvalue <- rep(0,n.raw) 
subject.file <- "/gpfs/gsfs4/users/NC_BW/icogs_onco/genotype/imputed2/onco_order.txt.gz"
onco.order <- read.table(gzfile(subject.file))
setwd("/data/zhangh24/breast_cancer_data_analysis/")
data2 <- fread("./data/Onco_euro_v10_10232017.csv",header=T)
data2 <- as.data.frame(data2)

Onc_ID <- data2$Onc_ID


idx.fil <- onco.order[,1]%in%Onc_ID
n <- length(idx.fil)
snpvalue <- rep(0,n)
idx.match <- match(Onc_ID,onco.order[idx.fil,1])

library(bc2)

#extract.num <- nrow(extract.list)
snpid.result <- rep("c",extract.num)
n.sub <- nrow(data2)
snpvalue.result <- matrix(0,n.sub,extract.num)


total <- 0

for(i in 397:397){
  
  print(i)  
  geno.file <- paste0("/data/zhangh24/breast_cancer_data_analysis/discovery_SNP/CIMBA_conditional_analysis/result/ONCO/CIMBA_onco397.txt")
  temp.out <- system(paste0('wc -l ',geno.file),intern=T)
  temp.out <- gsub("/data/zhangh24/breast_cancer_data_analysis/discovery_SNP/CIMBA_conditional_analysis/result/ONCO/CIMBA_onco397.txt","",temp.out)
  num <- as.numeric(temp.out)
  
  
  if(num!=0){
    
    
    con <- file(geno.file)
    temp <- 0
    open(con)
    for(k in 1:num){
      
      oneLine <- readLines(con,n=1)
      myVector <- strsplit(oneLine," ")
      snpid <- as.character(myVector[[1]][3])
      
      
      temp <- temp+1
      print(temp)
      snpid.result[temp+total] <- snpid
      snpvalue <- rep(0,n)
      
      
      snppro <- as.numeric(unlist(myVector)[7:length(myVector[[1]])])
      
      snpvalue <- convert(snppro,n.raw)
      snpvalue <- snpvalue[idx.fil][idx.match]
      snpvalue.result[,temp+total] <- snpvalue
      
      
    }
    close(con)
    
    total <- total+num
    
    
  }
  
  
}

snpid.result <- snpid.result[1:total]
snpvalue.result <- snpvalue.result[,1:total]
load(paste0("./discovery_SNP/CIMBA_conditional_analysis/result/CIMBA_snp_name_match.Rdata"))
#load("./whole_genome_age/ICOG/ERPRHER2GRADE_fixed_baseline/result/Julie_snp_name_match.Rdata")



idx.match <- match(CIMBA_snp$SNP.ONCO,snpid.result)
snpid.result <- snpid.result[idx.match]
all.equal(snpid.result,CIMBA_snp$SNP.ONCO)
snpvalue.result <- snpvalue.result[,idx.match]
extract.result <- list(snpid.result,snpvalue.result)
colnames(snpvalue.result) <- snpid.result
snpvalue.result.onco <- snpvalue.result
save(snpvalue.result.onco,file="./discovery_SNP/CIMBA_conditional_analysis/result/ONCO/CIMBA_onco_snpvalue.Rdata")













