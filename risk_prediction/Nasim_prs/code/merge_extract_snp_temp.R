#merge extracted oncoarray data for 313 nasim snps
extract_snp <-  "rs1416885:100010765:C:G"

n.raw <- 142273
snpvalue <- rep(0,n.raw) 
subject.file <- "/gpfs/gsfs4/users/NC_BW/icogs_onco/genotype/imputed2/onco_order.txt.gz"
onco.order <- read.table(gzfile(subject.file))


library(data.table)
library(bc2, lib.loc ="/home/zhangh24/R/x86_64-pc-linux-gnu-library/3.6/")


extract.num <- 1
snpid.result <- rep("c",extract.num)
snpvalue.result <- matrix(0,n.raw,extract.num)


total <- 0

for(i in 1:1){
  
  print(i)  
  geno.file <- paste0("/data/zhangh24/breast_cancer_data_analysis/risk_prediction/Nasim_prs/result/extracted_snp_onco_temp.txt")
  num <- as.numeric(system(paste0('cat ',geno.file,' | wc -l '),intern=T))
  
  if(num!=0){
    
    
    con <- file(geno.file)
    temp <- 0
    open(con)
    for(i in 1:num){
      
      oneLine <- readLines(con,n=1)
      myVector <- strsplit(oneLine," ")
      snpid <- as.character(myVector[[1]][3])
      
      
      temp <- temp+1
      snpid.result[temp+total] <- snpid
      snpvalue <- rep(0,n.raw)
      
      
      snppro <- as.numeric(unlist(myVector)[7:length(myVector[[1]])])
      
      snpvalue <- convert(snppro,n.raw)
      #snpvalue <- snpvalue[idx.fil][idx.match]
      snpvalue.result[,temp+total] <- snpvalue
      
      
    }
    close(con)
    
    total <- total+num
    
    
  }
  
  
}

snpid.result <- snpid.result[1:total]
snpvalue.result <- snpvalue.result[,1:total]
idx.match <- match(extract_snp,snpid.result)
snpid.result <- snpid.result[idx.match]
all.equal(snpid.result,extract_snp)
snpvalue.result <- matrix(snpvalue.result,ncol=1)
colnames(snpvalue.result) <- snpid.result
onco.nasim.snp <- cbind(onco.order,snpvalue.result)
colnames(onco.nasim.snp)[1] <- "ID"
save(onco.nasim.snp, file = "/data/zhangh24/breast_cancer_data_analysis/risk_prediction/Nasim_prs/result/onco.nasim.snp_temp.rdata")




#merge extracted icogs data for 313 nasim snps
setwd("/data/zhangh24/breast_cancer_data_analysis/")
load(paste0("./risk_prediction/Nasim_prs/result/nasim_snp_odds.rdata"))
load(paste0("./risk_prediction/Nasim_prs/result/nasim_snp_id.rdata"))
#extract oncoarray data
library(dplyr)
extract_snp_temp <-  snp_id %>%
  select(SNP.ICOGS) 
extract.num <- nrow(extract_snp_temp)
extract_snp <- rep("c",extract.num)
for(i in 1:extract.num){
  extract_snp[i] <- as.character(extract_snp_temp[i,1])
}
n.raw <- 109713
snpvalue <- rep(0,n.raw)
subject.file <- "/gpfs/gsfs4/users/NC_BW/icogs_onco/genotype/imputed2/icogs_order.txt.gz"
Icog.order <- read.table(gzfile(subject.file))


library(data.table)


snpid.result <- rep("c",extract.num)
snpvalue.result <- matrix(0,n.raw,extract.num)


total <- 0

for(i in 1:564){
  
  print(i)  
  geno.file <- paste0("/data/zhangh24/breast_cancer_data_analysis/risk_prediction/Nasim_prs/result/extracted_snp_icog_313_",i,".txt")
  num <- as.numeric(system(paste0('cat ',geno.file,' | wc -l '),intern=T))
  
  if(num!=0){
    
    
    con <- file(geno.file)
    temp <- 0
    open(con)
    for(i in 1:num){
      
      oneLine <- readLines(con,n=1)
      myVector <- strsplit(oneLine," ")
      snpid <- as.character(myVector[[1]][3])
      
      
      temp <- temp+1
      snpid.result[temp+total] <- snpid
      snpvalue <- rep(0,n.raw)
      
      
      snppro <- as.numeric(unlist(myVector)[7:length(myVector[[1]])])
      
      snpvalue <- convert(snppro,n.raw)
      #snpvalue <- snpvalue[idx.fil][idx.match]
      snpvalue.result[,temp+total] <- snpvalue
      
      
    }
    close(con)
    
    total <- total+num
    
    
  }
  
  
}

snpid.result <- snpid.result[1:total]
snpvalue.result <- snpvalue.result[,1:total]
idx.match <- match(extract_snp,snpid.result)
snpid.result <- snpid.result[idx.match]
all.equal(snpid.result,extract_snp)
snpvalue.result <- snpvalue.result[,idx.match]
colnames(snpvalue.result) <- snpid.result
icog.nasim.snp <- cbind(Icog.order,snpvalue.result)
colnames(icog.nasim.snp)[1] <- "ID"
save(icog.nasim.snp, file = "/data/zhangh24/breast_cancer_data_analysis/risk_prediction/Nasim_prs/result/icog.nasim.snp")
