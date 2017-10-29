
##########extract list was generated after we filter out all the SNPs with 1M around the known SNPs region###
##########all the SNPs with pvalue <= 5E-06 was token out

load("/spin1/users/zhangh24/breast_cancer_data_analysis/whole_genome/ICOG/ERPRHER2_fixed/result/extract_list.Rdata")

setwd("/spin1/users/zhangh24/breast_cancer_data_analysis/")

n.raw <- 109713
snpvalue <- rep(0,n.raw)
subject.file <- "/gpfs/gsfs4/users/NC_BW/icogs_onco/genotype/imputed2/icogs_order.txt.gz"
Icog.order <- read.table(gzfile(subject.file))
pheno.file <- "./data/pheno.Icog"
load(pheno.file)
n.sub = nrow(pheno)
y.pheno.mis1 <- cbind(pheno$Behaviour1,pheno$PR_status1,pheno$ER_status1,pheno$HER2_status1)
colnames(y.pheno.mis1) = c("Behaviour1","PR_status1",
                           "ER_status1","HER2_status1")
idx.fil <- Icog.order[,1]%in%pheno$SG_ID
idx.match <- match(pheno$SG_ID,Icog.order[idx.fil,1])
#Icog.order.match <- Icog.order[idx.fil,1][idx.match]
library(bc2)


extract.num <- nrow(extract.list)
snpid.result <- rep("c",extract.num)
n.sub <- 72411
snpvalue.result <- matrix(0,n.sub,extract.num)



total <- 0

for(i in 1:564){

print(i)  
 geno.file <- paste0("/spin1/users/zhangh24/breast_cancer_data_analysis/whole_genome/ICOG/ERPRHER2_fixed/result/ERPRHER2_extract",i,".txt")
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
      snpvalue <- rep(0,n)
      
      
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

extract.result <- list(snpid.result,snpvalue.result)
save(extract.result,file="/spin1/users/zhangh24/breast_cancer_data_analysis/whole_genome/ICOG/ERPRHER2_fixed/result/extract_result.Rdata")







