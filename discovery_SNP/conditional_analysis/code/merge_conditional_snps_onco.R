args = commandArgs(trailingOnly = T)
i1 <- as.numeric(args[[1]])
setwd("/spin1/users/zhangh24/breast_cancer_data_analysis/")


n.raw <- 142273
snpvalue <- rep(0,n.raw) 
subject.file <- "/gpfs/gsfs4/users/NC_BW/icogs_onco/genotype/imputed2/onco_order.txt.gz"
onco.order <- read.table(gzfile(subject.file))


library(data.table)
data2 <- fread("./data/Onco_euro_v10_10232017.csv",header=T)
data2 <- as.data.frame(data2)
y.pheno.mis2 <- cbind(data2$Behaviour1,data2$ER_status1,data2$PR_status1,data2$HER2_status1,data2$Grade1)
#y.pheno.mis2 <- cbind(data2$Behaviour1,data2$PR_status1,data2$ER_status1,data2$HER2_status1)
colnames(y.pheno.mis2) = c("Behaviour","ER",
                           "PR","HER2","Grade")



Onc_ID <- data2$Onc_ID



idx.fil <- onco.order[,1]%in%Onc_ID
n <- length(idx.fil)
snpvalue <- rep(0,n)
idx.match <- match(Onc_ID,onco.order[idx.fil,1])

load("/spin1/users/zhangh24/breast_cancer_data_analysis/discovery_SNP/conditional_analysis/result/all.conditional.snps.Rdata")

library(bc2)

extract.num <- nrow(all.conditional.snps)
snpid.result <- rep("c",extract.num)
n.sub <- nrow(data2)
library(bigmemory)

text <- system(paste0("cat | wc -l /spin1/users/zhangh24/breast_cancer_data_analysis/discovery_SNP/conditional_analysis/result/conditional_extract_onco",i1,".txt"),intern=T)
extract.num <- as.integer(gsub(paste0(" /spin1/users/zhangh24/breast_cancer_data_analysis/discovery_SNP/conditional_analysis/result/conditional_extract_onco",i1,".txt"),"",text))



if(extract.num==0){
  conditional.snp.list.onco <- NULL 
  save(conditional.snp.list.onco,file=paste0("/spin1/users/zhangh24/breast_cancer_data_analysis/discovery_SNP/conditional_analysis/result/conditional.snp.list.onco",i1,".Rdata"))
  
}else{
  snpvalue.result <- big.matrix(n.sub,extract.num,init=0)
  total <- 0  
  for(i in i1:i1){
    
    print(i)  
    geno.file <- paste0("/spin1/users/zhangh24/breast_cancer_data_analysis/discovery_SNP/conditional_analysis/result/conditional_extract_onco",i,".txt")
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
    
    gc()
  }
  
  snpid.result <- snpid.result[1:total]
  snpvalue.result <- snpvalue.result[,1:total]
  
  
  conditional.snp.list.onco <- list(snpid.result,snpvalue.result)
  save(conditional.snp.list.onco,file=paste0("/spin1/users/zhangh24/breast_cancer_data_analysis/discovery_SNP/conditional_analysis/result/conditional.snp.list.onco",i1,".Rdata"))
  
}




