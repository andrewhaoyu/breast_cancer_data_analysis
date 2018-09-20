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

load("/spin1/users/zhangh24/breast_cancer_data_analysis/whole_genome_age/ICOG/ERPRHER2GRADE_casecase/result/extract_list.Rdata")
library(bc2)

extract.num <- nrow(extract.list)
snpid.result <- rep("c",extract.num)
n.sub <- nrow(data2)
snpvalue.result <- matrix(0,n.sub,extract.num)


total <- 0

for(i in 1:567){
  
  print(i)  
  geno.file <- paste0("/spin1/users/zhangh24/breast_cancer_data_analysis/whole_genome_age/ONCO/ERPRHER2GRADE_casecase/result/ERPRHER2_extract",i,".txt")
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
  
  
}

snpid.result <- snpid.result[1:total]
snpvalue.result <- snpvalue.result[,1:total]

extract.list.shared <- extract.list[1:total,]
idx.match <- match(extract.list.shared$SNP.ONCO,snpid.result)
snpid.result <- snpid.result[idx.match]
all.equal(snpid.result,extract.list.shared$SNP.ONCO)
snpvalue.result <- snpvalue.result[,idx.match]
extract.result <- list(snpid.result,snpvalue.result)
colnames(snpvalue.result) <- snpid.result



extract.result <- list(snpid.result,snpvalue.result)
save(extract.result,file="/spin1/users/zhangh24/breast_cancer_data_analysis/whole_genome_age/ONCO/ERPRHER2GRADE_casecase/result/extract_result_shared.Rdata")


