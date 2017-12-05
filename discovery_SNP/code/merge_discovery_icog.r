
discovery_snp <- read.csv("/spin1/users/zhangh24/breast_cancer_data_analysis/data/discovery_snp_summary.csv",header=T)


setwd("/spin1/users/zhangh24/breast_cancer_data_analysis/")
snp.icogs.extract.id <- as.character(discovery_snp$SNP.ICOGS)
snp.onco.extract.id <- as.character(discovery_snp$SNP.ONCO)

setwd("/spin1/users/zhangh24/breast_cancer_data_analysis/")

n.raw <- 109713
snpvalue <- rep(0,n.raw)
subject.file <- "/gpfs/gsfs4/users/NC_BW/icogs_onco/genotype/imputed2/icogs_order.txt.gz"
library(data.table)
Icog.order <- read.table(gzfile(subject.file))
setwd("/spin1/users/zhangh24/breast_cancer_data_analysis/")
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

###one snp missing in ICOG
extract.num <- length(snp.icogs.extract.id)-1
snpid.result <- rep("c",extract.num)
n.sub <- 72411
snpvalue.result <- matrix(0,n.sub,extract.num)



total <- 0

for(i in 1:564){
  
  print(i)  
  geno.file <- paste0("/spin1/users/zhangh24/breast_cancer_data_analysis/discovery_SNP/result/discovery_icog",i,".txt"
  )
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




idx.match <- match(snp.icogs.extract.id,snpid.result)
snpid.result <- snpid.result[idx.match]
all.equal(snpid.result,snp.icogs.extract.id)
snpvalue.result <- snpvalue.result[,idx.match]
extract.result <- list(snpid.result,snpvalue.result)
colnames(snpvalue.result) <- snpid.result

write.csv(snpvalue.result,file="/spin1/users/zhangh24/breast_cancer_data_analysis/discovery_SNP/result/discovery_icog_data.csv",row.names = F,quote=F)






