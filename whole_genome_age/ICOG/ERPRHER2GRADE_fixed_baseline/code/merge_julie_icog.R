
##########extract list was generated after we filter out all the SNPs with 1M around the known SNPs region###
##########all the SNPs with pvalue <= 5E-06 was token out
new_filter <- read.csv("/data/zhangh24/breast_cancer_data_analysis/whole_genome/ICOG/ERPRHER2_fixed/result/Filter_based_on_Montse.csv",header=T,stringsAsFactors = F)
new_filter[,2] <- as.numeric(gsub(",","",new_filter[,2]))


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

for(i in 1:564){
  
  print(i)  
  geno.file <- paste0("/data/zhangh24/breast_cancer_data_analysis/whole_genome_age/ICOG/ERPRHER2GRADE_fixed_baseline/result/Julie_Icog",i,".txt"
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
load("./whole_genome_age/ICOG/ERPRHER2GRADE_fixed_baseline/result/Julie_snp_name_match.Rdata")



idx.match <- match(Julie_snp$SNP.ICOGS,snpid.result)
snpid.result <- snpid.result[idx.match]
all.equal(snpid.result,Julie_snp$SNP.ICOGS)
snpvalue.result <- snpvalue.result[,idx.match]
extract.result <- list(snpid.result,snpvalue.result)
colnames(snpvalue.result) <- snpid.result

write.csv(snpvalue.result,file="/data/zhangh24/breast_cancer_data_analysis/whole_genome_age/ICOG/ERPRHER2GRADE_fixed_baseline/result/Julie_snp_icog.csv",row.names = F,quote=F)







