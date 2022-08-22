#merge the extracted icogs_results into one file

#load your extract list
#in your case, it's data_clean_select_sig
load("/data/zhangh24/breast_cancer_data_analysis/whole_genome_age/ICOG/ERPRHER2GRADE_fixed_baseline/result/extract_list.Rdata")
#change the directory to your directory
setwd("/data/zhangh24/breast_cancer_data_analysis/")


subject.file <- "/gpfs/gsfs4/users/NC_BW/icogs_onco/genotype/imputed2/icogs_order.txt.gz"
Icog.order <- read.table(gzfile(subject.file))
n.raw <- nrow(Icog.order)
snpvalue <- rep(0,n.raw)
setwd("/data/zhangh24/breast_cancer_data_analysis/")
data1 <- fread("./data/iCOGS_euro_v10_10232017.csv",header=T)
data1 <- as.data.frame(data1)
y.pheno.mis1 <- cbind(data1$Behaviour1,data1$ER_status1,data1$PR_status1,data1$HER2_status1,data1$Grade1)
colnames(y.pheno.mis1) = c("Behavior","ER","PR","HER2","Grade")
#x.test.all.mis1 <- data1[,c(27:206)]
SG_ID <- data1$SG_ID
idx.fil <- Icog.order[,1]%in%SG_ID
idx.match <- match(SG_ID,Icog.order[idx.fil,1])
#Icog.order.match <- Icog.order[idx.fil,1][idx.match]
#library(bc2)
extract.num <- nrow(extract.list)
snpid.result <- rep("c",extract.num)
n.sub <- nrow(data1)
snpvalue.result <- matrix(0,n.sub,extract.num)


#run through all the files to merge the extracted snps
total <- 0

for(i in 1:564){
  #change the output to your folder which was extracted in the last step
  print(i)  
  geno.file <- paste0("/data/zhangh24/breast_cancer_data_analysis/whole_genome_age/ICOG/ERPRHER2GRADE_fixed_baseline/result/ERPRHER2_extract",i,".txt"
)
  num <- as.numeric(system(paste0('cat ',geno.file,' | wc -l '),intern=T))
  #if num is 0, if means there is no extracted snps in the file[i]
  #we only need the file with num!=0
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
#it's a large matrix; each row is a subject; each column is a SNP with value ranging from 0 to 2
snpvalue.result <- snpvalue.result[,1:total]


##########match the extracted results to the original order
##########extract.list.shared in your case is data_clean_select_sig
idx.match <- match(extract.list.shared$SNP.ICOGS,snpid.result)
snpid.result <- snpid.result[idx.match]
all.equal(snpid.result,extract.list.shared$SNP.ICOGS)
snpvalue.result <- snpvalue.result[,idx.match]
extract.result <- list(snpid.result,snpvalue.result)
colnames(snpvalue.result) <- snpid.result
extract.result <- list(snpid.result,snpvalue.result)
#save the extracted results to your local folder
save(extract.result,file="/data/zhangh24/breast_cancer_data_analysis/whole_genome_age/ICOG/ERPRHER2GRADE_fixed_baseline/result/extract_result_shared.Rdata")







