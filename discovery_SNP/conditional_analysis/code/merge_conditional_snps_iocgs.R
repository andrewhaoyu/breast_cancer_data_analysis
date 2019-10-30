args = commandArgs(trailingOnly = T)
i1 <- as.numeric(args[[1]])


load("/data/zhangh24/breast_cancer_data_analysis/discovery_SNP/conditional_analysis/result/all.conditional.snps.Rdata")
library(data.table)
setwd("/data/zhangh24/breast_cancer_data_analysis/")

n.raw <- 109713
snpvalue <- rep(0,n.raw)
subject.file <- "/gpfs/gsfs4/users/NC_BW/icogs_onco/genotype/imputed2/icogs_order.txt.gz"
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


extract.num <- nrow(all.conditional.snps)
snpid.result <- rep("c",extract.num)
n.sub <- 72411
#library(bigmemory)






text <- system(paste0("cat | wc -l /data/zhangh24/breast_cancer_data_analysis/discovery_SNP/conditional_analysis/result/conditional_extract_icog",i1,".txt"),intern=T)
extract.num <- as.integer(gsub(paste0(" /data/zhangh24/breast_cancer_data_analysis/discovery_SNP/conditional_analysis/result/conditional_extract_icog",i1,".txt"),"",text))
#snpvalue.result <- matrix(0,n.sub,extract.num)

if(extract.num==0){
  conditional.snp.list.icog <- NULL
  save(conditional.snp.list.icog,file=paste0("/data/zhangh24/breast_cancer_data_analysis/discovery_SNP/conditional_analysis/result/conditional.snp.list.icog",i1,".Rdata"))
}else{
  snpvalue.result <- matrix(0,n.sub,extract.num)
  
  total <- 0
  
  for(i in i1:i1){
    
    print(i)  
    geno.file <- paste0("/data/zhangh24/breast_cancer_data_analysis/discovery_SNP/conditional_analysis/result/conditional_extract_icog",i,".txt"
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
  
  
  conditional.snp.list.icog <- list(snpid.result,snpvalue.result)
  save(conditional.snp.list.icog,file=paste0("/data/zhangh24/breast_cancer_data_analysis/discovery_SNP/conditional_analysis/result/conditional.snp.list.icog",i1,".Rdata"))
  
}


# ########## five snps only exists in ONCO array
# 
# extract.list.shared <- extract.list[1:total,]
# 
# save(extract.list.shared,file="/data/zhangh24/breast_cancer_data_analysis/whole_genome_age/ICOG/ERPRHER2GRADE_casecase/result/extract_list_shared.Rdata")
# extract.list.onco.only <- extract.list[(total+1):nrow(extract.list),]
# save(extract.list.onco.only,file="/data/zhangh24/breast_cancer_data_analysis/whole_genome_age/ICOG/ERPRHER2GRADE_casecase/result/extract_list_onco_only.Rdata")
# 
# 
# ##########
# idx.match <- match(extract.list.shared$SNP.ICOGS,snpid.result)
# snpid.result <- snpid.result[idx.match]
# all.equal(snpid.result,extract.list.shared$SNP.ICOGS)
# snpvalue.result <- snpvalue.result[,idx.match]
# extract.result <- list(snpid.result,snpvalue.result)
# colnames(snpvalue.result) <- snpid.result
# 
# 
# 
# extract.result <- list(snpid.result,snpvalue.result)
# save(extract.result,file="/data/zhangh24/breast_cancer_data_analysis/whole_genome_age/ICOG/ERPRHER2GRADE_casecase/result/extract_result_shared.Rdata")
# 
# 



















