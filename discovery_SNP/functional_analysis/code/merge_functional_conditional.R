
##########extract list was generated after we filter out all the SNPs with 1M around the four BCAC discovery SNPs###


#new_filter[,2] <- as.numeric(gsub(",","",new_filter[,2]))


setwd("/spin1/users/zhangh24/breast_cancer_data_analysis/")
load(paste0("./discovery_SNP/functional_analysis/result/functional_snp_conditional.Rdata"))
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


extract.num <- nrow(functional_snp_conditional)
# snpid.result <- rep("c",extract.num)
# n.sub <- 72411
# snpvalue.result <- matrix(0,n.sub,extract.num)
# 
# 
# 
# total <- 0
# #temp1 = rep(0,564)
# for(i in 1:564){
#   
#   print(i)  
#   geno.file <- paste0("/spin1/users/zhangh24/breast_cancer_data_analysis/discovery_SNP/functional_analysis/result/ICOG/funcitonal_conditional_icog",i,".txt"
#   )
#   temp.out <- system(paste0('wc -l ',geno.file),intern=T)
#   temp.out <- gsub(paste0("/spin1/users/zhangh24/breast_cancer_data_analysis/discovery_SNP/functional_analysis/result/ICOG/funcitonal_conditional_icog",i,".txt"),"",temp.out)
#   num <- as.numeric(temp.out)
#   #temp1[i] <- num
#  # }
#   
#   if(num!=0){
#     
#     
#     con <- file(geno.file)
#     temp <- 0
#     open(con)
#     for(i in 1:num){
#       
#       oneLine <- readLines(con,n=1)
#       myVector <- strsplit(oneLine," ")
#       snpid <- as.character(myVector[[1]][3])
#       
#       
#       temp <- temp+1
#       print(temp)
#       snpid.result[temp+total] <- snpid
#       #snpvalue <- rep(0,n)
#       
#       
#       snppro <- as.numeric(unlist(myVector)[7:length(myVector[[1]])])
#       
#       snpvalue <- convert(snppro,n.raw)
#       snpvalue <- snpvalue[idx.fil][idx.match]
#       snpvalue.result[,temp+total] <- snpvalue
#       
#       
#     }
#     close(con)
#     
#     total <- total+num
#     
#     
#   }
#   
#   # if(is.null(result[[1]])==0){
#   #   temp <- length(result[[1]])
#   #   snpid.result[total+(1:temp)] <- result[[1]]
#   #   snpvalue.result[,total+(1:temp)] <- result[[2]]
#   #   total <- temp+total
#   # }
# }
# snpid.result <- snpid.result[1:total]
# snpvalue.result <- snpvalue.result[,1:total]
# #load("./whole_genome_age/ICOG/ERPRHER2GRADE_fixed_baseline/result/Julie_snp_name_match.Rdata")
# 
# 
# #idx.remove <- which(is.na(snpid.result))
# #snpid.result <- snpid.result[-idx.remove]
# #snpvalue.result <- snpvalue.result[,-idx.remove]
# 
# idx.match <- match(functional_snp_conditional$SNP.ICOGS,snpid.result)
# 
# 
# 
# snpid.result <- snpid.result[idx.match]
# 
# # idx.temp <- which(is.na(snpid.result))
# # head(snpvalue.result[,idx.temp[1]])
# #functional_snp_conditional[idx.temp,c(13,14,16)]
# 
# all.equal(snpid.result,functional_snp_conditional$SNP.ICOGS)
# snpvalue.result <- snpvalue.result[,idx.match]
# extract.result <- list(snpid.result,snpvalue.result)
# colnames(snpvalue.result) <- snpid.result
# snpvalue.result.icog <- snpvalue.result
# save(snpvalue.result.icog,file="./discovery_SNP/functional_analysis/result/ICOG/functional_conditional_icog_snpvalue.Rdata")




########merge CIMBA onco
n.raw <- 142273
snpvalue <- rep(0,n.raw) 
subject.file <- "/gpfs/gsfs4/users/NC_BW/icogs_onco/genotype/imputed2/onco_order.txt.gz"
onco.order <- read.table(gzfile(subject.file))
setwd("/spin1/users/zhangh24/breast_cancer_data_analysis/")
library(data.table)
data2 <- fread("./data/Onco_euro_v10_10232017.csv",header=T)
data2 <- as.data.frame(data2)

Onc_ID <- data2$Onc_ID


idx.fil <- onco.order[,1]%in%Onc_ID
n <- length(idx.fil)
snpvalue <- rep(0,n)
idx.match <- match(Onc_ID,onco.order[idx.fil,1])
extract.num <- nrow(functional_snp_conditional)
library(bc2)

#extract.num <- nrow(extract.list)
snpid.result <- rep("c",extract.num)
n.sub <- nrow(data2)
snpvalue.result <- matrix(0,n.sub,extract.num)


total <- 0

for(i in 1:567){
  
  print(i)  
  geno.file <- paste0("/spin1/users/zhangh24/breast_cancer_data_analysis/discovery_SNP/functional_analysis/result/ONCO/functional_conditional_onco",i,".txt")
  temp.out <- system(paste0('wc -l ',geno.file),intern=T)
  temp.out <- gsub(paste0("/spin1/users/zhangh24/breast_cancer_data_analysis/discovery_SNP/functional_analysis/result/ONCO/functional_conditional_onco",i,".txt"),"",temp.out)
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
#load(paste0("./discovery_SNP/CIMBA_conditional_analysis/result/CIMBA_snp_name_match.Rdata"))
#load("./whole_genome_age/ICOG/ERPRHER2GRADE_fixed_baseline/result/Julie_snp_name_match.Rdata")



idx.match <- match(functional_snp_conditional$SNP.ONCO,snpid.result)
snpid.result <- snpid.result[idx.match]
all.equal(snpid.result,functional_snp_conditional$SNP.ONCO)
snpvalue.result <- snpvalue.result[,idx.match]
extract.result <- list(snpid.result,snpvalue.result)
colnames(snpvalue.result) <- snpid.result
snpvalue.result.onco <- snpvalue.result
save(snpvalue.result.onco,file="./discovery_SNP/functional_analysis/result/ONCO/functional_conditional_onco_snpvalue.Rdata")













