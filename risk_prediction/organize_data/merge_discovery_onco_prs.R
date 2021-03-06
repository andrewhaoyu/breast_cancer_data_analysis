
discovery_snp <- read.csv("/data/zhangh24/breast_cancer_data_analysis/data/discovery_snp_summary_new.csv",header=T)


setwd("/data/zhangh24/breast_cancer_data_analysis/")
snp.icogs.extract.id <- as.character(discovery_snp$SNP.ICOGS)
snp.onco.extract.id <- as.character(discovery_snp$SNP.ONCO)

setwd("/data/zhangh24/breast_cancer_data_analysis/")


n.raw <- 142273
snpvalue <- rep(0,n.raw) 
subject.file <- "/gpfs/gsfs4/users/NC_BW/icogs_onco/genotype/imputed2/onco_order.txt.gz"
onco.order <- read.table(gzfile(subject.file))
library(data.table)
setwd("/data/zhangh24/breast_cancer_data_analysis/")
data2 <- fread("./data/PRS_subtype_Onco_euro_v10_08012018.csv",
               header=T)
data2 <- as.data.frame(data2)
y.pheno.mis2 <- cbind(data2$Behaviour1,data2$ER_status1,data2$PR_status1,data2$HER2_status1,data2$Grade1)
colnames(y.pheno.mis2) = c("Behavior","ER","PR","HER2","Grade")


Onc_ID <- data2$Onc_ID


idx.fil <- onco.order[,1]%in%Onc_ID
n <- length(idx.fil)
snpvalue <- rep(0,n)
idx.match <- match(Onc_ID,onco.order[idx.fil,1])

library(bc2, lib.loc ="/home/zhangh24/R/x86_64-pc-linux-gnu-library/3.6/")


extract.num <- length(snp.onco.extract.id)
snpid.result <- rep("c",extract.num)
n.sub <- nrow(data2)
snpvalue.result <- matrix(0,n.sub,extract.num)


total <- 0

for(i in 1:567){
  
  print(i)  
  geno.file <- paste0("/data/zhangh24/breast_cancer_data_analysis/discovery_SNP/result/discovery_onco",i,".txt")
  num <- as.numeric(system(paste0('cat ',geno.file,' | wc -l '),intern=T))
  
  if(num!=0){
    
    
    con <- file(geno.file)
    temp <- 0
    open(con)
    for(k in 1:num){
      
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



idx.match <- match(snp.onco.extract.id,snpid.result)
snpid.result <- snpid.result[idx.match]
all.equal(snpid.result,snp.onco.extract.id)
snpvalue.result <- snpvalue.result[,idx.match]
extract.result <- list(snpid.result,snpvalue.result)
colnames(snpvalue.result) <- snpid.result

write.csv(snpvalue.result,file="/data/zhangh24/breast_cancer_data_analysis/discovery_SNP/result/discovery_onco_data_prs.csv",row.names = F,quote=F)



data1 <- fread("./data/PRS_subtype_icgos_pheno_v10_euro.csv",header=T)
#data1 <- fread("./data/iCOGS_euro_v10_10232017.csv",header=T)
data1 <- as.data.frame(data1)

names1 <- colnames(data1)[27:203]
names2 <- colnames(data2)[27:204]
idx <- which(names2%in%names1==F)
all.equal(names2[c(1:(idx-1),(idx+1):178)],names1)


pc2 <- data2[5:14]
snpvalue2 <- cbind(data2[,c(1:(idx-1),(idx+1):178)+26],snpvalue.result)
all.equal(colnames(snpvalue2)[1:177],names1)
age <- data2[,205]
special.snp <- data2[,26+idx]
ID <- data2[,1]
sig_snp_onco_prs <- cbind(ID,y.pheno.mis2,
                          pc2,
                          snpvalue2,
                          age,
                          special.snp)
write.csv(sig_snp_onco_prs,file = "/data/zhangh24/breast_cancer_data_analysis/data/sig_snp_onco_prs.csv")






