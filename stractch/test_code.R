#-------------------------------------------------------------------
# Update Date: 11/21/2018
# Create Date: 11/2012018
# Goal: test tranforming code
# Author: Haoyu Zhang
#-------------------------------------------------------------------
#---------------------------------------#---------------------------------------
library(data.table)


code <- paste0("zcat /data/zhangh24/BCAC/impute_onco/chr_test.gz | awk '{print $2\t$5}' > snp.file")
system(code)

system("zcat /data/zhangh24/BCAC/impute_onco/chr_test.gz | awk '{print $2\"\\t\"$5}' > snp.file")

data <- read.table("/data/zhangh24/BCAC/impute_onco/snp.file",header=F)
n <- nrow(data)
beta <- rnorm(n)
data <- cbind(data,beta)
data <- data[1:78,]
colnames(data) <- c("SNP","effect_allele","beta")
write.table(data[,1,drop=F],file = "/data/zhangh24/BCAC/impute_onco/snp.file",row.names=F,col.names=F,quote=F)
write.table(data,file = "/data/zhangh24/BCAC/impute_onco/prs.file",row.names=F,col.names=T,quote=F)


"/data/zhangh24/plink --score /data/zhangh24/BCAC/impute_onco/prs.file no-sum no-mean-imputation --dosage /data/zhangh24/BCAC/impute_onco/dosage_try noheader skip0=1 skip1=1 format=1 --fam /data/zhangh24/BCAC/impute_onco/onco_plink.fam"





dosage <- read.table(gzfile("/data/zhangh24/BCAC/impute_onco/dosage_try"),header=F)

snp <- as.matrix(dosage[,c(6:ncol(dosage))])
snp <- t(snp)
prs <- snp%*%beta[1:78]/(78*2)
#dosage <- t(dosage)
