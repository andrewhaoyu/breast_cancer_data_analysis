#-------------------------------------------------------------------
# Update Date: 11/20/2018
# Create Date: 11/20/2018
# Goal: prepare assoc dataset for plink to perform LD clumping
# Author: Haoyu Zhang
#-------------------------------------------------------------------
setwd('/data/zhangh24/breast_cancer_data_analysis/risk_prediction/')
load("./intrinsic_subtypes_whole_genome/ICOG/result/whole_genome_threeadd.rdata")

#load("./EB_whole_genome/result/whole_gonome.rdata")


library(dplyr)
n <- nrow(whole_genome)



#take out all the SNPs that are within +-500kb of the 313 SNPs
#load 313 Nasim SNPs
setwd('/data/zhangh24/breast_cancer_data_analysis/')
library(dplyr)
load("./data/Nasim_313SNPs_complete_information.Rdata")
remove.ind <- NULL
#remove all the SNPs +-500kb of the 313 Nasim SNPs
for(k in 1:nrow(snp.new)){
  #find all the SNPs +-500kb in it and remove them
  print(k)
  ldx <- which(whole_genome$CHR==snp.new$Chromosome[k]&
                 whole_genome$position>=snp.new$Positionb[k]-500000&
                 whole_genome$position<=snp.new$Positionb[k]+500000&whole_genome$var_name%in%snp.new$var_name==F)
  print(k)
  remove.ind = c(remove.ind,ldx)
}
whole_genome <- whole_genome[-remove.ind,]
#check whether the 313 SNPs are in the list
head(snp.new)
idx <- which(snp.new$var_name%in%
               whole_genome$var_name==F)
length(idx)



n <- nrow(whole_genome)

assoc <- whole_genome %>% mutate(p.min = p.min,
                        TEST = rep("ADD",n),
                        NMISS = rep(0,n),
                        OR = exp(stan_logodds),
                        STAT = rnorm(n)) %>%
                        select(CHR,SNP.ONCO,position,
                               effect_allele,
                               TEST,NMISS,
                               OR,STAT,p.min,
                               var_name)

colnames(assoc) <- c("CHR",
                     "SNP",
                     "BP",
                     "A1",
                     "TEST",
                     "NMISS",
                     "OR",
                     "STAT",
                     "P",
                     "var_name")
#load 313 Nasim SNPs
setwd('/data/zhangh24/breast_cancer_data_analysis/')
library(dplyr)
load("./data/Nasim_313SNPs_complete_information.Rdata")
idx <- which((snp.new$var_name%in%assoc$var_name)==T)
length(idx)
min(assoc$P)
#find the 313 SNPs and put the p.value as 0 to make sure they stay in the LD_clumping procedure
jdx <- which((assoc$var_name)%in%snp.new$var_name==T)
length(jdx)
assoc$P[jdx] = 0



assoc <- assoc %>%
  select(CHR,
         SNP,
         BP,
         A1,
         TEST,
         NMISS,
         OR,
         STAT,
         P)

write.table(assoc,file="/data/zhangh24/BCAC/impute_plink_onco/LD_assoc",col.names = T, row.names=F, quote = F)

#library(data.table)
#data <- fread("/data/zhangh24/BCAC/impute_plink_onco/LD_assoc")
#idx <- which(data$SNP=="rs56305903:52041531:T:C")
#data[idx,]
#idx <- which(whole_genome$SNP.ONCO=="rs56305903:52041531:T:C")
#whole_genome[idx,]