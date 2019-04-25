#Goal: prepare the data from 1000 genome for the simulations
#data preprocess
#1. download the vcf data from 1000 genome
#2. use plink to do the LD pruning with 1000 genome 0.05 r2 and 1Mb

result <- matrix("c",22,1)
for(i in 1:22){
  result[i,1] <- paste0("/spin1/users/zhangh24/plink --vcf ALL.chr",i,".phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz --indep-pairwise 1000'kb' 1 0.1 --out /spin1/users/zhangh24/KG.vcf/prunded_result/chr_",i)
}
write.table(result,
  file="/spin1/users/zhangh24/KG.vcf/LD_pruned.sh",
  row.names = F,
  col.names = F,
  quote=F)

#####calculate the MAF for all of the SNPs
result <- matrix("c",22,1)
for(i in 1:22){
  result[i,1] <- paste0("/spin1/users/zhangh24/plink --vcf ALL.chr",i,".phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz --freq --out /spin1/users/zhangh24/KG.vcf/MAF_result/chr_",i)
}
write.table(result,
            file="/spin1/users/zhangh24/KG.vcf/MAF_cal.sh",
            row.names = F,
            col.names = F,
            quote=F)
#####merge all the maf files into one
##### head -1 chr_1.frq > all.freq 
##### tail -n +2 -q chr_*.frq >> all.freq


#####read in the LD pruned SNPs
num.snp = rep(0,22)
for(i in 1:22){
  geno.file = paste0("/spin1/users/zhangh24/KG.vcf/prunded_result/chr_",i,".prune.in")
  temp = system(paste0("wc -l ",geno.file),intern=T)
  temp = as.numeric(gsub(geno.file,"",temp))
  num.snp[i] = temp
}
snp.pruned <- data.frame(rep("c",sum(num.snp)),stringsAsFactors = F)
library(data.table)
total <- 0
for(i in 1:22){
  snp.temp <- as.data.frame(fread(paste0("/spin1/users/zhangh24/KG.vcf/prunded_result/chr_",i,".prune.in"),header=F))
  snp.pruned[(total+1):(total+num.snp[i]),1] <- snp.temp
  total <- total+num.snp[i]
  
}


colnames(snp.pruned) <- "rs_id"

all.snp <- as.data.frame(fread("/spin1/users/zhangh24/KG.vcf/MAF_result/all.freq",header=T))
dim(all.snp)
pruned.snp.infor <- merge(snp.pruned,all.snp,
                          by.x = "rs_id",
                          by.y = "SNP")
library(dplyr)
pruned.snp.clean= pruned.snp.infor %>% 
  filter(MAF>=0.05)
