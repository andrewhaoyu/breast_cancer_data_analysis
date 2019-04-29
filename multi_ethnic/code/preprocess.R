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
########read in phenotype information
phenotype <- as.data.frame(fread("/spin1/users/zhangh24/KG.vcf/integrated_call_samples_v2.20130502.ALL.ped"))



AFR <- c("GWD","MSL","ESN","GWJ","YRI","LWK","GWF","ASW","ACB","GWW")
AMR <- c("MXL","CLM","PEL","PUR")
EUR <- c("TSI","IBS","GBR","CEU","FIN")
SAS <- c("PJL","ITU","STU","GIH","BEB")
EAS <- c("JPT","CDX","CHB","CHS","KHV","CHD")


ID.EUR <- phenotype %>% 
  filter(Population%in%EUR) %>% 
  select(c("Family ID","Individual ID"))


ID.AFR <- phenotype %>% 
  filter(Population%in%AFR) %>% 
  select(c("Family ID","Individual ID"))


write.table(ID.EUR,file = "/spin1/users/zhangh24/KG.vcf/ID.EUR",row.names = F,
            col.names = F,quote=F)
write.table(ID.AFR,file = "/spin1/users/zhangh24/KG.vcf/ID.AFR",row.names = F,
            col.names = F,quote=F)


#####calculate the MAF for all of the SNPs for EUR and AFR
result <- matrix("c",44,1)
for(i in 1:22){
  result[i,1] <- paste0("/spin1/users/zhangh24/plink --vcf ALL.chr",i,".phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz --keep /spin1/users/zhangh24/KG.vcf/ID.EUR --freq --out /spin1/users/zhangh24/KG.vcf/MAF_result/EUR_chr_",i)
}

for(i in 1:22){
  result[i+22,1] <- paste0("/spin1/users/zhangh24/plink --vcf ALL.chr",i,".phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz --keep /spin1/users/zhangh24/KG.vcf/ID.AFR --freq --out /spin1/users/zhangh24/KG.vcf/MAF_result/AFR_chr_",i)
}
write.table(result,
            file="/spin1/users/zhangh24/KG.vcf/MAF_cal.sh",
            row.names = F,
            col.names = F,
            quote=F)
#####merge all the maf files into one
##### head -1 chr_1.frq > all.freq 
##### tail -n +2 -q chr_*.frq >> all.freq
##### head -1 EUR_chr_1.frq >> all_EUR.freq
##### tail -n +2 -q EUR_chr_*.frq >> all_EUR.freq
##### head -1 AFR_chr_1.frq >> all_AFR.freq
##### tail -n +2 -q AFR_chr_*.frq >> all_AFR.freq

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

all.snp.EUR <- as.data.frame(fread("/spin1/users/zhangh24/KG.vcf/MAF_result/all_EUR.freq",header=T))
colnames(all.snp.EUR)[5] <- "MAF.EUR"

all.snp.AFR <-  as.data.frame(fread("/spin1/users/zhangh24/KG.vcf/MAF_result/all_AFR.freq",header=T))
colnames(all.snp.AFR)[5] <- "MAF.AFR"
all.equal(all.snp.AFR$SNP,all.snp.EUR$SNP)

all.snp <- data.frame(all.snp.EUR,all.snp.AFR$MAF.AFR)

dim(all.snp)
pruned.snp.infor <- merge(snp.pruned,all.snp,
                          by.x = "rs_id",
                          by.y = "SNP")
library(dplyr)
pruned.snp.clean= pruned.snp.infor %>% 
  filter(MAF.EUR>=0.05&
        rs_id!=".")
save(pruned.snp.clean,file= "/spin1/users/zhangh24/KG.vcf/MAF_result/pruned_MAF.Rdata")


