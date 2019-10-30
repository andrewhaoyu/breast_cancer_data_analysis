setwd("/data/zhangh24/breast_cancer_data_analysis/")

load( "./genetic_correlation/ONCO/result/result.clean.Rdata")

total <- 0
size =1000
for(i in 1:size){
  print(i)
  load(paste0("./genetic_correlation/ONCO/result/result.sub",i,".Rdata"))
  temp <- nrow(result.sub)
  
#  result.all[total+(1:temp),] <- result.sub
  total <- total+temp
}



result.all <- matrix(0,total,30)
total <- 0
size =1000
for(i in 1:size){
  print(i)
  load(paste0("./genetic_correlation/ONCO/result/result.sub",i,".Rdata"))
  temp <- nrow(result.sub)
  
  result.all[total+(1:temp),] <- result.sub
  total <- total+temp
}



ONCO.result.transform <- ONCO.result.clean
ONCO.result.transform[,11:40] <- result.all

save(ONCO.result.transform,file=paste0("./genetic_correlation/ONCO/result/ONCO.result.transfrom.Rdata"))





log.odds <- result.all[,1:5]
var.odds <-  result.all[,c(6:30)]
id <- ONCO.result.clean[,c(3,6,7)]
freq.onco <- ONCO.result.clean[,41]

alleles.ICOG <- as.character(ONCO.result.clean$SNP.ICOGS)

alleles1 <- rep("c",total)
alleles2 <- rep("c",total)
alleles.split.icog <- strsplit(alleles.ICOG,split=":")

alleles.ONCO <- as.character(ONCO.result.clean$SNP.ONCO)
alleles3 <- rep("c",total)
alleles4 <- rep("c",total)
alleles.split.onco <- strsplit(alleles.ONCO,split=":")


for(i in 1:total){
  print(i)
  alleles1[i] <- alleles.split.icog[[i]][3]
  alleles2[i] <- alleles.split.icog[[i]][4]
  alleles3[i] <- alleles.split.onco[[i]][3]
  alleles4[i] <- alleles.split.onco[[i]][4]
}

alleles.data <- data.frame(alleles1,alleles2,alleles3,alleles4)


idx <- which(!is.na(alleles1)&is.na(alleles3))
alleles3[idx] <- alleles1[idx]
alleles4[idx] <- alleles2[idx] 

snpinfor <- data.frame(id,alleles3,alleles4)

colnames(log.odds) <- c("Triple Negative",
                        "Luminial A",
                        "HER2 Enriched",
                        "Luminal B",
                        "Luminal B HER2Neg")
# colnames(sd.odds) <- c("Triple Negative",
#                        "Luminial A",
#                        "HER2 Enriched",
#                        "Luminal B",
#                        "Luminal B HER2Neg")



ONCO.result <- list(snpinfor,log.odds,var.odds,freq.onco)
save(ONCO.result,file=paste0("./genetic_correlation/ONCO/result/ONCO.result.Rdata"))

idx <- which(is.na(alleles3)&is.na(alleles4))
length(idx)
