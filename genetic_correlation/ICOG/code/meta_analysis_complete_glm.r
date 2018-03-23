args = commandArgs(trailingOnly = T)
i1 = as.numeric(args[[1]])
setwd("/spin1/users/zhangh24/breast_cancer_data_analysis/")
load(paste0("./genetic_correlation/ICOG/result/icog.onco.merge.completeglm.Rdata"))
library(bc2)
size =1000
n <- nrow(icog.onco.merge)
start.end<- startend(n,size,i1)
start <- start.end[1]
end <- start.end[2]

result.sub <- matrix(0,end-start+1,30)
total <- 0
for(i in start:end){
  print(i)
  logodds.icog <- as.numeric(as.vector(icog.onco.merge[i,11:15]))
  sigma.icog <- matrix(as.numeric(icog.onco.merge[i,16:40]),5,5)
  logodds.onco <- as.numeric(as.vector(icog.onco.merge[i,50:54]))
  sigma.onco <- matrix(as.numeric(icog.onco.merge[i,55:79]),5,5)
  temp.result <- LogoddsMetaAnalysis(logodds.icog,
                                     sigma.icog,
                                     logodds.onco,
                                     sigma.onco)
  total <- total+1
  result.sub[total,] <- c(temp.result[[1]],
                          as.vector(temp.result[[2]]))
}

save(result.sub,file=paste0("./genetic_correlation/ICOG/result/result.sub.meta.completeglm",i1,".Rdata"))




# setwd("/spin1/users/zhangh24/breast_cancer_data_analysis/")
# total <- nrow(icog.onco.merge)
# 
# result.all <- matrix(0,total,30)
# total <- 0
# for(i in 1:size){
#   print(i)
#   load(paste0("./genetic_correlation/ICOG/result/result.sub.meta",i,".Rdata"))
#   temp <- nrow(result.sub)
#   result.all[total+(1:temp),] <- result.sub
#   total <- temp+total
# }
# 
# 
# 
# 
# 
# log.odds <- result.all[,1:5]
# var.odds <-  result.all[,c(6:30)]
# id <- icog.onco.merge[,c(3,6,7)]
# freq.meta <- cbind(icog.onco.merge$freq.icog,icog.onco.merge$freq.onco)
# 
# alleles.ICOG <- as.character(icog.onco.merge$SNP.ICOGS.x)
# 
# alleles1 <- rep("c",total)
# alleles2 <- rep("c",total)
# alleles.split.icog <- strsplit(alleles.ICOG,split=":")
# 
# alleles.ONCO <- as.character(icog.onco.merge$SNP.ONCO.x)
# alleles3 <- rep("c",total)
# alleles4 <- rep("c",total)
# alleles.split.onco <- strsplit(alleles.ONCO,split=":")
# 
# 
# for(i in 1:total){
#   print(i)
#   alleles1[i] <- alleles.split.icog[[i]][3]
#   alleles2[i] <- alleles.split.icog[[i]][4]
#   alleles3[i] <- alleles.split.onco[[i]][3]
#   alleles4[i] <- alleles.split.onco[[i]][4]
# }
# 
# alleles.data <- data.frame(alleles1,alleles2,alleles3,alleles4)
# 
# 
# idx <- which(!is.na(alleles1)&is.na(alleles3))
# alleles3[idx] <- alleles1[idx]
# alleles4[idx] <- alleles2[idx] 
# 
# snpinfor <- data.frame(id,alleles3,alleles4)
# 
# colnames(log.odds) <- c("Triple Negative",
#                         "Luminial A",
#                         "HER2 Enriched",
#                         "Luminal B",
#                         "Luminal B HER2Neg")
# # colnames(sd.odds) <- c("Triple Negative",
# #                        "Luminial A",
# #                        "HER2 Enriched",
# #                        "Luminal B",
# #                        "Luminal B HER2Neg")
# # 
# 
# 
# meta.result <- list(snpinfor,log.odds,var.odds,freq.meta)
# colnames(meta.result[[4]]) <- c("freq.icog","freq.onco")
# save(meta.result,file=paste0("./genetic_correlation/ICOG/result/meta.result.Rdata"))
# 
# 
# all.result <- list(meta.result,ICOG.result,ONCO.result)
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# idx <- which(is.na(alleles3)&is.na(alleles4))
# length(idx)
# 
