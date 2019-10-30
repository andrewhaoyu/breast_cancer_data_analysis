arg <- commandArgs(trailingOnly=T)
i1 <- as.numeric(arg[[1]])
print(i1)
z.design <- matrix(c(
  c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1),
  c(0,0,0,0,0,1,1,0,0,0,0,0,1,1,1,0,0,0,0,0,1,1,1),
  c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0),
  c(0,0,0,0,1,0,0,0,0,0,0,1,0,0,0,0,0,0,0,1,0,0,0),
  c(1,0,0,0,0,0,0,1,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0)
),ncol=5)

transfunction <- function(lododds,sigma){
  beta <- z.design%*%logodds
  beta.sigma <- z.design%*%sigma%*%t(z.design)
  p <- c(1,2,5,6,19)
  beta <- beta[p]
  beta.sigma <- beta.sigma[p,p]
  return(c(beta,as.vector(beta.sigma)))
}
setwd("/data/zhangh24/breast_cancer_data_analysis/")
##############transform the results into intrinsic subtypes results
##############ICOG.result.clean uses Luminal a as reference and other intrinsic subtypes versus luminal A as second stage parameters
load("./genetic_correlation/result/ICOG.result.clean.Rdata")
n <- nrow(ICOG.result.clean)
library(bc2)
size = 1000
start.end<- startend(n,size,i1)
start <- start.end[1]
end <- start.end[2]

result.sub <- matrix(0,end-start+1,30)
total <- 0
for(i in start:end){
print(i)
    logodds <- as.numeric(as.vector(ICOG.result.clean[i,11:15]))
  sigma <- matrix(as.numeric(ICOG.result.clean[i,16:40]),5,5)  
  total <- total+1
  result.sub[total,] <- transfunction(logodds,sigma)
}

save(result.sub,file=paste0("./genetic_correlation/ICOG/result/result.sub",i1,".Rdata"))



# result.all <- matrix(0,n,30)
# total <- 0
# for(i in 1:size){
#   print(i)
#   load(paste0("./genetic_correlation/ICOG/result/result.sub",i,".Rdata"))
#   temp <- nrow(result.sub)
# 
#     result.all[total+(1:temp),] <- result.sub
#   total <- total+temp
# }
# 
# 
# log.odds <- result.all[,1:5]
# sd.odds <-  sqrt(result.all[,c(6,12,18,24,30)])
# id <- ICOG.result.clean[,c(3,6,7)]
# 
# alleles.ICOG <- as.character(ICOG.result.clean$SNP.ICOGS)
# 
# alleles1 <- rep("c",n)
# alleles2 <- rep("c",n)
# alleles.split.icog <- strsplit(alleles.ICOG,split=":")
# 
# alleles.ONCO <- as.character(ICOG.result.clean$SNP.ONCO)
# alleles3 <- rep("c",n)
# alleles4 <- rep("c",n)
# alleles.split.onco <- strsplit(alleles.ONCO,split=":")
# 
# 
# for(i in 1:n){
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
# idx <- which(is.na(alleles1)&!is.na(alleles3))
# alleles1[idx] <- alleles3[idx]
# alleles2[idx] <- alleles4[idx]
# 
# snpinfor <- data.frame(id,alleles1,alleles2)
# 
# colnames(log.odds) <- c("Triple Negative",
#   "Luminial A",
#   "HER2 Enriched",
#   "Luminal B",
#   "Luminal B HER2Neg")
# colnames(sd.odds) <- c("Triple Negative",
#                        "Luminial A",
#                        "HER2 Enriched",
#                        "Luminal B",
#                        "Luminal B HER2Neg")
# 
# 
# 
# ICOG.result <- list(snpinfor,log.odds,sd.odds)
# save(ICOG.result,file=paste0("./genetic_correlation/ICOG/result/ICOG.result.Rdata"))
# 
# idx <- which(is.na(alleles1)&is.na(alleles2))
