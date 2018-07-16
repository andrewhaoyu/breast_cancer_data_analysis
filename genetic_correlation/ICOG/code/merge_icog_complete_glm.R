###########merge the intrinsic subtype results based on only complete data
setwd("/spin1/users/zhangh24/breast_cancer_data_analysis/")

total <- 0
sig <- c(1:564)
temp1 <- rep(0,564)



for(i1 in sig){
  print(i1)
  load(paste0("./genetic_correlation/ICOG/result/complete_glm",i1))  
  temp1[i1] <- length(result[[1]])
  total <- total+length(result[[1]])
}

sigma <- matrix(0,total,25)
logodds <- matrix(0,total,5)
snpid <- rep("c",total)
freq.icog <- rep(0,total)
total <- 0
for(i1 in sig){
  print(i1)
  load(paste0("./genetic_correlation/ICOG/result/complete_glm",i1)) 
  temp <- length(result[[1]])
  snpid[total+(1:temp)] <- result[[1]]
  logodds[total+(1:temp),] <- result[[2]]
  sigma[total+(1:temp),] <- result[[3]]
  freq.icog[total+(1:temp)] <- result[[4]]
  total <- total+ temp
}
snpid.temp <- snpid

load("./genetic_correlation/result/hapmap3list.Rdata")

ICOG.result <- data.frame(SNP.ICOGS=snpid,logodds,sigma,freq.icog)


ICOG.result.clean <- merge(shared.data,ICOG.result,by.x="SNP.ICOGS",
                           by.y = "SNP.ICOGS")

total <- nrow(ICOG.result.clean)


alleles.ICOG <- as.character(ICOG.result.clean$SNP.ICOGS)

alleles1 <- rep("c",total)
alleles2 <- rep("c",total)
alleles.split.icog <- strsplit(alleles.ICOG,split=":")

alleles.ONCO <- as.character(ICOG.result.clean$SNP.ONCO)
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

#alleles.data <- data.frame(alleles1,alleles2,alleles3,alleles4)


idx <- which(is.na(alleles1)&!is.na(alleles3))
alleles1[idx] <- alleles3[idx]
alleles2[idx] <- alleles4[idx]
ICOG.result.clean$A1 <- alleles1
ICOG.result.clean$A2 <- alleles2
#snpinfor <- data.frame(id,alleles1,alleles2)









# load(paste0("./genetic_correlation/ICOG/result/ICOG.result.Rdata"))
# #load(paste0("./genetic_correlation/ONCO/result/ONCO.result.transfrom.Rdata"))
# head(ICOG.result[[1]])
# head(ICOG.result.clean)
# ICOG.result.clean$A1 <- ICOG.result[[1]]$alleles1
# ICOG.result.clean$A2 <- ICOG.result[[1]]$alleles2

save(ICOG.result.clean,file= "./genetic_correlation/result/ICOG.result.clean.completeglm.Rdata")
#load("./genetic_correlation/result/hapmap3list.Rdata")




#dele  = c(56,160,281,291,292,299,311,350,351,416,421,422,435)