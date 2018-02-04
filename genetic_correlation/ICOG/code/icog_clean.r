setwd("/spin1/users/zhangh24/breast_cancer_data_analysis/")

total <- 0
sig <- c(1:564)




for(i1 in sig){
  print(i1)
  load(paste0("./genetic_correlation/ICOG/result/intrinsic_i1",i1))  
  total <- total+length(result[[1]])
}

sigma <- matrix(n,total,25)
logodds <- matrix(n,total,5)
snpid <- rep("c",total)
total <- 0
for(i1 in sig){
  print(i1)
  load(paste0("./genetic_correlation/ICOG/result/intrinsic_i1",i1))   
  temp <- length(result[[1]])
  snpid[total+(1:temp)] <- result[[1]]
  logodds[total+(1:temp),] <- result[[2]]
  sigma[total+(1:temp),] <- result[[3]]
  total <- total+ temp
}




ICOG.result <- data.frame(SNP.ICOGS=snpid,logodds,sigma)


ICOG.result.clean <- merge(shared.data,ICOG.result,by.x="SNP.ICOGS",
                           by.y = "SNP.ICOGS")

save(ICOG.result.clean,file= "./genetic_correlation/result/ICOG.result.clean.Rdata")

#load("./genetic_correlation/result/hapmap3list.Rdata")




#dele  = c(56,160,281,291,292,299,311,350,351,416,421,422,435)