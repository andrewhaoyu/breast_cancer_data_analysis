load("/spin1/users/zhangh24/breast_cancer_data_analysis/discovery_SNP/conditional_analysis/result/all.conditional.snps.Rdata")
n.snp <- nrow(all.conditional.snps)
p.value <- rep(0,n.snp)
total <- 0
for(i1 in 1:3000){
  print(i1)
  load(paste0("/spin1/users/zhangh24/breast_cancer_data_analysis/discovery_SNP/conditional_analysis/result/psub.6th",i1,".Rdata") )
  temp <- length(p.value.all)
  p.value[total+(1:temp)] <- p.value.all
  total <- total + temp
}

conditional.results.6th.all <- cbind(all.conditional.snps,p.value)



known.flag.all <- conditional.results.6th.all$known.flag
conditional.results.6th <- NULL
for(i in 1:207){
  print(i)
  idx <- which(known.flag.all==i)
  print(length(idx))
  p.value.temp <- p.value[idx]
  p.value.min <- min(p.value.temp)
  if(p.value.min<=1E-05){
    idx.min <- which.min(p.value.temp)
    conditional.results.6th <- rbind(conditional.results.6th,
                                     conditional.results.6th.all[idx,][idx.min,]               )
  }
  
}


library(tidyverse)
conditional.results.6th <- mutate(conditional.results.6th,
                                  cat.known.flag = ifelse(known.flag%in%1:178,1,
                                                          ifelse(known.flag%in%179:188,2,3))
)
table(conditional.results.6th$cat.known.flag)



save(conditional.results.6th,file="/spin1/users/zhangh24/breast_cancer_data_analysis/discovery_SNP/conditional_analysis/result/conditional.results.6th.Rdata")
load("/spin1/users/zhangh24/breast_cancer_data_analysis/discovery_SNP/conditional_analysis/result/conditional.results.6th.Rdata")
write.csv(conditional.results.6th,file="/spin1/users/zhangh24/breast_cancer_data_analysis/discovery_SNP/conditional_analysis/result/conditional.results.6th.csv",quote=F,row.names = F)



#conditional.results.first <- conditional.results.first[1:207,]
#conditional.results.first[,1] <- NA
#conditional.results.first[,2] <- NA
#conditional.results.first[,3] <- NA
