load("/data/zhangh24/breast_cancer_data_analysis/discovery_SNP/conditional_analysis/result/all.conditional.snps.Rdata")
n.snp <- nrow(all.conditional.snps)
p.value <- rep(0,n.snp)
total <- 0
for(i1 in 1:3000){
  print(i1)
  load(paste0("/data/zhangh24/breast_cancer_data_analysis/discovery_SNP/conditional_analysis/result/psub",i1,".Rdata"))  
  temp <- length(p.value.all)
  p.value[total+(1:temp)] <- p.value.all
  total <- total + temp
}

conditional.results.first.all <- cbind(all.conditional.snps,p.value)



known.flag.all <- conditional.results.first.all$known.flag
conditional.results.first <- NULL
for(i in 1:207){
  print(i)
  idx <- which(known.flag.all==i)
  print(length(idx))
  p.value.temp <- p.value[idx]
  p.value.min <- min(p.value.temp)
  if(p.value.min<=1E-05){
    idx.min <- which.min(p.value.temp)
    conditional.results.first <- rbind(conditional.results.first,
                                       conditional.results.first.all[idx,][idx.min,]               )
  }
  
}


library(tidyverse)
conditional.results.first <- mutate(conditional.results.first,
  cat.known.flag = ifelse(known.flag%in%1:178,1,
  ifelse(known.flag%in%179:188,2,3))
  )
table(conditional.results.first$cat.known.flag)



save(conditional.results.first,file="/data/zhangh24/breast_cancer_data_analysis/discovery_SNP/conditional_analysis/result/conditional.results.first.Rdata")




write.csv(conditional.results.first,file="/data/zhangh24/breast_cancer_data_analysis/discovery_SNP/conditional_analysis/result/conditional.results.first.csv",quote=F,row.names = F)


#conditional.results.first <- conditional.results.first[1:207,]
#conditional.results.first[,1] <- NA
#conditional.results.first[,2] <- NA
#conditional.results.first[,3] <- NA
