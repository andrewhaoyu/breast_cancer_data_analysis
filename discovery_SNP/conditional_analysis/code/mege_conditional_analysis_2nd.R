load("/spin1/users/zhangh24/breast_cancer_data_analysis/discovery_SNP/conditional_analysis/result/all.conditional.snps.Rdata")
n.snp <- nrow(all.conditional.snps)
p.value <- rep(0,n.snp)
total <- 0
for(i1 in 1:3000){
  print(i1)
 load(paste0("/spin1/users/zhangh24/breast_cancer_data_analysis/discovery_SNP/conditional_analysis/result/psub.2nd",i1,".Rdata") )
  temp <- length(p.value.all)
  p.value[total+(1:temp)] <- p.value.all
  total <- total + temp
}

conditional.results.2nd.all <- cbind(all.conditional.snps,p.value)



known.flag.all <- conditional.results.2nd.all$known.flag
conditional.results.2nd <- NULL
for(i in 1:207){
  print(i)
  idx <- which(known.flag.all==i)
  print(length(idx))
  p.value.temp <- p.value[idx]
  p.value.min <- min(p.value.temp)
  if(p.value.min<=1E-05){
    idx.min <- which.min(p.value.temp)
    conditional.results.2nd <- rbind(conditional.results.2nd,
                                       conditional.results.2nd.all[idx,][idx.min,]               )
  }
  
}


library(tidyverse)
conditional.results.2nd <- mutate(conditional.results.2nd,
                                    cat.known.flag = ifelse(known.flag%in%1:178,1,
                                                            ifelse(known.flag%in%179:188,2,3))
)
table(conditional.results.2nd$cat.known.flag)



save(conditional.results.2nd,file="/spin1/users/zhangh24/breast_cancer_data_analysis/discovery_SNP/conditional_analysis/result/conditional.results.2nd.Rdata")





#conditional.results.first <- conditional.results.first[1:207,]
#conditional.results.first[,1] <- NA
#conditional.results.first[,2] <- NA
#conditional.results.first[,3] <- NA
