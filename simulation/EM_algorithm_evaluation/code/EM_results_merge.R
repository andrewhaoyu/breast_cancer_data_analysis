setwd('/data/zhangh24/breast_cancer_data_analysis/')
stimes <- 100

odds1 <- matrix(0,nrow=100*1000,5)
sigma1 <- matrix(0,nrow=100*1000,5)
odds2 <- matrix(0,nrow=100*1000,5)
sigma2 <- matrix(0,nrow=100*1000,5)

total <- 0
for(i1 in 1:1000){
  print(i1)
  load(paste0("./simulation/EM_algorithm_evaluation/result/simu_result",i1,".Rdata"))
  temp <- stimes
  odds1[total+(1:temp),] <- result[[1]]
  sigma1[total+(1:temp),] <- result[[2]]
  odds2[total+(1:temp),] <- result[[3]]
  sigma2[total+(1:temp),] <- result[[4]]
  total <- total+temp
}

round(apply(odds1,2,mean),3)
(apply(odds1,2,sd))
sqrt(apply(sigma1,2,mean))

round(apply(odds2,2,mean),3)
sqrt(apply(odds2,2,var))
sqrt(apply(sigma2,2,mean))
