setwd('/data/zhangh24/breast_cancer_data_analysis/')
stimes <- 100



final_result <- NULL
for(i2 in 1:3){
  odds1 <- matrix(0,nrow=100*100,5)
  sigma1 <- matrix(0,nrow=100*100,5)
  total <- 0  
  for(i1 in 1:100){
    print(i1)
    load(paste0("./simulation/EM_algorithm_evaluation/result/simu_result_high",i1,"_",i2,".Rdata"))
    temp <- nrow(result[[1]])
    odds1[total+(1:temp),] <- result[[1]]
    sigma1[total+(1:temp),] <- result[[2]]
    total <- total+temp
  }
  result_temp <- round(apply(odds1,2,mean),3)-  c(0.08, 0.08, 0.05, 0.05, 0.05)
  final_result <- rbind(final_result,result_temp)
}
colnames(final_result) <- c("baseline","ER","PR","HER2",
                            "grade")
write.csv(final_result,file = "./simulation/EM_algorithm_evaluation/result/bias_result.csv")
