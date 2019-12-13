setwd('/data/zhangh24/breast_cancer_data_analysis/')
stimes <- 1000



final_result <- NULL
filedir <- "./simulation/EM_algorithm_evaluation/result/"
files <- dir(filedir,pattern="simu_result_high_",full.names=T)
for(i2 in 1:3){
  odds1 <- matrix(0,nrow=1000*10,7)
  sigma1 <- matrix(0,nrow=1000*10,7)
  total <- 0  
  for(i1 in 1:1000){
    print(i1)
    file <- paste0("./simulation/EM_algorithm_evaluation/result//simu_result_high_",i1,"_",i2,".Rdata")
   if(file%in%files){
      load(file)
      temp <- nrow(result[[1]])
      odds1[total+(1:temp),] <- result[[1]]
      sigma1[total+(1:temp),] <- result[[2]]
      total <- total+temp  
    }
  }
  odds1 <- odds1[1:total,]
  sigma1 <- sigma1[1:total,]
  result_temp <- round(apply(odds1,2,mean)-  c(0.08, 0.08, 0.05, 0.05, 0.05,0.05,0.05),3)
  final_result <- rbind(final_result,result_temp)
}
colnames(final_result) <- c("baseline","ER","PR","HER2",
                            "grade","T1","T2")
write.csv(final_result,file = "./simulation/EM_algorithm_evaluation/result/bias_result_high.csv")
