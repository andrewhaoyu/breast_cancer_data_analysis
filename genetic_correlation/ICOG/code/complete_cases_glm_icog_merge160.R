setwd("/data/zhangh24/breast_cancer_data_analysis/")
total <- 0
i1 = 160
for(i2 in 1:1000){
  print(i2)
  load(paste0("./genetic_correlation/ICOG/result/complete_glm",i1,"_",i2))
  temp <- length(result[[1]])
  total <- total+temp
}
rs_id <- rep("c",total)
score_result <- matrix(0,total,5)
infor_result <- matrix(0,total,25)
freq <- rep(0,total)
total <- 0
for(i2 in 1:1000){
  print(i2)
  load(paste0("./genetic_correlation/ICOG/result/complete_glm",i1,"_",i2))
  temp <- length(result[[1]])
  rs_id[total+(1:temp)] <- result[[1]]
  score_result[total+(1:temp),] <- result[[2]]
  infor_result[total+(1:temp),] <- result[[3]]
  freq[total+(1:temp)] <- result[[4]]
  total <- total+temp
}
result <- list(rs_id,score_result,infor_result,freq)
idx <- which(is.na(score_result)==T)
idx
score_result[399,]
infor_result[399,]
score_result[399,] <- rep(0,5)
infor_result[399,] <- as.vector(diag(1,5))
freq[399] <- 0
result <- list(rs_id,score_result,infor_result,freq)
save(result,file=paste0("./genetic_correlation/ICOG/result/complete_glm",i1))
