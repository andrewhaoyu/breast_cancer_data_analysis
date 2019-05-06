#Goal: merge summary level statistics
n.snp <- 0
for(i1 in 1:500){
  if(i1%%50==0){
    print(i1)  
  }
  load(paste0("./multi_ethnic/result/pruned_geno/beta_estimate_",i1,".Rdata"))
  temp <- nrow(beta_result[[1]])
  n.snp = n.snp +temp
}

beta_summary_train <- matrix(0,n.snp,9)
beta_summary_test <- matrix(0,n.snp,9)
beta_summary_vad <- matrix(0,n.snp,9)


colnames(beta_summary_train) <- c("beta_EUR","sd_EUR","p_EUR",
                                  "beta_AFR","sd_AFR","p_AFR", "beta_LAC","sd_LAC","p_LAC")
n.snp <- 0
for(i1 in 1:500){
  if(i1%%50==0){
    print(i1)  
  }
  load(paste0("./multi_ethnic/result/pruned_geno/beta_estimate_",i1,".Rdata"))
  temp <- nrow(beta_result[[1]])
  beta_summary_train[n.snp+(1:temp),] <- beta_result[[1]]
  beta_summary_test[n.snp+(1:temp),] <- beta_result[[2]]
  beta_summary_vad[n.snp+(1:temp),] <- beta_result[[3]]
  n.snp = n.snp +temp
}

beta_result <- list(beta_summary_train,
                    beta_summary_test,
                    beta_summary_vad)
save(beta_result,file = paste0("./multi_ethnic/result/pruned_geno/beta_all_",1,".Rdata"))



