#plug in ldsc result into R
load("./genetic_correlation/result/ldsc_result_meta_082919.rda")

ldsc_result[[1]][1,1] <- 0.4917
ldsc_result[[3]][1,1] <- 0.0719
ldsc_result[[1]][2,2] <- 0.6201
ldsc_result[[3]][2,2] <- 0.0564
ldsc_result[[1]][3,3] <- 0.6885
ldsc_result[[3]][3,3] <- 0.154
ldsc_result[[1]][4,4] <- 0.7398
ldsc_result[[3]][4,4] <- 0.0927
ldsc_result[[1]][5,5] <- 0.5971
ldsc_result[[3]][5,5] <- 0.0768
ldsc_result[[1]][6,6] <- 0.3094
ldsc_result[[3]][6,6] <- 0.0813
for(i in 1:6){
  for(j in 1:6){
    if(i!=j){
      ldsc_result[[1]][i,j] <- NA
      ldsc_result[[3]][i,j] <- NA
    }
  }
}
ldsc_result[[2]][2,1] <- 0.4557
ldsc_result[[4]][2,1] <- 0.0521
ldsc_result[[2]][3,1] <- 0.5572
ldsc_result[[4]][3,1] <- 0.1256
ldsc_result[[2]][4,1] <- 0.6027
ldsc_result[[4]][4,1] <- 0.0840
ldsc_result[[2]][5,1] <- 0.4041
ldsc_result[[4]][5,1] <- 0.0750
ldsc_result[[2]][6,1] <- 0.8401
ldsc_result[[4]][6,1] <- 0.154
ldsc_result[[2]][3,2] <- 0.5677
ldsc_result[[4]][3,2] <- 0.0700
ldsc_result[[2]][4,2] <- 0.7389
ldsc_result[[4]][4,2] <- 0.0490
ldsc_result[[2]][5,2] <- 0.7951
ldsc_result[[4]][5,2] <- 0.0461
ldsc_result[[2]][6,2] <- 0.3865
ldsc_result[[4]][6,2] <- 0.0853
ldsc_result[[2]][4,3] <- 0.3499
ldsc_result[[4]][4,3] <- 0.1043
ldsc_result[[2]][5,3] <- 0.5935
ldsc_result[[4]][5,3] <- 0.1073
ldsc_result[[2]][6,3] <- 0.8049
ldsc_result[[4]][6,3] <- 0.2400
ldsc_result[[2]][5,4] <- 0.6932
ldsc_result[[4]][5,4] <- 0.0722
ldsc_result[[2]][6,4] <- 0.3814
ldsc_result[[4]][6,4] <- 0.1567
ldsc_result[[2]][6,5] <- 0.3109
ldsc_result[[4]][6,5] <- 0.1236

for(i in 1:6){
  for(j in 1:6){
    if(j>i){
      ldsc_result[[2]][i,j] <- ldsc_result[[2]][j,i]
      ldsc_result[[4]][i,j] <- ldsc_result[[4]][j,i]
    }
  }
}
save(ldsc_result, file = "./genetic_correlation/result/ldsc_result_meta_091219.rda")
