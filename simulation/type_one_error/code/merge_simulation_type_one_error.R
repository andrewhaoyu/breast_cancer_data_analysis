setwd("/spin1/users/zhangh24/breast_cancer_data_analysis")
s_times <- 1000
p_global_result <- matrix(0,s_times*100,3)
p_heter_result <- matrix(0,s_times*100,3)
p_mglobal_result <- matrix(0,s_times*100,3)
p_mheter_result <- matrix(0,s_times*100,3)
p_indi_result <- matrix(0,s_times*100,3)


total <- 0
for(i in 1:s_times){
  print(i)
  load(paste0("./simulation/type_one_error/result/simu_result",i,".Rdata"))
  temp <- 100
  p_global_result[total+(1:temp),] <- matrix(result[[1]],temp,3)
  p_heter_result[total+(1:temp),] <- matrix(result[[2]],temp,3)
  p_mglobal_result[total+(1:temp),] <- matrix(result[[4]],temp,3)
  p_mheter_result[total+(1:temp),] <- matrix(result[[5]],temp,3)
  p_indi_result[total+(1:temp),] <- matrix(result[[3]],temp,3)
  total = total+temp
}

CountTypeOne <- function(p,alpha){
  n <- length(p)
  idx <- which(p<=alpha)
  return(length(idx)/n)
  
}
apply(p_global_result,2,function(x){CountTypeOne(x,10^-4)})
apply(p_heter_result,2,function(x){CountTypeOne(x,10^-4)})
apply(p_indi_result,2,function(x){CountTypeOne(x,10^-4)})
apply(p_mglobal_result,2,function(x){CountTypeOne(x,10^-4)})
apply(p_mheter_result,2,function(x){CountTypeOne(x,10^-4)})
CountTypeOne(p_global_result,10^-4)
