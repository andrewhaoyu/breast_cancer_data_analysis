setwd('/spin1/users/zhangh24/breast_cancer_data_analysis/')
filedir <- '/spin1/users/zhangh24/breast_cancer_data_analysis/simulation/type_one_error/result/'
files <- dir(filedir,pattern="simu_result",full.names=T)

total <- 0
for(i1 in 3001:6000){
  print(i1)
  file = paste0("/spin1/users/zhangh24/breast_cancer_data_analysis/simulation/type_one_error/result//simu_result",i1,".Rdata")
  if(file%in%files==T){
    load(paste0("./simulation/type_one_error/result/simu_result",i1,".Rdata"))
    temp1 = length(result.list[[1]][[1]])/3
    temp2 = length(result.list[[2]][[1]])/3
    total = total+ temp1 +temp2
       
    
  }
}



p_global_result <- matrix(0,total,3)
p_heter_result <- matrix(0,total,3)
p_mglobal_result <- matrix(0,total,3)
p_mheter_result <- matrix(0,total,3)
p_indi_result <- matrix(0,total,3)


total <- 0
for(i1 in 3001:6000){
  print(i1)
  file = paste0("/spin1/users/zhangh24/breast_cancer_data_analysis/simulation/type_one_error/result//simu_result",i1,".Rdata")
  if(file%in%files==T){
    load(paste0("./simulation/type_one_error/result/simu_result",i1,".Rdata"))
    temp1 = length(result.list[[1]][[1]])/3
    temp2 = length(result.list[[2]][[1]])/3
    
    
    if(temp1*temp2!=0){
      p_global_result[total+(1:(temp1+temp2)),] <- 
        rbind(matrix(result.list[[1]][[1]],temp1,3),
              matrix(result.list[[2]][[1]],temp2,3))
      p_heter_result[total+(1:(temp1+temp2)),] <- 
        rbind(matrix(result.list[[1]][[2]],temp1,3),
              matrix(result.list[[2]][[2]],temp2,3))
      p_mglobal_result[total+(1:(temp1+temp2)),] <- 
        rbind(matrix(result.list[[1]][[4]],temp1,3),
              matrix(result.list[[2]][[4]],temp2,3))
      p_mheter_result[total+(1:(temp1+temp2)),] <- 
        rbind(matrix(result.list[[1]][[5]],temp1,3),
              matrix(result.list[[2]][[5]],temp2,3))
      p_indi_result[total+(1:(temp1+temp2)),] <- 
        rbind(matrix(result.list[[1]][[3]],temp1,3),
              matrix(result.list[[2]][[3]],temp2,3))
      total = total+temp1+temp2
    }else if(temp1==0){
      p_global_result[total+(1:(temp1+temp2)),] <- 
        matrix(result.list[[2]][[1]],temp2,3)
      p_heter_result[total+(1:(temp1+temp2)),] <- 
        matrix(result.list[[2]][[2]],temp2,3)
      p_mglobal_result[total+(1:(temp1+temp2)),] <- 
        matrix(result.list[[2]][[4]],temp2,3)
      p_mheter_result[total+(1:(temp1+temp2)),] <- 
        matrix(result.list[[2]][[5]],temp2,3)
      p_indi_result[total+(1:(temp1+temp2)),] <- 
      matrix(result.list[[2]][[3]],temp2,3)
      total = total+temp1+temp2
    }else if(temp2 ==0){
      p_global_result[total+(1:(temp1+temp2)),] <- 
        matrix(result.list[[1]][[1]],temp1,3)
        
      p_heter_result[total+(1:(temp1+temp2)),] <- 
        matrix(result.list[[1]][[2]],temp1,3)
              
      p_mglobal_result[total+(1:(temp1+temp2)),] <- 
        matrix(result.list[[1]][[4]],temp1,3)
              
      p_mheter_result[total+(1:(temp1+temp2)),] <- 
        matrix(result.list[[1]][[5]],temp1,3)
              
      p_indi_result[total+(1:(temp1+temp2)),] <- 
        matrix(result.list[[1]][[3]],temp1,3)
              
      total = total+temp1+temp2
      }
    
  }
}


# total <- 0
# for(i in 1:s_times){
#   print(i)
#   load(paste0("./simulation/type_one_error/result/simu_result",i,".Rdata"))
#   temp <- 100
#   p_global_result[total+(1:temp),] <- matrix(result[[1]],temp,3)
#   p_heter_result[total+(1:temp),] <- matrix(result[[2]],temp,3)
#   p_mglobal_result[total+(1:temp),] <- matrix(result[[4]],temp,3)
#   p_mheter_result[total+(1:temp),] <- matrix(result[[5]],temp,3)
#   p_indi_result[total+(1:temp),] <- matrix(result[[3]],temp,3)
#   total = total+temp
# }

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

apply(p_global_result,2,function(x){CountTypeOne(x,10^-5)})
apply(p_heter_result,2,function(x){CountTypeOne(x,10^-5)})
apply(p_indi_result,2,function(x){CountTypeOne(x,10^-5)})
apply(p_mglobal_result,2,function(x){CountTypeOne(x,10^-5)})
apply(p_mheter_result,2,function(x){CountTypeOne(x,10^-5)})

apply(p_global_result,2,function(x){CountTypeOne(x,10^-6)})
apply(p_heter_result,2,function(x){CountTypeOne(x,10^-6)})
apply(p_indi_result,2,function(x){CountTypeOne(x,10^-6)})
apply(p_mglobal_result,2,function(x){CountTypeOne(x,10^-6)})
apply(p_mheter_result,2,function(x){CountTypeOne(x,10^-6)})

#CountTypeOne(p_global_result,10^-4)
