setwd('/spin1/users/zhangh24/breast_cancer_data_analysis/')
filedir <- '/spin1/users/zhangh24/breast_cancer_data_analysis/simulation/power/result/'
files <- dir(filedir,pattern="simu_high",full.names=T)
total <- 0
for(i1 in 1:3200){
  print(i1)
  file = paste0("/spin1/users/zhangh24/breast_cancer_data_analysis/simulation/power/result//simu_high",i1,".Rdata")
  if(file%in%files==T){
    load(paste0("./simulation/power/result/simu_high",i1,".Rdata")) 
    total = total+ length(result.list[[1]][[1]])/6 + length(result.list[[2]][[1]])/6
    
  }
}

#total = total*2

p_global_result <- matrix(0,total,6)
p_mglobal_result <- matrix(0,total,6)
p_standard <- matrix(0,total,6)
p_global_complete <- matrix(0,total,6)
#p_poly <- matrix(0,total,9)

total <- 0

for(i1 in 1:3200){
  print(i1)
  file = paste0("/spin1/users/zhangh24/breast_cancer_data_analysis/simulation/power/result//simu_high",i1,".Rdata")
  if(file%in%files==T){
    load(paste0("./simulation/power/result/simu_high",i1,".Rdata")) 
    temp1 = length(result.list[[1]][[1]])/6
    temp2 = length(result.list[[2]][[1]])/6
    temp = temp1+temp2
    if(temp1==0){
      p_global_result[total+(1:temp2),] <- matrix(result.list[[2]][[1]],ncol=6)
      
      p_mglobal_result[total+(1:temp2),] <- matrix(result.list[[2]][[2]],ncol=6)
      p_standard[total+(1:temp2),] <- matrix(result.list[[2]][[3]],ncol=6)
      
      p_global_complete[total+(1:temp2),] <-matrix(result.list[[2]][[5]],ncol=6)  
    }else if(temp2==0){
      p_global_result[total+(1:temp1),] <- matrix(result.list[[1]][[1]],ncol=6)
      
      p_mglobal_result[total+(1:temp1),] <- matrix(result.list[[1]][[2]],ncol=6)
      p_standard[total+(1:temp1),] <- matrix(result.list[[1]][[3]],ncol=6)
      
      p_global_complete[total+(1:temp1),] <-matrix(result.list[[1]][[5]],ncol=6)  
    }else{
      p_global_result[total+(1:temp),] <-rbind(matrix(result.list[[1]][[1]],ncol=6),
                                               matrix(result.list[[2]][[1]],ncol=6))
      
      p_mglobal_result[total+(1:temp),] <- rbind(matrix(result.list[[1]][[2]],ncol=6),
                                                 matrix(result.list[[2]][[2]],ncol=6))
      p_standard[total+(1:temp),] <- rbind(matrix(result.list[[1]][[3]],ncol=6),
                                           matrix(result.list[[2]][[3]],ncol=6))
      
      p_global_complete[total+(1:temp),] <-rbind(matrix(result.list[[1]][[5]],ncol=6),
                                                 matrix(result.list[[2]][[5]],ncol=6))
      
      #p_poly[total+(1:temp),] <- rbind(matrix(result.list[[1]][[6]],ncol=9),
      #                                matrix(result.list[[2]][[6]],ncol=9))
    }
    
    
    #p_poly[total+(1:temp),] <- rbind(matrix(result.list[[1]][[6]],ncol=9),
    # matrix(result.list[[2]][[6]],ncol=9))
    
    total = total+ temp
    
  }
}

CountPower <- function(p,alpha){
  n <- length(p)
  idx <- which(p<=alpha)
  return(length(idx)/n)
  
}



setwd('/spin1/users/zhangh24/breast_cancer_data_analysis/')
filedir <- '/spin1/users/zhangh24/breast_cancer_data_analysis/simulation/power/result/'

files <- dir(filedir,pattern="poly_high",full.names=T)
total <- 0
for(i1 in 1:500){
  print(i1)
  file = paste0("/spin1/users/zhangh24/breast_cancer_data_analysis/simulation/power/result//poly_high",i1,".Rdata")
  if(file%in%files==T){
    load(paste0("./simulation/power/result/poly_high",i1,".Rdata")) 
    total = total+ length(p_poly)/6 
    
  }
}

#total = total*2

p_poly_result <- matrix(0,total,6)

total <- 0

for(i1 in 1:500){
  print(i1)
  file = paste0("/spin1/users/zhangh24/breast_cancer_data_analysis/simulation/power/result//poly_high",i1,".Rdata")
  if(file%in%files==T){
    load(paste0("./simulation/power/result/poly_high",i1,".Rdata"))
    temp = length(p_poly)/6
    p_poly_result[total+(1:temp),] = matrix(p_poly,ncol=6)
    #p_poly[total+(1:temp),] <- rbind(matrix(result.list[[1]][[6]],ncol=9),
    # matrix(result.list[[2]][[6]],ncol=9))
    
    total = total+ temp
    
  }
}








apply(p_global_result,2,function(x){CountPower(x,10^-3)})
apply(p_mglobal_result,2,function(x){CountPower(x,10^-3)})
apply(p_standard,2,function(x){CountPower(x,10^-3)})
apply(p_global_complete,2,function(x){CountPower(x,10^-3)})
apply(p_poly_result,2,function(x){CountPower(x,10^-3)})


result <- cbind(apply(p_global_result,2,function(x){CountPower(x,10^-3)}),
                apply(p_mglobal_result,2,function(x){CountPower(x,10^-3)}),
                apply(p_standard,2,function(x){CountPower(x,10^-3)}),
                apply(p_global_complete,2,function(x){CountPower(x,10^-3)}),
                apply(p_poly_result,2,function(x){CountPower(x,10^-3)}))
result.1 <- result

# 
# result <- cbind(apply(p_global_result,2,function(x){CountPower(x,10^-3)}),
#                 apply(p_mglobal_result,2,function(x){CountPower(x,10^-3)}),
#                 apply(p_standard,2,function(x){CountPower(x,10^-3)}),
#                 apply(p_global_complete,2,function(x){CountPower(x,10^-3)}),
#                 apply(p_poly,2,function(x){CountPower(x,10^-3)}))













setwd('/spin1/users/zhangh24/breast_cancer_data_analysis/')
filedir <- '/spin1/users/zhangh24/breast_cancer_data_analysis/simulation/power/result/'
files <- dir(filedir,pattern="simu_high",full.names=T)
total <- 0
for(i1 in 4001:5000){
  print(i1)
  file = paste0("/spin1/users/zhangh24/breast_cancer_data_analysis/simulation/power/result//simu_high",i1,".Rdata")
  if(file%in%files==T){
    load(paste0("./simulation/power/result/simu_high",i1,".Rdata")) 
    total = total+ length(result.list[[1]][[1]])/3 + length(result.list[[2]][[1]])/3
    
  }
}

#total = total*2

p_global_result <- matrix(0,total,3)
p_mglobal_result <- matrix(0,total,3)
p_standard <- matrix(0,total,3)
p_global_complete <- matrix(0,total,3)
#p_poly <- matrix(0,total,9)

total <- 0

for(i1 in 4001:5000){
  print(i1)
  file = paste0("/spin1/users/zhangh24/breast_cancer_data_analysis/simulation/power/result//simu_high",i1,".Rdata")
  if(file%in%files==T){
    load(paste0("./simulation/power/result/simu_high",i1,".Rdata")) 
    temp1 = length(result.list[[1]][[1]])/3
    temp2 = length(result.list[[2]][[1]])/3
    temp = temp1+temp2
    if(temp1==0){
      p_global_result[total+(1:temp2),] <- matrix(result.list[[2]][[1]],ncol=3)
      
      p_mglobal_result[total+(1:temp2),] <- matrix(result.list[[2]][[2]],ncol=3)
      p_standard[total+(1:temp2),] <- matrix(result.list[[2]][[3]],ncol=3)
      
      p_global_complete[total+(1:temp2),] <-matrix(result.list[[2]][[5]],ncol=3)  
    }else if(temp2==0){
      p_global_result[total+(1:temp1),] <- matrix(result.list[[1]][[1]],ncol=3)
      
      p_mglobal_result[total+(1:temp1),] <- matrix(result.list[[1]][[2]],ncol=3)
      p_standard[total+(1:temp1),] <- matrix(result.list[[1]][[3]],ncol=3)
      
      p_global_complete[total+(1:temp1),] <-matrix(result.list[[1]][[5]],ncol=3)  
    }else{
      p_global_result[total+(1:temp),] <-rbind(matrix(result.list[[1]][[1]],ncol=3),
                                               matrix(result.list[[2]][[1]],ncol=3))
      
      p_mglobal_result[total+(1:temp),] <- rbind(matrix(result.list[[1]][[2]],ncol=3),
                                                 matrix(result.list[[2]][[2]],ncol=3))
      p_standard[total+(1:temp),] <- rbind(matrix(result.list[[1]][[3]],ncol=3),
                                           matrix(result.list[[2]][[3]],ncol=3))
      
      p_global_complete[total+(1:temp),] <-rbind(matrix(result.list[[1]][[5]],ncol=3),
                                                 matrix(result.list[[2]][[5]],ncol=3))
      
      #p_poly[total+(1:temp),] <- rbind(matrix(result.list[[1]][[6]],ncol=9),
      #                                matrix(result.list[[2]][[6]],ncol=9))
    }
    
    
    #p_poly[total+(1:temp),] <- rbind(matrix(result.list[[1]][[6]],ncol=9),
    # matrix(result.list[[2]][[6]],ncol=9))
    
    total = total+ temp
    
  }
}

CountPower <- function(p,alpha){
  n <- length(p)
  idx <- which(p<=alpha)
  return(length(idx)/n)
  
}



setwd('/spin1/users/zhangh24/breast_cancer_data_analysis/')
filedir <- '/spin1/users/zhangh24/breast_cancer_data_analysis/simulation/power/result/'

files <- dir(filedir,pattern="poly_high",full.names=T)
total <- 0
for(i1 in 501:1000){
  print(i1)
  file = paste0("/spin1/users/zhangh24/breast_cancer_data_analysis/simulation/power/result//poly_high",i1,".Rdata")
  if(file%in%files==T){
    load(paste0("./simulation/power/result/poly_high",i1,".Rdata")) 
    total = total+ length(p_poly)/3
    
  }
}

#total = total*2

p_poly_result <- matrix(0,total,3)

total <- 0

for(i1 in 501:1000){
  print(i1)
  file = paste0("/spin1/users/zhangh24/breast_cancer_data_analysis/simulation/power/result//poly_high",i1,".Rdata")
  if(file%in%files==T){
    load(paste0("./simulation/power/result/poly_high",i1,".Rdata"))
    temp = length(p_poly)/3
    p_poly_result[total+(1:temp),] = matrix(p_poly,ncol=3)
    #p_poly[total+(1:temp),] <- rbind(matrix(result.list[[1]][[6]],ncol=9),
    # matrix(result.list[[2]][[6]],ncol=9))
    
    total = total+ temp
    
  }
}


apply(p_global_result,2,function(x){CountPower(x,10^-3)})
apply(p_mglobal_result,2,function(x){CountPower(x,10^-3)})
apply(p_standard,2,function(x){CountPower(x,10^-3)})
apply(p_global_complete,2,function(x){CountPower(x,10^-3)})
apply(p_poly_result,2,function(x){CountPower(x,10^-3)})


result <- cbind(apply(p_global_result,2,function(x){CountPower(x,10^-3)}),
                apply(p_mglobal_result,2,function(x){CountPower(x,10^-3)}),
                apply(p_standard,2,function(x){CountPower(x,10^-3)}),
                apply(p_global_complete,2,function(x){CountPower(x,10^-3)}),
                apply(p_poly_result,2,function(x){CountPower(x,10^-3)}))


result.2 <- result

result <- matrix(0,9,5)
result[c(2,3,5,6,8,9),] <- result.1
result[c(1,4,7),] <- result.2
write.csv(result,file=paste0("./simulation/power/result/power_high.result.csv") )







#apply(p_poly,2,function(x){CountPower(x,10^-3)})
#CountTypeOne(p_global_result,10^-4)




rbind(matrix(result.list[[1]][[3]],ncol=9),
      matrix(result.list[[2]][[3]],ncol=9))


#result <- list(p_global_result,p_mglobal_result,p_standard,p_mglobal_result,p_global_complete,p_poly)
