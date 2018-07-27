setwd('/spin1/users/zhangh24/breast_cancer_data_analysis/')
filedir <- '/spin1/users/zhangh24/breast_cancer_data_analysis/simulation/power/result/'
files <- dir(filedir,pattern="simu_result",full.names=T)
total <- 0
for(i1 in 1:5999){
  print(i1)
  file = paste0("/spin1/users/zhangh24/breast_cancer_data_analysis/simulation/power/result//simu_result",i1,".Rdata")
  if(file%in%files==T){
    load(paste0("./simulation/power/result//simu_result",i1,".Rdata")) 
    total = total+ length(result.list[[1]][[1]])/9 + length(result.list[[2]][[1]])/9
    
  }
  }

#total = total*2

p_global_result <- matrix(0,total,9)
p_mglobal_result <- matrix(0,total,9)
p_standard <- matrix(0,total,9)
p_global_complete <- matrix(0,total,9)
#p_poly <- matrix(0,total,9)

total <- 0

for(i1 in 1:5999){
  print(i1)
  file = paste0("/spin1/users/zhangh24/breast_cancer_data_analysis/simulation/power/result//simu_result",i1,".Rdata")
  if(file%in%files==T){
    load(paste0("./simulation/power/result//simu_result",i1,".Rdata"))
    temp1 = length(result.list[[1]][[1]])/9
    temp2 = length(result.list[[2]][[1]])/9
    temp = temp1+temp2
    if(temp1==0){
      p_global_result[total+(1:temp2),] <- matrix(result.list[[2]][[1]],ncol=9)
      
      p_mglobal_result[total+(1:temp2),] <- matrix(result.list[[2]][[2]],ncol=9)
      p_standard[total+(1:temp2),] <- matrix(result.list[[2]][[3]],ncol=9)
      
      p_global_complete[total+(1:temp2),] <-matrix(result.list[[2]][[5]],ncol=9)  
    }else if(temp2==0){
      p_global_result[total+(1:temp1),] <- matrix(result.list[[1]][[1]],ncol=9)
      
      p_mglobal_result[total+(1:temp1),] <- matrix(result.list[[1]][[2]],ncol=9)
      p_standard[total+(1:temp1),] <- matrix(result.list[[1]][[3]],ncol=9)
      
      p_global_complete[total+(1:temp1),] <-matrix(result.list[[1]][[5]],ncol=9)  
    }else{
      p_global_result[total+(1:temp),] <-rbind(matrix(result.list[[1]][[1]],ncol=9),
                                               matrix(result.list[[2]][[1]],ncol=9))
      
      p_mglobal_result[total+(1:temp),] <- rbind(matrix(result.list[[1]][[2]],ncol=9),
                                                 matrix(result.list[[2]][[2]],ncol=9))
      p_standard[total+(1:temp),] <- rbind(matrix(result.list[[1]][[3]],ncol=9),
                                           matrix(result.list[[2]][[3]],ncol=9))
      
      p_global_complete[total+(1:temp),] <-rbind(matrix(result.list[[1]][[5]],ncol=9),
                                                 matrix(result.list[[2]][[5]],ncol=9))
      
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
files <- dir(filedir,pattern="simu_result",full.names=T)
total <- 0
for(i1 in 4000:6000){
  print(i1)
  file = paste0("/spin1/users/zhangh24/breast_cancer_data_analysis/simulation/power/result//simu_result",i1,".Rdata")
  if(file%in%files==T){
    load(paste0("./simulation/power/result//simu_result",i1,".Rdata")) 
    total = total+ length(result.list[[1]][[1]])/9 + length(result.list[[2]][[1]])/9
    
  }
}

#total = total*2

# p_global_result <- matrix(0,total,9)
# p_mglobal_result <- matrix(0,total,9)
# p_standard <- matrix(0,total,9)
# p_global_complete <- matrix(0,total,9)
p_poly <- matrix(0,total,9)

total <- 0

for(i1 in 4000:6000){
  print(i1)
  file = paste0("/spin1/users/zhangh24/breast_cancer_data_analysis/simulation/power/result//simu_result",i1,".Rdata")
  if(file%in%files==T){
    load(paste0("./simulation/power/result//simu_result",i1,".Rdata"))
    temp1 = length(result.list[[1]][[1]])/9
    temp2 = length(result.list[[2]][[1]])/9
    temp = temp1+temp2
    if(temp1==0){
      p_poly[total+(1:temp2),] <- matrix(result.list[[2]][[6]],ncol=9)  
    }else if(temp2==0){
      p_poly[total+(1:temp1),] <- matrix(result.list[[1]][[6]],ncol=9)
                                       
    }else{
    
      
      p_poly[total+(1:temp),] <- rbind(matrix(result.list[[1]][[6]],ncol=9),
                                      matrix(result.list[[2]][[6]],ncol=9))
    }
    
    
    #p_poly[total+(1:temp),] <- rbind(matrix(result.list[[1]][[6]],ncol=9),
    # matrix(result.list[[2]][[6]],ncol=9))
    
    total = total+ temp
    
  }
}










apply(p_global_result,2,function(x){CountPower(x,10^-3)})
apply(p_mglobal_result,2,function(x){CountPower(x,10^-3)})
apply(p_standard,2,function(x){CountPower(x,10^-3)})
apply(p_global_complete,2,function(x){CountPower(x,10^-3)})
apply(p_poly,2,function(x){CountPower(x,10^-3)})


result <- cbind(apply(p_global_result,2,function(x){CountPower(x,10^-3)}),
                apply(p_mglobal_result,2,function(x){CountPower(x,10^-3)}),
                apply(p_standard,2,function(x){CountPower(x,10^-3)}),
                apply(p_global_complete,2,function(x){CountPower(x,10^-3)}),
                apply(p_poly,2,function(x){CountPower(x,10^-3)}))


result.1 <- result
#write.csv(result,file=paste0("./simulation/power/result/power.simulation.result.csv") )

















setwd('/spin1/users/zhangh24/breast_cancer_data_analysis/')
filedir <- '/spin1/users/zhangh24/breast_cancer_data_analysis/simulation/power/result/'
files <- dir(filedir,pattern="simu_result",full.names=T)
total <- 0
for(i1 in 6000:7000){
  print(i1)
  file = paste0("/spin1/users/zhangh24/breast_cancer_data_analysis/simulation/power/result//simu_result",i1,".Rdata")
  if(file%in%files==T){
    load(paste0("./simulation/power/result//simu_result",i1,".Rdata")) 
    total = total+ length(result.list[[1]][[1]])/9 + length(result.list[[2]][[1]])/9
    
  }
}

#total = total*2

p_global_result <- matrix(0,total,9)
p_mglobal_result <- matrix(0,total,9)
p_standard <- matrix(0,total,9)
p_global_complete <- matrix(0,total,9)
p_poly <- matrix(0,total,9)

total <- 0

for(i1 in 6000:7000){
  print(i1)
  file = paste0("/spin1/users/zhangh24/breast_cancer_data_analysis/simulation/power/result//simu_result",i1,".Rdata")
  if(file%in%files==T){
    load(paste0("./simulation/power/result//simu_result",i1,".Rdata"))
    temp1 = length(result.list[[1]][[1]])/9
    temp2 = length(result.list[[2]][[1]])/9
    temp = temp1+temp2
    if(temp1==0){
      p_global_result[total+(1:temp2),] <- matrix(result.list[[2]][[1]],ncol=9)
      
      p_mglobal_result[total+(1:temp2),] <- matrix(result.list[[2]][[2]],ncol=9)
      p_standard[total+(1:temp2),] <- matrix(result.list[[2]][[3]],ncol=9)
      
      p_global_complete[total+(1:temp2),] <-matrix(result.list[[2]][[5]],ncol=9)  
      p_poly[total+(1:temp2),] <- matrix(result.list[[2]][[6]],ncol=9)  
    }else if(temp2==0){
      p_global_result[total+(1:temp1),] <- matrix(result.list[[1]][[1]],ncol=9)
      
      p_mglobal_result[total+(1:temp1),] <- matrix(result.list[[1]][[2]],ncol=9)
      p_standard[total+(1:temp1),] <- matrix(result.list[[1]][[3]],ncol=9)
      
      p_global_complete[total+(1:temp1),] <-matrix(result.list[[1]][[5]],ncol=9)  
      p_poly[total+(1:temp1),] <- matrix(result.list[[1]][[6]],ncol=9)
    }else{
      p_global_result[total+(1:temp),] <-rbind(matrix(result.list[[1]][[1]],ncol=9),
                                               matrix(result.list[[2]][[1]],ncol=9))
      
      p_mglobal_result[total+(1:temp),] <- rbind(matrix(result.list[[1]][[2]],ncol=9),
                                                 matrix(result.list[[2]][[2]],ncol=9))
      p_standard[total+(1:temp),] <- rbind(matrix(result.list[[1]][[3]],ncol=9),
                                           matrix(result.list[[2]][[3]],ncol=9))
      
      p_global_complete[total+(1:temp),] <-rbind(matrix(result.list[[1]][[5]],ncol=9),
                                                 matrix(result.list[[2]][[5]],ncol=9))
      
      p_poly[total+(1:temp),] <- rbind(matrix(result.list[[1]][[6]],ncol=9),
                                     matrix(result.list[[2]][[6]],ncol=9))
    }
    
    
    #p_poly[total+(1:temp),] <- rbind(matrix(result.list[[1]][[6]],ncol=9),
    # matrix(result.list[[2]][[6]],ncol=9))
    
    total = total+ temp
    
  }
}

























apply(p_global_result,2,function(x){CountPower(x,10^-3)})
apply(p_mglobal_result,2,function(x){CountPower(x,10^-3)})
apply(p_standard,2,function(x){CountPower(x,10^-3)})
apply(p_global_complete,2,function(x){CountPower(x,10^-3)})
apply(p_poly,2,function(x){CountPower(x,10^-3)})


result <- cbind(apply(p_global_result,2,function(x){CountPower(x,10^-3)}),
                apply(p_mglobal_result,2,function(x){CountPower(x,10^-3)}),
                apply(p_standard,2,function(x){CountPower(x,10^-3)}),
                apply(p_global_complete,2,function(x){CountPower(x,10^-3)}),
                apply(p_poly,2,function(x){CountPower(x,10^-3)}))

result.2 <- result

result <- result.1
result[c(1,4,7),] <- result.2[1:3,]

write.csv(result,file=paste0("./simulation/power/result/power.simulation.result.csv") )


result.low <- result
result.low[,2]/result.low[,3]
result.low[,2]/result.low[,4]
#apply(p_poly,2,function(x){CountPower(x,10^-3)})
#CountTypeOne(p_global_result,10^-4)




# rbind(matrix(result.list[[1]][[3]],ncol=9),
#       matrix(result.list[[2]][[3]],ncol=9))


#result <- list(p_global_result,p_mglobal_result,p_standard,p_mglobal_result,p_global_complete,p_poly)
