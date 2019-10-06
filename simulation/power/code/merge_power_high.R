setwd('/spin1/users/zhangh24/breast_cancer_data_analysis/')
filedir <- '/spin1/users/zhangh24/breast_cancer_data_analysis/simulation/power/result/'
files <- dir(filedir,pattern="simu_high",full.names=T)
total <- 0
n.loop = 12
for(i1 in 1:12000){
  print(i1)
  file = paste0("/spin1/users/zhangh24/breast_cancer_data_analysis/simulation/power/result//simu_high",i1,".Rdata")
  if(file%in%files==T){
    load(paste0("./simulation/power/result/simu_high",i1,".Rdata")) 
    total = total+ length(result.list[[1]][[1]])/n.loop + length(result.list[[2]][[1]])/n.loop
    
  }
}

#total = total*2

p_global_result <- matrix(0,total,n.loop)
p_mglobal_result <- matrix(0,total,n.loop)
p_standard <- matrix(0,total,n.loop)
p_global_complete <- matrix(0,total,n.loop)
#p_poly <- matrix(0,total,9)

total <- 0

for(i1 in 1:12000){
  print(i1)
  file = paste0("/spin1/users/zhangh24/breast_cancer_data_analysis/simulation/power/result//simu_high",i1,".Rdata")
  if(file%in%files==T){
    load(paste0("./simulation/power/result/simu_high",i1,".Rdata")) 
    temp1 = length(result.list[[1]][[1]])/n.loop
    temp2 = length(result.list[[2]][[1]])/n.loop
    temp = temp1+temp2
    if(temp1==0){
      p_global_result[total+(1:temp2),] <- matrix(result.list[[2]][[1]],ncol=n.loop)
      
      p_mglobal_result[total+(1:temp2),] <- matrix(result.list[[2]][[2]],ncol=n.loop)
      p_standard[total+(1:temp2),] <- matrix(result.list[[2]][[3]],ncol=n.loop)
      
      p_global_complete[total+(1:temp2),] <-matrix(result.list[[2]][[5]],ncol=n.loop)  
    }else if(temp2==0){
      p_global_result[total+(1:temp1),] <- matrix(result.list[[1]][[1]],ncol=n.loop)
      
      p_mglobal_result[total+(1:temp1),] <- matrix(result.list[[1]][[2]],ncol=n.loop)
      p_standard[total+(1:temp1),] <- matrix(result.list[[1]][[3]],ncol=n.loop)
      
      p_global_complete[total+(1:temp1),] <-matrix(result.list[[1]][[5]],ncol=n.loop)  
    }else{
      p_global_result[total+(1:temp),] <-rbind(matrix(result.list[[1]][[1]],ncol=n.loop),
                                               matrix(result.list[[2]][[1]],ncol=n.loop))
      
      p_mglobal_result[total+(1:temp),] <- rbind(matrix(result.list[[1]][[2]],ncol=n.loop),
                                                 matrix(result.list[[2]][[2]],ncol=n.loop))
      p_standard[total+(1:temp),] <- rbind(matrix(result.list[[1]][[3]],ncol=n.loop),
                                           matrix(result.list[[2]][[3]],ncol=n.loop))
      
      p_global_complete[total+(1:temp),] <-rbind(matrix(result.list[[1]][[5]],ncol=n.loop),
                                                 matrix(result.list[[2]][[5]],ncol=n.loop))
      
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
for(i1 in 1:2000){
  print(i1)
  file = paste0("/spin1/users/zhangh24/breast_cancer_data_analysis/simulation/power/result//poly_high",i1,".Rdata")
  if(file%in%files==T){
    load(paste0("./simulation/power/result/poly_high",i1,".Rdata")) 
    total = total+ length(p_poly)/n.loop
    
  }
}

#total = total*2

p_poly_result <- matrix(0,total,n.loop)

total <- 0
for(i1 in 1:2000){
  print(i1)
  file = paste0("/spin1/users/zhangh24/breast_cancer_data_analysis/simulation/power/result//poly_high",i1,".Rdata")
  if(file%in%files==T){
    load(paste0("./simulation/power/result/poly_high",i1,".Rdata"))
    temp = length(p_poly)/n.loop
    #remove polytmous regression unstable converge
    p_poly[p_poly==0] = NA
    p_poly_result[total+(1:temp),] = matrix(p_poly,ncol=n.loop)
    #p_poly[total+(1:temp),] <- rbind(matrix(result.list[[1]][[6]],ncol=9),
    # matrix(result.list[[2]][[6]],ncol=9))
    
    total = total+ temp
    
  }
}


thres = 5E-08
#remove standard polytomous function 
#unstable outliers
#idx <- which(p_poly[,4]==0)
#p_poly = p_poly[-idx,,drop=F]

result <- cbind(apply(p_global_result,2,function(x){CountPower(x,thres)}),
                apply(p_mglobal_result,2,function(x){CountPower(x,thres)}),
                apply(p_standard,2,function(x){CountPower(x,thres)}),
                apply(p_global_complete,2,function(x){CountPower(x,thres)}),
                apply(p_poly_result,2,function(x){CountPower(x,thres)}))

result.1 <- result

write.csv(result,file=paste0("./simulation/power/result/power_high.result.csv") )











#load results for effect with 0.25
setwd('/spin1/users/zhangh24/breast_cancer_data_analysis/')
filedir <- '/spin1/users/zhangh24/breast_cancer_data_analysis/simulation/power/result/'
files <- dir(filedir,pattern="simu_high_0.25_",full.names=T)
total <- 0
n.loop = 3
for(i1 in 1:12000){
  print(i1)
  file = paste0("/spin1/users/zhangh24/breast_cancer_data_analysis/simulation/power/result//simu_high_0.25_",i1,".Rdata")
  if(file%in%files==T){
    load(file) 
    total = total+ length(result.list[[1]][[1]])/n.loop + length(result.list[[2]][[1]])/n.loop
    
  }
}

#total = total*2

p_global_result <- matrix(0,total,n.loop)
p_mglobal_result <- matrix(0,total,n.loop)
p_standard <- matrix(0,total,n.loop)
p_global_complete <- matrix(0,total,n.loop)
#p_poly <- matrix(0,total,9)

total <- 0

for(i1 in 1:12000){
  print(i1)
  file = paste0("/spin1/users/zhangh24/breast_cancer_data_analysis/simulation/power/result//simu_high_0.25_",i1,".Rdata")
  if(file%in%files==T){
    load(file) 
    temp1 = length(result.list[[1]][[1]])/n.loop
    temp2 = length(result.list[[2]][[1]])/n.loop
    temp = temp1+temp2
    if(temp1==0){
      p_global_result[total+(1:temp2),] <- matrix(result.list[[2]][[1]],ncol=n.loop)
      
      p_mglobal_result[total+(1:temp2),] <- matrix(result.list[[2]][[2]],ncol=n.loop)
      p_standard[total+(1:temp2),] <- matrix(result.list[[2]][[3]],ncol=n.loop)
      
      p_global_complete[total+(1:temp2),] <-matrix(result.list[[2]][[5]],ncol=n.loop)  
    }else if(temp2==0){
      p_global_result[total+(1:temp1),] <- matrix(result.list[[1]][[1]],ncol=n.loop)
      
      p_mglobal_result[total+(1:temp1),] <- matrix(result.list[[1]][[2]],ncol=n.loop)
      p_standard[total+(1:temp1),] <- matrix(result.list[[1]][[3]],ncol=n.loop)
      
      p_global_complete[total+(1:temp1),] <-matrix(result.list[[1]][[5]],ncol=n.loop)  
    }else{
      p_global_result[total+(1:temp),] <-rbind(matrix(result.list[[1]][[1]],ncol=n.loop),
                                               matrix(result.list[[2]][[1]],ncol=n.loop))
      
      p_mglobal_result[total+(1:temp),] <- rbind(matrix(result.list[[1]][[2]],ncol=n.loop),
                                                 matrix(result.list[[2]][[2]],ncol=n.loop))
      p_standard[total+(1:temp),] <- rbind(matrix(result.list[[1]][[3]],ncol=n.loop),
                                           matrix(result.list[[2]][[3]],ncol=n.loop))
      
      p_global_complete[total+(1:temp),] <-rbind(matrix(result.list[[1]][[5]],ncol=n.loop),
                                                 matrix(result.list[[2]][[5]],ncol=n.loop))
      
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

files <- dir(filedir,pattern="poly_high_0.25_",full.names=T)
total <- 0
for(i1 in 1:2000){
  print(i1)
  file = paste0("/spin1/users/zhangh24/breast_cancer_data_analysis/simulation/power/result//poly_high_0.25_",i1,".Rdata")
  if(file%in%files==T){
    load(file) 
    total = total+ length(p_poly)/n.loop
    
  }
}

#total = total*2

p_poly_result <- matrix(0,total,n.loop)

total <- 0
for(i1 in 1:2000){
  print(i1)
  file = paste0("/spin1/users/zhangh24/breast_cancer_data_analysis/simulation/power/result//poly_high_0.25_",i1,".Rdata")
  if(file%in%files==T){
    load(file)
    temp = length(p_poly)/n.loop
    #remove polytmous regression unstable converge
    p_poly[p_poly==0] = NA
    p_poly_result[total+(1:temp),] = matrix(p_poly,ncol=n.loop)
    #p_poly[total+(1:temp),] <- rbind(matrix(result.list[[1]][[6]],ncol=9),
    # matrix(result.list[[2]][[6]],ncol=9))
    
    total = total+ temp
    
  }
}


thres = 5E-08
#remove standard polytomous function 
#unstable outliers
#idx <- which(p_poly[,4]==0)
#p_poly = p_poly[-idx,,drop=F]

result <- cbind(apply(p_global_result,2,function(x){CountPower(x,thres)}),
                apply(p_mglobal_result,2,function(x){CountPower(x,thres)}),
                apply(p_standard,2,function(x){CountPower(x,thres)}),
                apply(p_global_complete,2,function(x){CountPower(x,thres)}),
                apply(p_poly_result,2,function(x){CountPower(x,thres)}))

result.1 <- result

write.csv(result,file=paste0("./simulation/power/result/power_high_0.25_result.csv") )
