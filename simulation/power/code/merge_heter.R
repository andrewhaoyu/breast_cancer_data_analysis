#merge power results with the effect size as 0.08 level
setwd('/spin1/users/zhangh24/breast_cancer_data_analysis/')
filedir <- '/spin1/users/zhangh24/breast_cancer_data_analysis/simulation/power/result/'
files <- dir(filedir,pattern="heter_result_",full.names=T)
total <- 0
#n.loop = number of simulation setting* number of sample size setting
n.loop <- 9
for(i1 in 1:2000){
  print(i1)
  file = paste0("/spin1/users/zhangh24/breast_cancer_data_analysis/simulation/power/result//heter_result_",i1,".Rdata")
  if(file%in%files==T){
    load(file) 
    total = total+ length(result.list[[1]][[1]])/n.loop + length(result.list[[2]][[1]])/n.loop
    
  }
}

#total = total*2

p_global_result <- matrix(0,total,n.loop)
p_mglobal_result <- matrix(0,total,n.loop)
p_global_complete <- matrix(0,total,n.loop)
#p_poly <- matrix(0,total,9)

total <- 0
#args 1:2000 contains the simulation results for FTOP, MTOP, standard logistic regressionn, complete FTOP
n.sizes <- 3
for(i1 in 1:2000){
  print(i1)
  file = paste0("/spin1/users/zhangh24/breast_cancer_data_analysis/simulation/power/result//heter_result_",i1,".Rdata")
  if(file%in%files==T){
    load(file)
    #due to deletion of setting 3
    temp1 = length(result.list[[1]][[1]])/n.loop
    temp2 = length(result.list[[2]][[1]])/n.loop
    temp = temp1+temp2
      p_global_result[total+(1:temp),] <-rbind(matrix(result.list[[1]][[1]],ncol=n.loop),
                                               matrix(result.list[[2]][[1]],ncol=n.loop))
      
      p_mglobal_result[total+(1:temp),] <- rbind(matrix(result.list[[1]][[2]],ncol=n.loop),
                                                 matrix(result.list[[2]][[2]],ncol=n.loop))
      
      p_global_complete[total+(1:temp),] <-rbind(matrix(result.list[[1]][[3]],ncol=n.loop),
                                                 matrix(result.list[[2]][[3]],ncol=n.loop))

    total = total+ temp
    
  }
}

CountPower <- function(p,alpha){
  n <- length(p)
  idx <- which(p<=alpha)
  return(length(idx)/n)
  
}






# apply(p_global_result,2,function(x){CountPower(x,10^-3)})
# apply(p_mglobal_result,2,function(x){CountPower(x,10^-3)})
# apply(p_standard,2,function(x){CountPower(x,10^-3)})
# apply(p_global_complete,2,function(x){CountPower(x,10^-3)})
# apply(p_poly,2,function(x){CountPower(x,thres)})

thres = 5E-08
#remove standard polytomous function 
#unstable outliers
#idx <- which(p_poly[,4]==0)
#p_poly = p_poly[-idx,,drop=F]

result <- cbind(apply(p_global_result,2,function(x){CountPower(x,thres)}),
                apply(p_mglobal_result,2,function(x){CountPower(x,thres)}),
                apply(p_global_complete,2,function(x){CountPower(x,thres)}))

#due to extra simulation parameter, the row 4-6 and row 7-9 are under the same setting.
result.temp <- (result[4:6,]+result[7:9,])/2
result[]
result.1 <- rbind(result[1:3,],result.temp)
#write.csv(result,file=paste0("./simulation/power/result/power.simulation.result.csv") )


write.csv(result.1,file=paste0("./simulation/power/result/heter_simulation.result.csv") )




setwd('/spin1/users/zhangh24/breast_cancer_data_analysis/')
filedir <- '/spin1/users/zhangh24/breast_cancer_data_analysis/simulation/power/result/'
files <- dir(filedir,pattern="heter_high_",full.names=T)
total <- 0
#n.loop = number of simulation setting* number of sample size setting
n.loop <- 9
for(i1 in 1:4000){
  print(i1)
  file = paste0("/spin1/users/zhangh24/breast_cancer_data_analysis/simulation/power/result//heter_high_",i1,".Rdata")
  if(file%in%files==T){
    load(file) 
    total = total+ length(result.list[[1]][[1]])/n.loop + length(result.list[[2]][[1]])/n.loop
    
  }
}

#total = total*2

p_global_result <- matrix(0,total,n.loop)
p_mglobal_result <- matrix(0,total,n.loop)
p_global_complete <- matrix(0,total,n.loop)
#p_poly <- matrix(0,total,9)

total <- 0
#args 1:2000 contains the simulation results for FTOP, MTOP, standard logistic regressionn, complete FTOP
n.sizes <- 3
for(i1 in 1:4000){
  print(i1)
  file = paste0("/spin1/users/zhangh24/breast_cancer_data_analysis/simulation/power/result//heter_high_",i1,".Rdata")
  if(file%in%files==T){
    load(file)
    #due to deletion of setting 3
    temp1 = length(result.list[[1]][[1]])/n.loop
    temp2 = length(result.list[[2]][[1]])/n.loop
    temp = temp1+temp2
    p_global_result[total+(1:temp),] <-rbind(matrix(result.list[[1]][[1]],ncol=n.loop),
                                             matrix(result.list[[2]][[1]],ncol=n.loop))
    
    p_mglobal_result[total+(1:temp),] <- rbind(matrix(result.list[[1]][[2]],ncol=n.loop),
                                               matrix(result.list[[2]][[2]],ncol=n.loop))
    
    p_global_complete[total+(1:temp),] <-rbind(matrix(result.list[[1]][[3]],ncol=n.loop),
                                               matrix(result.list[[2]][[3]],ncol=n.loop))
    
    total = total+ temp
    
  }
}

CountPower <- function(p,alpha){
  n <- length(p)
  idx <- which(p<=alpha)
  return(length(idx)/n)
  
}






# apply(p_global_result,2,function(x){CountPower(x,10^-3)})
# apply(p_mglobal_result,2,function(x){CountPower(x,10^-3)})
# apply(p_standard,2,function(x){CountPower(x,10^-3)})
# apply(p_global_complete,2,function(x){CountPower(x,10^-3)})
# apply(p_poly,2,function(x){CountPower(x,thres)})

thres = 5E-08
#remove standard polytomous function 
#unstable outliers
#idx <- which(p_poly[,4]==0)
#p_poly = p_poly[-idx,,drop=F]

result <- cbind(apply(p_global_result,2,function(x){CountPower(x,thres)}),
                apply(p_mglobal_result,2,function(x){CountPower(x,thres)}),
                apply(p_global_complete,2,function(x){CountPower(x,thres)}))

#due to extra simulation parameter, the row 4-6 and row 7-9 are under the same setting.
result.temp <- (result[4:6,]+result[7:9,])/2

result.2 <- rbind(result[1:3,],result.temp)
#write.csv(result,file=paste0("./simulation/power/result/power.simulation.result.csv") )


write.csv(result.2,file=paste0("./simulation/power/result/heter_high.result.csv") )











