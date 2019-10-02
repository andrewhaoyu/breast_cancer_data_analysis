setwd('/spin1/users/zhangh24/breast_cancer_data_analysis/')
filedir <- '/spin1/users/zhangh24/breast_cancer_data_analysis/simulation/power/result/'
files <- dir(filedir,pattern="simu_result",full.names=T)
total <- 0
n.loop <- 12
for(i1 in 1:2000){
  print(i1)
  file = paste0("/spin1/users/zhangh24/breast_cancer_data_analysis/simulation/power/result//simu_result",i1,".Rdata")
  if(file%in%files==T){
    load(paste0("./simulation/power/result/simu_result",i1,".Rdata")) 
    #result.list is a list of pvalue
    #three different simulation settings: 1. no heterogneity 2. one tumor heter 3. multiple tumor heterogeneity
    #3 different sample size were implemented 5000, 50,000 and 100,000
    #
    #its the output of a foreach parallele result
    #[[1]] and [[2]] share the same structure
    #[[1]] [[1]] is the vector of p_value from FTOP
    #[[1]] [[1]] is a long vector looped by first simulation setting, then sample size (inner loop), 9 different sections
    #[[1]] [[2]] is the vector of p_value from MTOP
    #[[1]] [[3]] is the vector of p_value from standard logistoc regression
    #[[1]] [[4]] is the vector of p_value from MTOP this is because of a previous typo
    #[[1]] [[5]] is the vector of FTOP from complete analysis
    #[[1]] [[6]] is the vector of polytomous model from complete analysis
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
#args 1:2000 contains the simulation results for FTOP, MTOP, standard logistic regressionn, complete FTOP
for(i1 in 1:2000){
  print(i1)
  file = paste0("/spin1/users/zhangh24/breast_cancer_data_analysis/simulation/power/result//simu_result",i1,".Rdata")
  if(file%in%files==T){
    load(paste0("./simulation/power/result//simu_result",i1,".Rdata"))
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

#load results for polytomous model



setwd('/spin1/users/zhangh24/breast_cancer_data_analysis/')
filedir <- '/spin1/users/zhangh24/breast_cancer_data_analysis/simulation/power/result/'
files <- dir(filedir,pattern="poly",full.names=T)
total <- 0
#args 1:2000 contains the results for polytomous
for(i1 in 1:2000){
  print(i1)
  file = paste0("/spin1/users/zhangh24/breast_cancer_data_analysis/simulation/power/result//poly",i1,".Rdata")
  if(file%in%files==T){
    load(file) 
    total = total+ length(result.list[[1]][[1]])/n.loop + length(result.list[[2]][[1]])/n.loop
    
  }
}


p_poly <- matrix(0,total,n.loop)

total <- 0

for(i1 in 1:2000){
  print(i1)
  file = paste0("/spin1/users/zhangh24/breast_cancer_data_analysis/simulation/power/result//poly",i1,".Rdata")
  if(file%in%files==T){
    load(paste0(file))
    temp1 = length(result.list[[1]][[1]])/n.loop
    temp2 = length(result.list[[2]][[1]])/n.loop
    temp = temp1+temp2
    if(temp1==0){
      p_poly[total+(1:temp2),] <- matrix(result.list[[2]][[1]],ncol=n.loop)  
    }else if(temp2==0){
      p_poly[total+(1:temp1),] <- matrix(result.list[[1]][[1]],ncol=n.loop)
                                       
    }else{
    
      
      p_poly[total+(1:temp),] <- rbind(matrix(result.list[[1]][[1]],ncol=n.loop),
                                      matrix(result.list[[2]][[1]],ncol=n.loop))
    }
    
    
    #p_poly[total+(1:temp),] <- rbind(matrix(result.list[[1]][[6]],ncol=9),
    # matrix(result.list[[2]][[6]],ncol=9))
    
    total = total+ temp
    
  }
}










# apply(p_global_result,2,function(x){CountPower(x,10^-3)})
# apply(p_mglobal_result,2,function(x){CountPower(x,10^-3)})
# apply(p_standard,2,function(x){CountPower(x,10^-3)})
# apply(p_global_complete,2,function(x){CountPower(x,10^-3)})
# apply(p_poly,2,function(x){CountPower(x,thres)})

thres = 1E-03
#remove standard polytomous function 
#unstable outliers
#idx <- which(p_poly[,4]==0)
#p_poly = p_poly[-idx,,drop=F]

result <- cbind(apply(p_global_result,2,function(x){CountPower(x,thres)}),
                apply(p_mglobal_result,2,function(x){CountPower(x,thres)}),
                apply(p_standard,2,function(x){CountPower(x,thres)}),
                apply(p_global_complete,2,function(x){CountPower(x,thres)}),
                apply(p_poly,2,function(x){CountPower(x,thres)}))



result.1 <- result
#write.csv(result,file=paste0("./simulation/power/result/power.simulation.result.csv") )


write.csv(result,file=paste0("./simulation/power/result/power.simulation.result.csv") )














# setwd('/spin1/users/zhangh24/breast_cancer_data_analysis/')
# filedir <- '/spin1/users/zhangh24/breast_cancer_data_analysis/simulation/power/result/'
# files <- dir(filedir,pattern="simu_result",full.names=T)
# total <- 0
# for(i1 in 6000:7000){
#   print(i1)
#   file = paste0("/spin1/users/zhangh24/breast_cancer_data_analysis/simulation/power/result//simu_result",i1,".Rdata")
#   if(file%in%files==T){
#     load(paste0("./simulation/power/result//simu_result",i1,".Rdata")) 
#     total = total+ length(result.list[[1]][[1]])/9 + length(result.list[[2]][[1]])/9
#     
#   }
# }
# 
# #total = total*2
# 
# p_global_result <- matrix(0,total,9)
# p_mglobal_result <- matrix(0,total,9)
# p_standard <- matrix(0,total,9)
# p_global_complete <- matrix(0,total,9)
# p_poly <- matrix(0,total,9)
# 
# total <- 0
# 
# for(i1 in 6000:7000){
#   print(i1)
#   file = paste0("/spin1/users/zhangh24/breast_cancer_data_analysis/simulation/power/result/simu_result",i1,".Rdata")
#   if(file%in%files==T){
#     load(paste0("./simulation/power/result//simu_result",i1,".Rdata"))
#     temp1 = length(result.list[[1]][[1]])/9
#     temp2 = length(result.list[[2]][[1]])/9
#     temp = temp1+temp2
#     if(temp1==0){
#       p_global_result[total+(1:temp2),] <- matrix(result.list[[2]][[1]],ncol=9)
#       
#       p_mglobal_result[total+(1:temp2),] <- matrix(result.list[[2]][[2]],ncol=9)
#       p_standard[total+(1:temp2),] <- matrix(result.list[[2]][[3]],ncol=9)
#       
#       p_global_complete[total+(1:temp2),] <-matrix(result.list[[2]][[5]],ncol=9)  
#       p_poly[total+(1:temp2),] <- matrix(result.list[[2]][[6]],ncol=9)  
#     }else if(temp2==0){
#       p_global_result[total+(1:temp1),] <- matrix(result.list[[1]][[1]],ncol=9)
#       
#       p_mglobal_result[total+(1:temp1),] <- matrix(result.list[[1]][[2]],ncol=9)
#       p_standard[total+(1:temp1),] <- matrix(result.list[[1]][[3]],ncol=9)
#       
#       p_global_complete[total+(1:temp1),] <-matrix(result.list[[1]][[5]],ncol=9)  
#       p_poly[total+(1:temp1),] <- matrix(result.list[[1]][[6]],ncol=9)
#     }else{
#       p_global_result[total+(1:temp),] <-rbind(matrix(result.list[[1]][[1]],ncol=9),
#                                                matrix(result.list[[2]][[1]],ncol=9))
#       
#       p_mglobal_result[total+(1:temp),] <- rbind(matrix(result.list[[1]][[2]],ncol=9),
#                                                  matrix(result.list[[2]][[2]],ncol=9))
#       p_standard[total+(1:temp),] <- rbind(matrix(result.list[[1]][[3]],ncol=9),
#                                            matrix(result.list[[2]][[3]],ncol=9))
#       
#       p_global_complete[total+(1:temp),] <-rbind(matrix(result.list[[1]][[5]],ncol=9),
#                                                  matrix(result.list[[2]][[5]],ncol=9))
#       
#       p_poly[total+(1:temp),] <- rbind(matrix(result.list[[1]][[6]],ncol=9),
#                                      matrix(result.list[[2]][[6]],ncol=9))
#     }
#     
#     
#     #p_poly[total+(1:temp),] <- rbind(matrix(result.list[[1]][[6]],ncol=9),
#     # matrix(result.list[[2]][[6]],ncol=9))
#     
#     total = total+ temp
#     
#   }
# }
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# apply(p_global_result,2,function(x){CountPower(x,10^-3)})
# apply(p_mglobal_result,2,function(x){CountPower(x,10^-3)})
# apply(p_standard,2,function(x){CountPower(x,10^-3)})
# apply(p_global_complete,2,function(x){CountPower(x,10^-3)})
# apply(p_poly,2,function(x){CountPower(x,10^-3)})
# 
# thres = 5E-08
# result <- cbind(apply(p_global_result,2,function(x){CountPower(x,thres)}),
#                 apply(p_mglobal_result,2,function(x){CountPower(x,thres)}),
#                 apply(p_standard,2,function(x){CountPower(x,thres)}),
#                 apply(p_global_complete,2,function(x){CountPower(x,thres)}),
#                 apply(p_poly,2,function(x){CountPower(x,thres)}))
# 
# result.2 <- result
# 
# result <- result.1
# result[c(1,4,7),] <- result.2[1:3,]
# 
# write.csv(result,file=paste0("./simulation/power/result/power.simulation.result.csv") )
# 
# 
# result.low <- result
# result.low[,2]/result.low[,3]
# result.low[,2]/result.low[,4]
#apply(p_poly,2,function(x){CountPower(x,10^-3)})
#CountTypeOne(p_global_result,10^-4)




# rbind(matrix(result.list[[1]][[3]],ncol=9),
#       matrix(result.list[[2]][[3]],ncol=9))


#result <- list(p_global_result,p_mglobal_result,p_standard,p_mglobal_result,p_global_complete,p_poly)
