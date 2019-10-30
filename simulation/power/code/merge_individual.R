#merge power results two-way interactions with the effect size as 0.08 level
setwd('/data/zhangh24/breast_cancer_data_analysis/')
filedir <- '/data/zhangh24/breast_cancer_data_analysis/simulation/power/result/'
files <- dir(filedir,pattern="simu_indi_",full.names=T)
total <- 0
#n.loop = number of simulation setting* number of sample size setting
n.loop <- 9
for(i1 in 1:2000){
  print(i1)
  file = paste0("/data/zhangh24/breast_cancer_data_analysis/simulation/power/result//simu_indi_",i1,".Rdata")
  if(file%in%files==T){
    load(file) 
    #result.list is a list of pvalue
    #three different simulation settings: 1. no heterogneity 2. one tumor heter 3. multiple tumor heterogeneity
    #3 different sample size were implemented 5000, 50,000 and 100,000
    #
    #its the output of a foreach parallele result
    #[[1]] and [[2]] share the same structure
    #[[1]] [[1]] is the vector of p_value from FTOP with additive structure
    #[[1]] [[1]] is a long vector looped by first simulation setting (only 1 under this case), then sample size (inner loop 3), then the number of replicates (third loop)
    #[[1]] [[2]] is the vector of p_value from MTOP with additive structure
    #[[1]] [[3]] is the vector of p_value from ftop with all interaction
    #[[1]] [[4]] is the vector of p_value from MTOP with all interactions
    total = total+ nrow(result.list[[1]][[1]])/n.loop + nrow(result.list[[2]][[1]])/n.loop
    
  }
}

#total = total*2

p_all_ER <- matrix(0,total,n.loop)
p_all_PR <- matrix(0,total,n.loop)
p_all_HER2 <- matrix(0,total,n.loop)
p_all_grade <- matrix(0,total,n.loop)
p_complete_ER <- matrix(0,total,n.loop)
p_complete_PR <- matrix(0,total,n.loop)
p_complete_HER2 <- matrix(0,total,n.loop)
p_complete_grade <- matrix(0,total,n.loop)
#p_poly <- matrix(0,total,9)

total <- 0
#args 1:2000 contains the simulation results for FTOP, MTOP, standard logistic regressionn, complete FTOP
for(i1 in 1:2000){
  print(i1)
  file = paste0("/data/zhangh24/breast_cancer_data_analysis/simulation/power/result//simu_indi_",i1,".Rdata")
  if(file%in%files==T){
    load(file) 
    temp1 = nrow(result.list[[1]][[1]])/n.loop
    temp2 = nrow(result.list[[2]][[1]])/n.loop
    temp = temp1+temp2
  
    p_all_ER[total+(1:temp),] <-rbind(matrix(result.list[[1]][[1]][,1],ncol=n.loop),
                                          matrix(result.list[[2]][[1]][,1],ncol=n.loop))
    
    p_all_PR[total+(1:temp),] <- rbind(matrix(result.list[[1]][[1]][,2],ncol=n.loop),
                                         matrix(result.list[[2]][[1]][,2],ncol=n.loop))  
    
    p_all_HER2[total+(1:temp),] <- rbind(matrix(result.list[[1]][[1]][,3],ncol=n.loop),
                                           matrix(result.list[[2]][[1]][,3],ncol=n.loop))
   
    
    p_all_grade[total+(1:temp),] <- rbind(matrix(result.list[[1]][[1]][,4],ncol=n.loop),
                                         matrix(result.list[[2]][[1]][,4],ncol=n.loop))
    
    p_complete_ER[total+(1:temp),] <-rbind(matrix(result.list[[1]][[1]][,5],ncol=n.loop),
                                      matrix(result.list[[2]][[1]][,5],ncol=n.loop))
    
    p_complete_PR[total+(1:temp),] <-rbind(matrix(result.list[[1]][[1]][,6],ncol=n.loop),
                                           matrix(result.list[[2]][[1]][,6],ncol=n.loop))
    
    p_complete_HER2[total+(1:temp),] <-rbind(matrix(result.list[[1]][[1]][,7],ncol=n.loop),
                                           matrix(result.list[[2]][[1]][,7],ncol=n.loop))
    
    p_complete_grade[total+(1:temp),] <-rbind(matrix(result.list[[1]][[1]][,8],ncol=n.loop),
                                             matrix(result.list[[2]][[1]][,8],ncol=n.loop))
    
    
    
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
thres2 = 0.001

#remove standard polytomous function 
#unstable outliers
#idx <- which(p_poly[,4]==0)
#p_poly = p_poly[-idx,,drop=F]


result1 <- cbind(apply(p_all_ER,2,function(x){CountPower(x,thres)}),
              apply(p_complete_ER,2,function(x){CountPower(x,thres)})
               )


result2 <- cbind(  apply(p_all_PR,2,function(x){CountPower(x,thres2)}),
        apply(p_all_HER2,2,function(x){CountPower(x,thres2)}),
        apply(p_all_grade,2,function(x){CountPower(x,thres2)}),
        apply(p_complete_PR,2,function(x){CountPower(x,thres2)}),
        apply(p_complete_HER2,2,function(x){CountPower(x,thres2)}),
        apply(p_complete_grade,2,function(x){CountPower(x,thres2)}))



result.1 <- result
#write.csv(result,file=paste0("./simulation/power/result/power.simulation.result.csv") )


write.csv(result1,file=paste0("./simulation/power/result/indi.simulation.result.csv") )
write.csv(result2,file = paste0("./simulation/power/result/indi.simulation.result_others.csv"))







