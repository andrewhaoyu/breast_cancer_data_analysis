#goal: weighted combination of PRS methods
setwd('/spin1/users/zhangh24/breast_cancer_data_analysis')
for(pop.ind in 1:3){
  load(paste0("./multi_ethnic/result/LDP.result_",pop.ind))  
}

p.thr <- c(10^-8,10^-7,10^-6,10^-5,10^-4,10^-3,10^-2,0.1,0.3,0.5)
