setwd("/spin1/users/zhangh24/breast_cancer_data_analysis/simulation/type_one_error/result")
simulation <- 1000
simulation.result <- matrix(0,simulation*100,6)
total <- 0
for(i in 1:1000){
  print(i)
  load(paste0("pvalue",i,".Rdata"))
  temp <- nrow(p.value.simulation)
  simulation.result[(1:temp)+total,] <- p.value.simulation
  total <- temp+total
}
dim(simulation.result)
hist(simulation.result[,1])
library(qqman)
p.value <- simulation.result[,1]
qqman::qq(simulation.result[,6])

library(lattice)
qqmath(~-log10(p.value),
       distribution=function(x){-log10(qunif(1-x))}
       )
save(simulation.result,file="/spin1/users/zhangh24/breast_cancer_data_analysis/simulation/type_one_error/result/type_one_error_combine.Rdata")

load("breast_cancer_data_analysis/simulation/type_one_error/result/type_one_error_combine.Rdata")


my.pvalues<-runif(10000)

library(lattice);
qqmath(~-log10(my.pvalues),
       distribution=function(x){-log10(qunif(1-x))}
);
