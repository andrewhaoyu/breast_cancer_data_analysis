args = commandArgs(trailingOnly = T)
i1 = as.numeric(args[[1]])
setwd("/spin1/users/zhangh24/breast_cancer_data_analysis/")
paste0("./genetic_correlation/ICOG/result/icog.onco.merge.Rdata")
library(bc2)
size =1000
start.end<- startend(n,size,i1)
start <- start.end[1]
end <- start.end[2]

result.sub <- matrix(0,end-start+1,30)
total <- 0
for(i in start:end){
  print(i)
  logodds.icog <- as.numeric(as.vector(icog.onco.merge[i,11:15]))
  sigma.icog <- matrix(as.numeric(icog.onco.merge[i,16:40]),5,5)
  logodds.onco <- as.numeric(as.vector(icog.onco.merge[i,50:54]))
  sigma.onco <- matrix(as.numeric(icog.onco.merge[i,55:79]),5,5)
  temp.result <- LogoddsMetaAnalysis(logodds.icog,
                                     sigma.icog,
                                     logodds.onco,
                                     sigma.onco)
  total <- total+1
  result.sub[total,] <- c(temp.result[[1]],
                          as.vector(temp.result[[2]]))
}

save(result.sub,file=paste0("./genetic_correlation/ICOG/result/result.sub.meta",i1,".Rdata"))
