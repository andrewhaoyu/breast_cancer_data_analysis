arg <- commandArgs(trailingOnly=T)
i1 <- as.numeric(arg[[1]])
z.design <- matrix(c(
  c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1),
  c(0,0,0,0,0,1,1,0,0,0,0,0,1,1,1,0,0,0,0,0,1,1,1),
  c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0),
  c(0,0,0,0,1,0,0,0,0,0,0,1,0,0,0,0,0,0,0,1,0,0,0),
  c(1,0,0,0,0,0,0,1,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0)
),ncol=5)

transfunction <- function(lododds,sigma){
  beta <- z.design%*%logodds
  beta.sigma <- z.design%*%sigma%*%t(z.design)
  p <- c(1,2,5,6,19)
  beta <- beta[p]
  beta.sigma <- beta.sigma[p,p]
  return(c(beta,as.vector(beta.sigma)))
}
setwd("/spin1/users/zhangh24/breast_cancer_data_analysis/")

load("./genetic_correlation/result/ICOG.result.clean.Rdata")
n <- nrow(ICOG.result.clean)
library(bc2)
size = 1000
start.end<- startend(n,size,i1)
start <- start.end[1]
end <- start.end[2]

result.sub <- matrix(0,(end-start+1),30)
for(i in start:end){
print(i)
    logodds <- as.numeric(as.vector(ICOG.result.clean[i,11:15]))
  sigma <- matrix(as.numeric(ICOG.result.clean[i,16:40]),5,5)  
  result.sub[i,] <- transfunction(logodds,sigma)
}

save(result.sub,paste0("./genetic_correlation/ICOG/result.sub",i1,".Rdata"))


