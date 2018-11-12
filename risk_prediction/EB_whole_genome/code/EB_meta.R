arg <- commandArgs(trailingOnly=T)
i1 <- as.numeric(arg[[1]])
library(bc2)
library(bcutility)
setwd('/spin1/users/zhangh24/breast_cancer_data_analysis/risk_prediction/')
load('./intrinsic_subtypes_whole_genome/ICOG/result/meta_result_shared_1p.Rdata')
meta_intrin <- meta_result_shared_1p
load("./standard_whole_genome/ICOG/result/meta_result_shared_1p.Rdata")
meta_stan <- meta_result_shared_1p

stan_result <- meta_stan[,c(11,12,17)]
intrin_result <- meta_intrin[,c(11:40,45)]
second.num <- 5
n <- nrow(meta_stan)
size <- 1000
start.end <- startend(n,size,i1)
start <- start.end[1]
end <- start.end[2]
eblogodds_sub <- matrix(0,end-start+1,second.num)
hetervar.sub <- rep(0,end-start+1)
temp = 1
for(j in start:end){
  
  print(j)
  intri_odds <-  as.numeric(intrin_result[j,1:second.num])
  intri_covar <-   matrix(
    as.numeric(intrin_result[j,(second.num+1):
              (second.num+second.num^2)]),
    second.num,second.num)
  hetervar.sub[temp] <- HeterVarianceEstimate(
    intri_odds,
    intri_covar)
  eblogodds_sub[temp,] <- ebestimate(intri_odds,
                                 intri_covar,
                                 stan_result[j,1],
                                 hetervar.sub[j])
  temp = temp+1
}

result_sub <- list(eblogodds_sub,
                   hetervar.sub)

save(result_sub,file=paste0("./EB_whole_genome/result/result_sub",i1,".Rdata"))

