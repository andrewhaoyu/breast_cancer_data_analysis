prior.sigma <- function(log.odds,sigma){
  p <- ncol(sigma)
  n <- nrow(log.odds)
  mean.log.odds <- apply(log.odds,1,mean)
  n.eff <- 1/sigma
  n.eff.sum <- apply(n.eff,1,sum)
  prior.sigma.result <- rep(0,205)
  
  for(i in 1:n){
    temp <- 0
    for(j in 1:p){
      temp <- (n.eff[i,j]/n.eff.sum[i])*(log.odds[i,j]-mean.log.odds[i])^2
      prior.sigma.result[i] <- temp+prior.sigma.result[i]
    }
  }
  prior.sigma.result <- prior.sigma.result*p/(p-1)
  return(prior.sigma.result)
}




library(bc2)
p.heter.add <- rep(0,205)
log.odds.add <- matrix(0,205,5)
sigma.add <- matrix(0,205,25)
log.odds.add.tri <- rep(0,205)
sigma.add.tri <- rep(0,205)
sigma.add.tri <- rep(0,205)
sigma.first <- matrix(0,205,23)
log.odds.add.first <- matrix(0,205,23)





load(paste0("/spin1/users/zhangh24/breast_cancer_data_analysis/risk_prediction/two_stage_model/result/z.standard.Rdata"))

load(paste0("/spin1/users/zhangh24/breast_cancer_data_analysis/risk_prediction/two_stage_model/result/grade.ratio.Rdata"))
z.grade <- as.matrix(cbind(1,z.standard[c(1,8,16),]))
z.design <- as.matrix(cbind(1,z.standard))
for(i1 in 1:205){
  print(i1)
  load(paste0("/spin1/users/zhangh24/breast_cancer_data_analysis/risk_prediction/two_stage_model/result/add.meta.result",i1,".Rdata"))
  log.odds.add[i1,] <- meta.result[[1]]
  sigma.add[i1,] <- as.vector(meta.result[[2]])
  log.odds.add.tri[i1] <- t(z.grade%*%meta.result[[1]])%*%grade.ratio
  sigma.add.tri[i1] <- grade.ratio%*%z.grade%*%meta.result[[2]]%*%t(z.grade)%*%grade.ratio
  log.odds.add.first[i1,] <- z.design%*%meta.result[[1]]
  sigma.first[i1,] <- diag(z.design%*%meta.result[[2]]%*%t(z.design))
  p.heter.add[i1] <- as.numeric(DisplaySecondStageTestResult(meta.result[[1]],meta.result[[2]])[12])
}

prior.sigma.add <- prior.sigma(log.odds.add.first, sigma.first)

save(p.heter.add,file="/spin1/users/zhangh24/breast_cancer_data_analysis/risk_prediction/two_stage_model/result/p.heter.add.Rdata")
save(log.odds.add.tri,file="/spin1/users/zhangh24/breast_cancer_data_analysis/risk_prediction/two_stage_model/result/log.odds.add.tri.Rdata")
save(prior.sigma.add ,file="/spin1/users/zhangh24/breast_cancer_data_analysis/risk_prediction/two_stage_model/result/prior.sigma.add.Rdata")
save(sigma.add.tri,file="/spin1/users/zhangh24/breast_cancer_data_analysis/risk_prediction/two_stage_model/result/sigma.add.tri.Rdata")
