




library(bc2)

load("/spin1/users/zhangh24/breast_cancer_data_analysis/risk_prediction/standard_analysis/result/log.odds.meta.Rdata")
load(paste0("/spin1/users/zhangh24/breast_cancer_data_analysis/risk_prediction/two_stage_model/result/z.standard.Rdata"))

finner_subtype_to_one <- function(logodds,sigma){
  log.odds.new <- sum(solve(sigma)%*%logodds)/sum(solve(sigma))
  sigma.new <- 1/sum(solve(sigma))
  return(list(log.odds.new,sigma.new))
  }
vectorMeta_analysis <- function(logodds,sigma){
  meta.logodds <- sum(1/sigma)^-1*sum(logodds/sigma)
  meta.var <- sum(1/sigma)^-1
  return(list(meta.logodds,meta.var))
}


heter.variance.estimate <- function(log.odds,sigma){
  M <- length(log.odds)
  result <- (sum((log.odds-mean(log.odds))^2)-sum(diag(sigma))+sum(sigma)/M)/(M-1)
  if(result <= 0){
    result <- 0
  }
  return(result)
}
z.grade <- as.matrix(cbind(1,z.standard[c(1,8,16),]))
z.design <- as.matrix(cbind(1,z.standard))
p.heter.add <- rep(0,205)
log.odds.add <- matrix(0,205,5)
sigma.add <- matrix(0,205,25)
log.odds.add.triple <- rep(0,205)
log.odds.add.triple.eb <- rep(0,205)

case.sigma <- rep(0,205)
baseline.sigma <- rep(0,205)


for(i1 in 1:205){
  print(i1)
  load(paste0("/spin1/users/zhangh24/breast_cancer_data_analysis/risk_prediction/two_stage_model/result/add.meta.result",i1,".Rdata"))
  log.odds.add[i1,] <- meta.result[[1]]
  sigma.add[i1,] <- as.vector(meta.result[[2]])
  temp <- vectorMeta_analysis(z.grade%*%meta.result[[1]],
                        diag(z.grade%*%meta.result[[2]]%*%t(z.grade)))
  log.odds.add.triple[i1] <- temp[[1]]
  
  case.sigma[i1] <- heter.variance.estimate(meta.result[[1]][2:5],
                          meta.result[[2]][2:5,2:5]) 
  baseline.sigma[i1] <- (meta.result[[1]][1]-log.odds.meta[i1])^2-meta.result[[2]][1,1]
  if(case.sigma[i1]<=0|baseline.sigma[i1]<=0){
    log.odds.add.triple.eb[i1] <- log.odds.meta[i1]
  }else{
    P <- length(meta.result[[1]])
    case.length <-length(meta.result[[1]])-1 
    K <- diag(c(baseline.sigma[i1],rep(case.sigma[i1] ,case.length )))
    log.odds.eb <- solve(solve(meta.result[[2]])+solve(K))%*%(solve(meta.result[[2]])%*%meta.result[[1]]+solve(K)%*%rep(log.odds.meta[i1],P))
    sigma.eb <- solve(solve(meta.result[[2]])+solve(K))
    log.odds.add.triple.eb[i1] <- vectorMeta_analysis(z.grade%*%log.odds.eb,
                                                      diag(z.grade%*%    sigma.eb%*%t(z.grade)))[[1]]
  }
  p.heter.add[i1] <- as.numeric(DisplaySecondStageTestResult(meta.result[[1]],meta.result[[2]])[12])
}

save(p.heter.add,file="/spin1/users/zhangh24/breast_cancer_data_analysis/risk_prediction/two_stage_model/result/p.heter.add.Rdata")
save(log.odds.add.triple,file="/spin1/users/zhangh24/breast_cancer_data_analysis/risk_prediction/two_stage_model/result/log.odds.add.triple.Rdata")
save(log.odds.add.triple.eb,file="/spin1/users/zhangh24/breast_cancer_data_analysis/risk_prediction/two_stage_model/result/log.odds.add.triple.eb.Rdata")
