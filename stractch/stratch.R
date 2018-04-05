load("/Users/zhangh24/GoogleDrive/breast_cancer_data_analysis/risk_prediction/two_stage_model/result/temp.Rdata")
log.odds <- temp[[1]]
sigma <- temp[[2]]
tree_function(log.odds,sigma)
library(mvtnorm)
tree_function <- function(log.odds,sigma){
  M <- length(log.odds)
  beta.mat <- matrix(0,M,M)
  beta.mat[1,] <- log.odds
  norm.mean <-  log.odds
  log.like <- dmvnorm(x= log.odds,mean= norm.mean,sigma=sigma,log=T)
  
  pairs <- combn(M,2)
  total.pair <- ncol(pairs)
  
  C.start <- diag(M)
  AIC <- rep(0,M)
  AIC[1] <- 2*M-2*log.like
  trans <- matrix(0,M,total.pair)
  distance <- rep(0,total.pair)
  for(i in 1:total.pair){
    trans[pairs[1,i],i] <- 1  
    trans[pairs[-1,i],i] <- -1  
    distance[i] <- abs(trans[,i]%*%log.odds/sqrt(t(trans[,i])%*%sigma%*%trans[,i]))
  }
  idx.min <- which.min(distance)
  min.pair <- pairs[,idx.min]
  new.col <- rowSums(C.start[,min.pair])
  C.update <- C.start[,-(min.pair),drop=F]
  C.update <- cbind(C.update,new.col)
  new.beta <- solve(t(C.update)%*%solve(sigma)%*%C.update)%*%t(C.update)%*%solve(sigma)%*%log.odds
  new.norm.mean <- C.update%*%new.beta
  beta.mat[2,] <- new.norm.mean
  new.sigma <- solve(t(C.update)%*%solve(sigma)%*%C.update)
  log.like <- dmvnorm(x= log.odds,mean= new.norm.mean,sigma=sigma,log=T)
  AIC[2] <- 2*length(new.beta)-2*log.like
  
  
  new.M <- M-1
  pairs <- combn(new.M,2)
  total.pair <- ncol(pairs)
  trans <- matrix(0,new.M,total.pair)
  distance <- rep(0,total.pair)
  for(i in 1:total.pair){
    trans[pairs[1,i],i] <- 1  
    trans[pairs[-1,i],i] <- -1  
    distance[i] <- abs(trans[,i]%*%new.beta/sqrt(t(trans[,i])%*%new.sigma%*%trans[,i]))
  }
  idx.min <- which.min(distance)
  min.pair <- pairs[,idx.min]
  new.col <- rowSums(C.update[,min.pair])
  C.update <- C.update[,-(min.pair),drop=F]
  C.update <- cbind(C.update,new.col)
  new.beta <- solve(t(C.update)%*%solve(sigma)%*%C.update)%*%t(C.update)%*%solve(sigma)%*%log.odds
  new.norm.mean <- C.update%*%new.beta
  beta.mat[3,] <- new.norm.mean
  new.sigma <- solve(t(C.update)%*%solve(sigma)%*%C.update)
  log.like <- dmvnorm(x= log.odds,mean= new.norm.mean,sigma=sigma,log=T)
  AIC[3] <- 2*length(new.beta)-2*log.like
  
  
  
  new.M <- M-2
  pairs <- combn(new.M,2)
  total.pair <- ncol(pairs)
  trans <- matrix(0,new.M,total.pair)
  distance <- rep(0,total.pair)
  for(i in 1:total.pair){
    trans[pairs[1,i],i] <- 1  
    trans[pairs[-1,i],i] <- -1  
    distance[i] <- abs(trans[,i]%*%new.beta/sqrt(t(trans[,i])%*%new.sigma%*%trans[,i]))
  }
  idx.min <- which.min(distance)
  min.pair <- pairs[,idx.min]
  new.col <- rowSums(C.update[,min.pair])
  C.update <- C.update[,-(min.pair),drop=F]
  C.update <- cbind(C.update,new.col)
  new.beta <- solve(t(C.update)%*%solve(sigma)%*%C.update)%*%t(C.update)%*%solve(sigma)%*%log.odds
  new.norm.mean <- C.update%*%new.beta
  beta.mat[4,] <- new.norm.mean
  new.sigma <- solve(t(C.update)%*%solve(sigma)%*%C.update)
  log.like <- dmvnorm(x= log.odds,mean= new.norm.mean,sigma=sigma,log=T)
  AIC[4] <- 2*length(new.beta)-2*log.like
  
  
  
  new.M <- M-3
  pairs <- combn(new.M,2)
  total.pair <- ncol(pairs)
  trans <- matrix(0,new.M,total.pair)
  distance <- rep(0,total.pair)
  for(i in 1:total.pair){
    trans[pairs[1,i],i] <- 1  
    trans[pairs[-1,i],i] <- -1  
    distance[i] <- abs(trans[,i]%*%new.beta/sqrt(t(trans[,i])%*%new.sigma%*%trans[,i]))
  }
  idx.min <- which.min(distance)
  min.pair <- pairs[,idx.min]
  new.col <- rowSums(C.update[,min.pair])
  C.update <- C.update[,-(min.pair),drop=F]
  C.update <- cbind(C.update,new.col)
  new.beta <- solve(t(C.update)%*%solve(sigma)%*%C.update)%*%t(C.update)%*%solve(sigma)%*%log.odds
  new.norm.mean <- C.update%*%new.beta
  beta.mat[5,] <- new.norm.mean
  new.sigma <- solve(t(C.update)%*%solve(sigma)%*%C.update)
  log.like <- dmvnorm(x= log.odds,mean= new.norm.mean,sigma=sigma,log=T)
  AIC[5] <- 2*length(new.beta)-2*log.like
  
  idx.AIC <- which.min(AIC)
  best.log.odds <- beta.mat[idx.AIC,]
return(best.log.odds)  
}

