#goal:merge TN interaction model results from multiple imputation
setwd('/data/zhangh24/breast_cancer_data_analysis/')

coef.all <- list()
sigma.all <- list()
n.imp <- 100
n.sub <- 23

for(i1 in 1:n.imp){
  load(paste0("./risk_factor/result/TN_interaction_imp",i1,".Rdata"))  
  coef.all[[i1]] <- model[[1]]
  sigma.all[[i1]] <- model[[2]]
}
result.temp <- RubinPool(coef.all,sigma.all,n.imp)
coef.temp <- matrix(result.temp[[1]][(n.sub+1):(n.sub+length(model[[3]]))],nrow=6)
colnames(coef.temp) <- colnames(model[[3]])
rownames(coef.temp) <- rownames(model[[3]])
rownames(coef.temp)[1] <- "baseline effect"
sigma.temp <- result.temp[[2]][(n.sub+1):(n.sub+length(model[[3]])),
                               (n.sub+1):(n.sub+length(model[[3]]))]
result <- list(coef.temp,sigma.temp)
save(result,file = paste0("./risk_factor/result/TN_interaction_imp_merge.Rdata"))
RubinPool <- function(coef.all,sigma.all,n.imp){
  coef.result.sum <- coef.all[[1]]
  sigma.result.sum <- sigma.all[[1]]
  for(i in 2:n.imp){
    coef.result.sum <- coef.result.sum+coef.all[[i]]
    sigma.result.sum <- sigma.result.sum+sigma.all[[i]]
  }
  coef.result <- coef.result.sum/n.imp
  ##
  U_bar <- sigma.result.sum/n.imp
  #take vector by row, just to match the sigma variance ordering
  Q_bar = as.vector(t(coef.result))
  Q_temp = as.vector(t(coef.all[[1]]))
  Q_sum= (Q_temp-Q_bar)%*%t(Q_temp-Q_bar)
  
  for(i in 2:n.imp){
    Q_temp = as.vector(t(coef.all[[i]]))
    Q_sum = Q_sum +(Q_temp-Q_bar)%*%t(Q_temp-Q_bar)
  }
  B = Q_sum/(n.imp-1)
  sigma.result <-   U_bar+B+B/n.imp
  return(list(coef.result,sigma.result))
}
