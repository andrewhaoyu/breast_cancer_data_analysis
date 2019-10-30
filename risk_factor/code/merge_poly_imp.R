#goal:merge polytomous model results from multiple imputation


setwd('/data/zhangh24/breast_cancer_data_analysis/')

coef.all <- list()
sigma.all <- list()
n.imp <- 100
for(i1 in 1:n.imp){
  load(paste0("./risk_factor/result/poly_imp",i1,".Rdata"))  
  coef.all[[i1]] <- result[[1]]
  sigma.all[[i1]] <- result[[2]]
}
result <- RubinPool(coef.all,sigma.all,n.imp)
save(result,file = paste0("./risk_factor/result/poly_imp_merge.Rdata"))
#####Rubin's formula
#####Q is the estimate for the lth imputation
#####U is the sigma for the lth imputation
#####Q_bar = mean(Q) would be the estimate
#####U_bar = mean(U) would be the within imputation variable
#####B = 1/(m-1)*sum(Q_l-Q_bar)*t(Q_l-Q_bar) would be the cross imputation variable
#####final variable = U_bar+B+B/m
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
