if(s==1){
  #theta_test <- c(0.05,0,0,0,0)
  theta_test <- c(0.07,0,0,0,0)
  #theta_test <- c(0.25,0,0,0,0)
}else if(s==2){
  #theta_test <- c(0,0.05,0,0,0)
  theta_test <- c(0,0.07,0,0,0)
}else{
  #theta_test <- c(c(0,0.05),rnorm(3,0,0.02))
  theta_test <- c(c(0,0.07),rnorm(3,0,0.02))
}

    temp.simu <- SimulateDataPower(theta_intercept,theta_test,theta_covar,n)
    y.pheno.mis <- temp.simu[[1]]
    G <- temp.simu[[2]]
    x_covar <- temp.simu[[3]]
model <- glm(y.pheno.mis[,1]~G+x_covar,family=binomial())    
summary(model)
