#-------------------------------------------------------------------
# Update Date: 01/20/2019
# Create Date: 01/20/2019
# Goal: implement fisher score and gradient descent of logisitic regression
# Author: Haoyu Zhang
#-------------------------------------------------------------------
#---------------------------------------

n <- 100
x1 = rnorm(n)
x2 = rnorm(n)

logit_inv <- function(x){
  exp(x)/(1+exp(x))
}

p = logit_inv(1+x1*2+x2*1)

y = rbinom(n,1,p)

library(ggplot2)
data <- data.frame(y=y,x1=x1,x2=x2)

 
model <- glm(y~x1+x2,family=binomial)
beta <- summary(model)$coefficient[,1]

ab_fun <- function(beta){
  a = -beta[1]/beta[2]
  b = -beta[2]/beta[3]
  return(c(a,b))
}
ab <- ab_fun(beta)
ggplot(data,aes(x1,x2)) + geom_point(aes(color = y))+
  geom_abline(intercept= ab[1],
              slope = ab[2])
y_predict <- ifelse(logit_inv(predict(model))>0.5,1,0)
error <- mean((y-y_predict)^2)






theta_old = c(0,1,1)
x  = cbind(1,x1,x2)
#############gradient descent
n.up = 100
alpha = 0.1
for(i in 1:n.up){
  print(i)
  ab <- ab_fun(theta_old)
  p1 = ggplot(data,aes(x1,x2)) + geom_point(aes(color = y))+
    geom_abline(intercept= ab[1],
                slope = ab[2])
  #print(p1)
  theta_new = theta_old + alpha*t(x)%*%(y-logit_inv(x%*%theta_old))
  error = max(abs(theta_new-theta_old))
  if(error<1/(10*n)){
    break
  }
  theta_old = theta_new
  print(theta_new)
  Sys.sleep(1)  
}


theta_old = c(0,1,1)
x  = cbind(1,x1,x2)
#############fisher scoring
n.up = 100
alpha = 0.1
#save matrix x_i * x_i^t
#save as a list
x_list <- list()
for(i in 1:n){
  x_list[[i]] <- crossprod(t(x[i,]),x[i,])
}
infor_fun <- function(x_list,w_vec,p){
  infor <- matrix(0,p,p)
  for(i in 1:length(w_vec)){
  infor <- infor + w_vec[i]*x_list[[i]]    
  }
  return(infor)
}
p <- length(theta_old)
for(i in 1:n.up){
  print(i)
  ab <- ab_fun(theta_old)
  p1 = ggplot(data,aes(x1,x2)) + geom_point(aes(color = y))+
    geom_abline(intercept= ab[1],
                slope = ab[2])
  mu = logit_inv(x%*%theta_old)
  w = mu*(1-mu)
  infor <- infor_fun(x_list,w,p)
    
  theta_new = theta_old + solve(infor)%*%t(x)%*%(y-mu)
  error = max(abs(theta_new-theta_old))
  if(error<1/(10*n)){
    break
  }
  theta_old = theta_new
  print(theta_new)
  Sys.sleep(1)  
}

