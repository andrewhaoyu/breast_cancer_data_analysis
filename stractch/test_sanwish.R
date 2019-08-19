#x1 and x2 are correlated
library(MASS)

p = 1000
beta_est <- rep(0,p)

sigma1 <- 1
sigma2 <- 1
cov1 <- 0
beta1 <- 0.3
beta2 <- 1.5
error = 1 
mean1 = 0
mean2 = 0
betam = 1
alpha1 = 0.2
alpha2 = 0.3
M_est = rep(0,p)
for(i in 1:p){
  Sigma <- matrix(c(sigma1,cov1,cov1,sigma2),2,2)
  n <- 10000
  x <- mvrnorm(n,mu = c(mean1,mean2),Sigma = Sigma)
  x1 = x[,1]
  x2 = x[,2]
  M = alpha1*x1+alpha2*x2+rnorm(n)
  y = x1*beta1+x2*beta2+M*betam+rnorm(n)
  
  model <- lm(y~x1+M-1)
  beta_est[i] <- coefficients(model)[1]
  M_est[i] <- coefficients(model)[2]
}
mean(beta_est)
mean(M_est)


part2 = matrix(0,2,2)
part2[1,1] = crossprod(x1)
part2[2,1]= part2[1,2] = crossprod(x1,M)
part2[2,2] = crossprod(M)



part1 = matrix(0,2,2)
part1[1,1] = crossprod(x1)
part1[2,1]= part1[1,2] = alpha1*crossprod(x1)
part1[2,2] = alpha1^2*crossprod(x1)+alpha2^2*crossprod(x2)+error*n

c(beta1,betam)+solve(part1)[,2]*alpha2*beta2*crossprod(x2)





beta1+sigma1^-1*cov1*beta2

x2Gx1 = mean2+cov1*sigma1^-1*(x1-mean1)
meanest = beta1+crossprod(x1)^-1*crossprod(x1,x2Gx1)*beta2
meanest2 = beta1+crossprod(x1)^-1*crossprod(x1,x2)*beta2
meanest
meanest2
varx2givenx1 = sigma2-cov1*sigma1^-1*cov1

#est3 = error*crossprod(x1)^-1+cross







p = 1000
beta_est <- rep(0,p)

sigma1 <- 1
sigma2 <- 1
cov1 <- 0
beta1 <- 0.3
beta2 <- 1.5
error = 1 
mean1 = 0
mean2 = 0
betam = 1
alpha1 = 0.2
alpha2 = 0.3
M_est = rep(0,p)
for(i in 1:p){
  Sigma <- matrix(c(sigma1,cov1,cov1,sigma2),2,2)
  n <- 10000
  x <- mvrnorm(n,mu = c(mean1,mean2),Sigma = Sigma)
  x1 = x[,1]
  x2 = x[,2]
  y = x1*beta1+x2*beta2+rnorm(n)
  
  model <- lm(y~x1-1)
  beta_est[i] <- coefficients(model)[1]
}
mean(beta_est)
est1 = var(beta_est)
est2 = crossprod(x1)^-1*error+crossprod(x1)^-1*sum(x1*x2*x2*x1)*beta2^2*crossprod(x1)^-1
est1/est2

sum(x1*x2*x2*x1)/crossprod(x1,x2)^2

part2 = matrix(0,2,2)
part2[1,1] = crossprod(x1)
part2[2,1]= part2[1,2] = crossprod(x1,M)
part2[2,2] = crossprod(M)



part1 = matrix(0,2,2)
part1[1,1] = crossprod(x1)
part1[2,1]= part1[1,2] = alpha1*crossprod(x1)
part1[2,2] = alpha1^2*crossprod(x1)+alpha2^2*crossprod(x2)+error*n

c(beta1,betam)+solve(part1)[,2]*alpha2*beta2*crossprod(x2)





beta1+sigma1^-1*cov1*beta2

x2Gx1 = mean2+cov1*sigma1^-1*(x1-mean1)
meanest = beta1+crossprod(x1)^-1*crossprod(x1,x2Gx1)*beta2
meanest2 = beta1+crossprod(x1)^-1*crossprod(x1,x2)*beta2
meanest
meanest2
varx2givenx1 = sigma2-cov1*sigma1^-1*cov1






n = 5000
G = rnorm(n)
M = G+rnorm(n)
Y = 2*M + rnorm(n)
lm(Y~M)
lm(Y~G+M)
lm(Y~M+G)

try_est = rep(0,n)
for(i in 1:p)
G= rnorm(n)
M = 2*G+rnorm(n)

try = t(M)%*%G%*%crossprod(G)^-1%*%t(G)%*%M
try2 = crossprod(M)/n
try3 = 
try/try2


total = 0
total2 = 0
p = 10000
for(i in 1:p){
  n = 100
  G = rnorm(n)
  M = 2*G+rnorm(n)
  total = total + t(M)%*%(G%*%crossprod(G)^-1%*%t(G))%*%M
  #total = total + t(M)%*%G%*%crossprod(G)^-1%*%t(G)%*%M
  total2 = total2 + t(M)%*%(diag(n)/n)%*%M
}
total/p
total2/p
sum((total/p-total2/p)^2)
