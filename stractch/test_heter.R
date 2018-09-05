M <- 5
A <- matrix(rnorm(M^2,mean=0,sd=0.2),M,M)
Sigma <- A%*%t(A)

heter.sigma <- 0.3
b <- sqrt(heter.sigma/2)
M <- 5
n <- 10000
beta0 = 1
Q <- matrix(runif(M*M),M,M)
Q <- Q%*%t(Q)
R <- cov2cor(Q)

beta <- mvrnorm(n=n,mu=rep(beta0,M),Sigma=(R*heter.sigma))
#beta <- matrix(rnorm(n*M,mean=beta0,sd=sqrt(heter.sigma)),n,M)
#beta <- matrix(rlaplace(n*M,m=beta0,s=b),n,M)
library(MASS)
#library(rmutil)

beta_hat <- beta
heter.sigma.est <- rep(0,n)
for(i in 1:n){
  beta_hat[i,] <- mvrnorm(1,beta[i,],Sigma) 
  heter.sigma.est[i] <- heter.variance.estimate(beta_hat[i,],Sigma,R)
}
mean(heter.sigma.est)
heter.variance.estimate <- function(log.odds,sigma,R){
  M <- length(log.odds)
  result <- (sum((log.odds-mean(log.odds))^2)-sum(diag(sigma))+sum(sigma)/M)/(M-sum(R)/M)
  #result <- (sum((log.odds-mean(log.odds))^2)-(M+1)*sum(diag(sigma))/M)/(M-1)
  #if(result <= 0){
  #  result <- 0
 # }
  return(result)
}
i <- 3
betahat <- beta_hat[i,]
hetersigma <- heter.sigma.est[i]
b <- sqrt(hetersigma/2)
data <- list(M=M,beta0=beta0,betahat=betahat,Sigma=Sigma,
             b=b)
library(rstan)
stan.model1 <- '
data{
//define data
int<lower=1> M;
real beta0;
vector[M] betahat;
matrix[M,M] Sigma;
real<lower=0> hetersigma;
}
parameters{
vector[M] beta;
}
model{
//prior
for(i in 1:M){
beta[i] ~ normal(beta0,hetersigma);
}
//data
betahat ~ multi_normal(beta,Sigma);
}
'
stan.model2 <- '
data{
//define data
int<lower=1> M;
real beta0;
vector[M] betahat;
matrix[M,M] Sigma;
real<lower=0> b;
}
parameters{
vector[M] beta;
}
model{
//prior
for(i in 1:M){
beta[i] ~ double_exponential(beta0,b);
}
//data
betahat ~ multi_normal(beta,Sigma);
}
'

smodel <- stan_model(model_code = stan.model2)
fit1 <- sampling(smodel,
                 data=data,
                 warmup=5000,
                 iter=10000,
                 control = list(adapt_delta=0.95),
                 chains=4)
traceplot(fit1,pars=c('beta'))
pred.y <- extract(fit1)
sum((colMeans(pred.y[[1]])-beta[i,])^2)
sum((beta_hat[i,]-beta[i,])^2)
ebestimate(betahat,Sigma,beta0,hetersigma)
data <- list(M=M,beta0=beta0,betahat=betahat,Sigma=Sigma,
             hetersigma = hetersigma)
ebestimate <- function(logodds.subtype,
                       sigma.subtype,
                       logodds.standard,
                       prior.sigma
){
  M <- length(logodds.subtype)
  if(prior.sigma==0){
    return(rep(logodds.standard,M))
  }else{
    result <- solve(solve(sigma.subtype)+(1/prior.sigma)*diag(M))%*%(solve(sigma.subtype)%*%logodds.subtype+
                                                                       (1/prior.sigma)*rep(logodds.standard,M))
    return(result)
  }
}