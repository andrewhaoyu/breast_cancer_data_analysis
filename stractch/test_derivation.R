lamda= 1
beta0 <- 4
sigma <- 2
n <- 10000
beta <- rnorm(n,beta0,sqrt(sigma))

betahat <- rnorm(n,beta,sqrt(lamda))
var(betahat)
mean(betahat)
var(betahat)-1
