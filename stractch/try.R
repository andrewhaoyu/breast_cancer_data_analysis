
try <- function(Sigma1,Sigma2){
  result <- c(0,0)
  result[1] <- (Sigma1[1,1]^-1+Sigma2[1,1]^-1)^-1
  result[2] <- (solve(solve(Sigma1)+solve(Sigma2)))[1,1]
  return(result)
}

n <- 1000
result <- matrix(0,n,2)
for(i in 1:n){
  A <- matrix(rnorm(4),2,2)
  Sigma1 <- A%*%t(A)
  
  B <- matrix(rnorm(4),2,2)
  Sigma2 <- B%*%t(B)
  result[i,] <- try(Sigma1,Sigma2)
}
