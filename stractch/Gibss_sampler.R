#-------------------------------------------------------------------
# Update Date: 01/29/2019
# Create Date: 01/29/2019
# Goal: implement Gibbs sampler
# Author: Haoyu Zhang
#-------------------------------------------------------------------
#---------------------------------------
GenA_G_B <- function(B){
  p1= 0.8
  p2 = 0.4
  if(B==0){
    A = rbinom(1,1,p1)
  }else if(B==1){
    A = rbinom(1,1,p2)
  }
  return(A)
}
GenB_G_A <- function(A){
  p1= 3/4
  p2 = 1/3
  if(A==0){
    B = rbinom(1,1,p1)
  }else if(A==1){
    B = rbinom(1,1,p2)
  }
  return(B)
}

n <- 10000
A_r <- rep(0,n)
B_r <- rep(0,n)

for(i in 1:(n-1)){
  A_r[i] <- GenA_G_B(B_r[i])
  B_r[i+1] <- GenB_G_A(A_r[i])
}
result <- cbind(A_r,B_r)
result <- result[(n/2):(n-1),]
#hist(result)
table(result[,1],result[,2])/(nrow(result))

