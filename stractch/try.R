p <- 0.5 
n <- 600000
result_m <- rep(0,n)
result_f <- rep(0,n)

for(i in 1:n){
  count <- 1
  temp <- rbinom(1,1,0.5)
  if(temp==1){
    result_m[i] = 1
    result_f[i] = 0
  }
  while(temp!=1){
    temp = rbinom(1,1,0.5)
    count = count + 1
  }
  result_m[i] = 1
  result_f[i] = count-1
  
  
}
sum(result_m)/sum(result_f)
