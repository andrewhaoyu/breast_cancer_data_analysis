n <- 10000

x1 = runif(n)
y1 = runif(n)

x2 = runif(n)
y2 = runif(n)

x3 = runif(n)
y3 = runif(n)

l1 = x1^2+y1^2
l2 = x2^2+y2^2
l3 = x3^2+y3^2
idx <- which(l1<=1&l2<=1&l3<=1)

GetTheta = function(x1,y1,x2,y2){
  acos((x1*x2+y1*y2)/(sqrt(x1^2+y1^2)*sqrt(x2^2+y2^2)))
}

GetTheta(x1[1],y1[1],x2[1],y2[1])
