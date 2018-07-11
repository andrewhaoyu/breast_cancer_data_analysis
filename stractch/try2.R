n <- 10000
x <- rnorm(n,2,1)
y <- 2*x+rnorm(n,0,0.5)
x2 <- (x-mean(x))/sd(x)
model1 <- lm(y~x)
result.1 <- summary(model1)$coefficients[2,1]
model2 <- lm(y~x2)
result.2 <- summary(model2)$coefficients[2,1]
sd(x)*result.1
