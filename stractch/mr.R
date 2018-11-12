n <- 50000
Z <- rnorm(n)
Z2 <- rnorm(n)
X <- 3*Z2+0.5*Z+rnorm(n)
X2 <- 0.3*Z2+rnorm(n)
Y = 0.5*X+3*X2+rnorm(n)

model <- lm(Y~X)
summary(model)
model1 <- lm(Y~Z)
summary(model1)
model2 <- lm(X~Z)
summary(model2)
coef(model1)[2]/coef(model2)[2]

model4 <- lm(Y~X+X2)
summary(model4)
