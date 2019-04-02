#goal: test mice
library(mice)
library(VIM)
head(nhanes)
md.pattern(nhanes)

p <- md.pairs(nhanes)
p
pbox(nhanes,pos=1,int=F,cex=0.7)

imp <- mice(nhanes)
imp$imp$chl
head(complete(imp))
head(complete(imp,2))
summary(lm(chl~age+hyp,data=nhanes))
confint(lm(chl~age+hyp,data=nhanes))
fit <- with(imp,lm(chl~age+hyp))
summary(pool(fit))





n <- 100
x1 <- rnorm(n)
x2 <- 0.5*x1+rnorm(n)
cor(x1,x2)
y <- 1*x1+0.5*x2+rnorm(n)
summary(lm(y~x1+x2))


idx1 <- sample(c(1:n),0.2*n)
idx2 <- sample(c(1:n),0.2*n)
x1.new <- x1
x1.new[idx1] <- NA
x2.new <- x2
x2.new[idx2] <- NA
data.new <- data.frame(y,x1.new,x2.new)

model2 <- lm(y~x1.new+x2.new)
summary(model2)
md.pattern(data.new)
md.pairs(data.new)

imp <- mice(data.new,m=10)
head(complete(imp))
fit <- with(imp,lm(y~x1.new+x2.new))
summary(pool(fit))


library(mice)
imp <- mice(airquality,method="mean",m=1,maxit=1)
fit <- lm(Ozone~Solar.R,data=airquality)
pred <- predict(fit,newdata=ic(airquality))


imp <- mice(airquality,seed=1,print=FALSE)
fit <- with(imp,lm(Ozone~Wind+Temp+Solar.R))
tab <- summary(pool(fit))
tab




logit_trans <- function(x){exp(x)/(1+exp(x))}




library(MASS)
n <- 1000
y <- mvrnorm(n,mu=c(0,0),Sigma= matrix(c(1,0.5,0.5,1),2,2))
y1 <- y[,1]
y2 <- y[,2]
alpha1 <- 0
alpha2 <- 1
alpha3 <- 0
p <- alpha1+logit_trans(y1)*alpha2+logit_trans(y2)*alpha3
r <- 1- rbinom(n,1,p)
y_new <- y1[r,]



