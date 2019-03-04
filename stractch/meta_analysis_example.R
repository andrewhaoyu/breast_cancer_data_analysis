########n studies with incidence rate and standard error
########beta the vector of incidence rate
########se the vector of standard error
MetaFunc <- function(beta,se){
  n <- length(beta)
  se_meta <- sqrt(sum(1/se^2)^-1)
  beta_meta <- se_meta^2*(sum(beta/se^2))
  return(c(beta_meta,se_meta))
}

######example
######1000 studies
n <- 1000
beta <- rnorm(n)
se <- runif(n,0,0.5)

result <- MetaFunc(beta,se)

beta_meta <- result[1]
se_meta <- result[2]
