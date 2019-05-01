#Goal: Generate phenotypes data for
#simulate phenotypes data for AFR, EUR, AMR
#MAF based on 1000KG
#sample size EUR n =120000
#sample size AFR n = 18000
#sample size AMR n = 18000
#heritability for EUR 0.8
#heritability for AFR 0.8
#heritability for AMR 0.8
#5000 causal SNPs for each population
#4000 shared causal SNPs
#Genetic correlation for the shared SNPs is 0.6
#1000 SNPs for each as independent causal
setwd('/spin1/users/zhangh24/breast_cancer_data_analysis')
load(paste0("./multi_ethnic/result/pruned_geno/geno_",1))

n.shared <- 4000
n.nonshare <- 1000
n.cau <- n.shared+n.nonshare
her.EUR <- 0.8/n.cau
her.AFR <- 0.8/n.cau
her.AMR <- 0.8/n.cau
gr <- 0.6
her.cov <- gr*sqrt(her.EUR*her.AFR)
library(mvtnorm)


vec.sigma <- rep(her.cov,9)
vec.sigma[c(1,5,9)] <- c(her.EUR,her.AFR,her.AMR)
Sigma <- matrix(vec.sigma,3,3)

#since the SNP MAF data is permutated
#for EUR and AFR, we will use 1:4000 as shared SNPs
#use 4000:5000 as nonshared SNPs for EUR
#5000:6000 as nonshared SNPs for AFR
#6000:7000 as nonshared SNPs for AMR

#shared beta coefcients
set.seed(666)
beta_shared <- rmvnorm(n.shared,c(0,0,0),
                       sigma=Sigma)
#non shared beta coefcients
all_nonshare <- 3*n.nonshare
beta_nonshared <- matrix(0,all_nonshare,3)
beta_nonshared[1:n.nonshare,1] <- rnorm(n.nonshare,0,sd=sqrt(her.EUR/n.cau))
beta_nonshared[(n.nonshare+1):(n.nonshare*2),2] <- rnorm(n.nonshare,0,sd=sqrt(her.AFR/n.cau))
beta_nonshared[(n.nonshare*2+1):(all_nonshare),3] <- rnorm(n.nonshare,0,sd=sqrt(her.AMR/n.cau))
beta <- rbind(beta_shared,beta_nonshared)

#######save effect size for the causal SNPs
save(beta,file = "./multi_ethnic/result/pruned_geno/effect_size.Rdata")









i1 = 1
set.seed(i1)
n.EUR <- nrow(genotype[[1]])
n.AFR <- nrow(genotype[[2]])
n.AMR <- nrow(genotype[[3]])

y_EUR <- genotype[[1]]%*%beta[,1]+rnorm(n.EUR,0,sd=sqrt(1-her.EUR))
y_AFR <- genotype[[2]]%*%beta[,2]+rnorm(n.AFR,0,sd=sqrt(1-her.AFR))
y_AMR <- genotype[[3]]%*%beta[,3]+rnorm(n.AMR,0,sd=sqrt(1-her.AMR))

y <- list(y_EUR,y_AFR,y_AMR)
save(y,file = paste0("./multi_ethnic/result/y_",i1))











y_EUR_train <- y_EUR[1:n.train.EUR]
y_AFR_train <- y_AFR[1:n.train.AFR]
y_EUR_test <- y_EUR[(n.train.EUR+1):(n.train.EUR+n.test.EUR)]
y_AFR_test <- y_AFR[(n.train.AFR+1):(n.train.AFR+n.test.AFR)]


####regression on the training and testing dataset
####record the summary level statistics
####first on the effect SNPs
####simulate the other SNPs and record the summary level statistics
n.snp <- nrow(pruned.snp.permu)
beta_summary_train <- matrix(0,n.snp,6)
beta_summary_test <- matrix(0,n.snp,6)
colnames(beta_summary_train) <- c("beta_EUR","sd_EUR","p_EUR",
                                  "beta_AFR","sd_AFR","p_AFR")
colnames(beta_summary_test) <- colnames(beta_summary_train)


FitLinearmodel <- function(y,x){
  model <- fastLm(X=cbind(1,x),y=y)
  if(is.na(coef(model)[2])){
    result <- c(0,1,1)
  }else{
    result <- coef(summary(model))[2,c(1,2,4)]  
  }
  
  return(result)
}
temp=1
for(i in 1:n.cau){
  if(i%%100==0){
    print(i)
  }
  beta_summary_train[temp,1:3] <- FitLinearmodel(y_EUR_train,G_EUR_effect[1:n.train.EUR,i])
  beta_summary_test[temp,1:3] <-  FitLinearmodel(y_EUR_test,G_EUR_effect[(n.train.EUR+1):(n.train.EUR+n.test.EUR),i])
  beta_summary_train[temp,4:6] <- FitLinearmodel(y_AFR_train,G_AFR_effect[1:n.train.AFR,i])
  beta_summary_test[temp,4:6] <-  FitLinearmodel(y_AFR_test,G_AFR_effect[(n.train.AFR+1):(n.train.AFR+n.test.AFR),i])
  temp = temp+1
}

for(i in 1:(n.snp-n.cau)){
  if(i%%100==0){
    print(i)
  }
  G_EUR <- rbinom(n.EUR,2,MAF.EUR[n.cau+i])
  G_AFR <- rbinom(n.AFR,2,MAF.AFR[n.cau+i])
  beta_summary_train[temp,1:3] <- FitLinearmodel(y_EUR_train,G_EUR[1:n.train.EUR])
  beta_summary_test[temp,1:3] <-  FitLinearmodel(y_EUR_test,G_EUR[(n.train.EUR+1):(n.train.EUR+n.test.EUR)])
  beta_summary_train[temp,4:6] <- FitLinearmodel(y_AFR_train,G_AFR[1:n.train.AFR])
  beta_summary_test[temp,4:6] <-  FitLinearmodel(y_AFR_test,G_AFR[(n.train.AFR+1):(n.train.AFR+n.test.AFR)])
  temp = temp+1
}


temp.result <- list(beta_summary_train,
                    beta_summary_test,
)














