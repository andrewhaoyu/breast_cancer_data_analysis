##################Generate Simulation Data###########################
#This is an example with 10000 samples;
#There are 3 covaraites in total;p_col means the number of covariates;
#Three second stage categories;
rm(list=ls())
set.seed(1234)
###########a,b,c represent three different tumor characteristics
###########under breast cancer example, they could be ER,PR,HER2,these three are binary tumor characteristics
a <- c(0,1)
b <- c(0,1)
c <- c(0,1)
p_col <- 3
z <- as.matrix(expand.grid(a,b,c)) # orig
#z <- as.matrix(expand.grid(a,b))


#this z_design matrix is the second stage matrix for a single covariate
z_design <- cbind(1,z)
M <- nrow(z_design)
#z_all is the second stage design matrix for all the covariates
z_all <- NULL
z_all_temp <- NULL
for(i in 1:M){
  z_all_temp <- rbind(z_all_temp,kronecker(diag(p_col),t(z_design[i,])))
}

z_all <- matrix(0,nrow = M*(p_col+1),ncol= M+p_col*ncol(z_design))
for(i in 1:M){
  z_all[1+(i-1)*(p_col+1),i] <- 1
}
for(i in 1:M){
  z_all[(2+(i-1)*(p_col+1)):(i*(p_col+1)),(M+1):ncol(z_all)] <- 
    z_all_temp[(1+(i-1)*p_col):(i*p_col),]
}
# for(i in 1:(M)){
#   temp <- rep(0,ncol(z_all))
#   temp[i] <- 1
#   z_all[1+(i-1)*(p_col+1),] = temp
# }
#K is the total number of second stage paramters for a single covariate
K <- ncol(z_design)

# z <- kronecker(diag(2),z)

theta_intercept <- rep(0.2,M)
theta_test <- rep(0,K)
theta_covar <- rep(c(1:K)/10,p_col-1)

#this theta is the true value
theta <- c(theta_intercept,theta_test,theta_covar)

#this is the true beta
beta <- z_all%*%theta
beta <- matrix(beta,nrow=p_col+1)
#alpha <- c(0,rep(1,length(beta)-1))

n <- 10000
x <-  matrix(rnorm(p_col*n),nrow = n)
##x_test represent the variale you are interested (SNP)
##x_covar represent the other covariates (Principle components,age,etc.)
x_test <- x[,1]
x_covar <- x[,2:ncol(x)]
x_all <- cbind(1,x) ##adding the intercept into the model

predictor <- x_all%*%(beta)
predictor <- cbind(0,predictor)
#predictor <- sweep(predictor,2,alpha,"+")
p <- exp(predictor)
sp <- rowSums(p)
#this standarize the probability sum into 1,since logit model:pr(D=1|predictor)=exp(predictor)/(1+exp(predictor))
p <- sweep(p,1,sp,"/")
y <- t(apply(p,1,function(x){rmultinom(1,1,x)}))

y <- y[,-1] 


#########this y matrix is a binary matrix
#########each row represent a person, each column represent a subtype
#########with ER,PR,HER2 three different tumor characteristics
#########there are 8 possible subtypes
#########there order are defined based on the z matrix
#########the first row of z matrix is 0,0,0, this represent ER-PR-HER2-
#########the second row of z matrix is 1,0,0, this represent ER-PR+HER2-


##transform y to the tumor characteristics type
y.pheno <- matrix(0,n,3)
for(i in 1:n){
  print(i)
  idx <- which(y[i,]==1)
  ##########length(idx)==0 represent control
  if(length(idx)==0){
    y.pheno[i,]=NA
  }else{
    y.pheno[i,] <- z[idx,]
  }
}


x <- cbind(x_test,x_covar)
colnames(x) <- c("SNP","PC1","PC2")

y.case.control <- rowSums(y)
y.pheno <- cbind(y.case.control,y.pheno)
colnames(y.pheno) <- c("case_control_status","ER","PR","HER2")

##################Finish Generating Simulation Data###########################


library(devtools)
##########this is my package implementing the two-stage model methods
install_github("andrewhaoyu/bc2")
library(bc2)


#############under this case, you need to put into twoparameters
#############first is y.pheno matrix, which contain the casecontrol status and the tumor characteristics
#############second is the link between the second stage and first stage parameters
#############here I will use additive model

result <- TwoStageModel(y.pheno,additive=x)

#######this is the two-stage model results
#######result[[1]] are the second stage parameters estimate
#######interect have 8 second stage parameters(case-control log odds for each subtype)
#######SNP,PC1 and PC2 have each have 4 second stage parameters (case control log odds ratio for a reference subtype and the casecase log odds ratio for ER,PR,HER2)

#######result[[2]] is the covariance matrix for the 20 second stage parameters
#######result[[3]] is the matrix form of the second stage parameter(we took out intersect to make it easier for people to look)
#######result[[4]]: transform the log odds ratio to odds ratio. We also showed the 95% CI and p-value
#######result[[5]]: the global association test p-value: test whether this snp will have an effect with the any subtypes
#######resutlt[[6]]: corresponding log odds ratio for 8 different subtypes
#######result[[7]]: the odds ratio for 8 different subtypes with their 95% CI and p-value
#######result[[8]]: likelihood
#######result[[9]]:AIC




######now we simulate the case where we have some missing tumor characteristics
######let's set ER,PR,HER2 each have 20% missing, we use 888 to represent missing
idx.case <- which(y.pheno[,1]==1)
for(i in 2:4){
  n.case <- length(idx.case)
  idx.mis <- rbinom(n.case,1,0.2)
  y.pheno[idx.case,i][idx.mis==1] <- 888 
}
result <- TwoStageModel(y.pheno,additive=x,missingTumorIndicator = 888)
