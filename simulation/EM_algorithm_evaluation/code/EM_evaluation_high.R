###########simulate data with six tumor characteristics containing five binary tumor characteristics and one oridinal tumor characteristics.
###########one genotype with MAF 0.25 is simulated
###########one covariate with rnorm(n) simulated

library(nnet)
SimulateData <- function(theta_intercept,theta_test,theta_covar,n){
  a <- c(0,1)
  b <- c(0,1)
  c <- c(0,1)
  d <- c(1,2,3)
  e <- c(0,1)
  f <- c(0,1)
  ##number of the other covarites
  p_col <- 1
  z.standard <- as.matrix(expand.grid(a,b,c,d,e,f)) # orig
  M <- nrow(z.standard)
  #z <- as.matrix(expand.grid(a,b))
  
  
  #this z_design matrix is the second stage matrix 
  z.design.additive <- cbind(1,z.standard)
  M <- nrow(z.design.additive)
  additive.number <- 2
  additive.second.cat <- ncol(z.design.additive)
  total.covar.number <- 1+ additive.number
  z.all <- matrix(0,nrow=(M*total.covar.number),ncol = (M+                                                        additive.second.cat*additive.number))
  for(i in c("intercept",
             "additive"
  )){
    ##we always keep intercept as saturated model and to simply, we always use diagnonal matrix for intercept
    if(i=="intercept"){
      ###row start and column start point for this category
      row.start <- 0
      column.start <- 0
      for(j in 1:M){
        z.all[row.start+1+(j-1)*total.covar.number,(column.start+j)] = 1
      }
    }else if(i=="additive"){
      column.start <- M
      row.start <- 1
      if(additive.number!=0){
        for(j in 1:M){
          for(k in 1:additive.number){
            z.all[row.start+k+(j-1)*total.covar.number,
                  (column.start+(k-1)*additive.second.cat+1):
                    (column.start+k*additive.second.cat)] <- as.matrix(z.design.additive[j,])
          }
        }
      }
    }
  }
  
  
  K <- ncol(z.design.additive)
  
  # z <- kronecker(diag(2),z)
  #this theta is the true value
  theta <- c(theta_intercept,theta_test,theta_covar)
  
  #this is the true beta
  beta <- z.all%*%theta
  beta <- matrix(beta,nrow=p_col+2)
  #alpha <- c(0,rep(1,length(beta)-1))
  
  G <-  rbinom(n,2,0.25)
  x_covar <- rnorm(n)
  x_all <- cbind(1,G,x_covar) ##adding the intercept into the model
  
  predictor <- x_all%*%(beta)
  predictor <- cbind(0,predictor)
  #predictor <- sweep(predictor,2,alpha,"+")
  p <- exp(predictor)
  sp <- rowSums(p)
  #this standarize the probability sum into 1,since logit model:pr(D=1|predictor)=exp(predictor)/(1+exp(predictor))
  p <- sweep(p,1,sp,"/")
  y <- t(apply(p,1,function(x){rmultinom(1,1,x)}))
  
  y <- y[,-1] # this is the y matrix for the model
  y.case.control <- rowSums(y)
  y.tumor <- y%*%z.standard
  
  idx.control <- which(y.case.control==0)
  y.tumor[idx.control,] <- rep(NA,6)
  idx.case <- which(y.case.control==1)
  rate <- c(0.17,0.25,0.42,0.27,0.05,0.05)
  for(i in 1:6){
    idx.mis <-  sample(idx.case,size=length(idx.case)*rate[i])
    y.tumor[idx.mis,i] <- 888
  }
  y.pheno.mis <- cbind(y.case.control,y.tumor)
  return(list(y.pheno.mis=y.pheno.mis,
              G=G,
              x_covar=x_covar))  
}






Generatesubtypes<- function(ER,PR,HER2,Grade,T5,T6){
  n <- length(ER)
  subtypes <- rep("control",n)
  temp = 1
  for(i in 0:1){
    for(j in 0:1){
      for(k in 0:1){
        for(l in 1:3){
          for(q in 0:1){
            for(w in 0:1)
              idx <- which(ER==i&PR==j&HER2==k&Grade==l&T5==q&T6==w)
            if(length(idx)!=0){
              subtypes[idx] <- temp
              temp = temp+1  
            }
            
          }
          
          
        }
      }
    }
  }
  
  subtypes <- factor(subtypes,levels=c("control",
                                       c(1:temp)))
  sum <- table(subtypes)
  idx.cat <- which(sum<=10)
  idx.remove <- which((subtypes%in%(unique(idx.cat)-1))==T)
  if(length(idx.remove)==0){
    return(list(subtypes,idx.remove))
  }else{
    subtypes <- subtypes[-idx.remove]
    return(list(subtypes,idx.remove))  
  }
  
}







args = commandArgs(trailingOnly = T)
i1 = as.numeric(args[[1]])
i2 = as.numeric(args[[2]])
print(i1)
#setwd("/dcl01/chatterj/data/hzhang1/breast_cancer_data_analysis/")
setwd('/data/zhangh24/breast_cancer_data_analysis/')
library(bc2, lib.loc ="/home/zhangh24/R/x86_64-pc-linux-gnu-library/3.6/")

theta_intercept <- c(-6.51, -3.64, -3.71, -3.93, -4.74, -3.43, -4.45, -2.40, -3.60, -5.85,-1.20,-3.50, -4.51, -2.39, -4.46, -3.53, -5.95,-4.00, -3.62,-2.14,-5.14, -2.65, -3.88,-2.91,
                     -6.51, -3.64, -3.71, -3.93, -4.74, -3.43, -4.45, -2.40, -3.60, -5.85,-1.20,-3.50, -4.51, -2.39, -4.46, -3.53, -5.95,-4.00, -3.62,-2.14,-5.14, -2.65, -3.88,-2.91,
                     -6.51, -3.64, -3.71, -3.93, -4.74, -3.43, -4.45, -2.40, -3.60, -5.85,-1.20,-3.50, -4.51, -2.39, -4.46, -3.53, -5.95,-4.00, -3.62,-2.14,-5.14, -2.65, -3.88,-2.91,
                     -6.51, -3.64, -3.71, -3.93, -4.74, -3.43, -4.45, -2.40, -3.60, -5.85,-1.20,-3.50, -4.51, -2.39, -4.46, -3.53, -5.95,-4.00, -3.62,-2.14,-5.14, -2.65, -3.88,-2.91)-1.45

theta_covar <- c(0.05,0,0,0,0,0,0)
theta_test <- c(0.08, 0.08, 0.05, 0.05, 0.05,0.05,0.05)
s_times <- 10
odds1 <- matrix(0,s_times,7)
sigma1 <- matrix(0,s_times,7)
sizes <- c(25000,50000,100000)
n <- sizes[i2]
for(i in 1:s_times){
  print(i)
  
  simulate_data <- SimulateData(theta_intercept,theta_test,theta_covar,n)
  y.pheno.mis <- simulate_data[[1]]
  G <- simulate_data[[2]]
  x_covar <- simulate_data[[3]]
  model1 <- TwoStageModel(y.pheno.mis,additive=cbind(G,x_covar),missingTumorIndicator = 888)
  z.standard <- model1[[12]]
  M <- nrow(z.standard)
  K <- ncol(z.standard)+1
  odds <- model1[[1]][M+(1:(K))]
  sigma <-  (model1[[2]][M+(1:K),M+(1:K)])
  odds1[i,] <- odds
  sigma1[i,] <-  sigma
  
}

result <- list(odds1,sigma1)
save(result,file=paste0("./simulation/EM_algorithm_evaluation/result/simu_result_high_",i1,"_",i2,".Rdata"))
