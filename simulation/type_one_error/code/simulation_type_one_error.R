#This is an example with 1000 samples;
#There are 3 covaraites in total;p_col means the number of covariates;
#Three second stage categories;

rm(list=ls())
args <- commandArgs(trailingOnly = T)
i1 <- as.numeric(args[[1]])

set.seed(i1)

n.simulation <- 100
p.value.simulation <- matrix(0,n.simulation,6)
for(i.simu in 1:n.simulation){
  a <- c(0,1)
  b <- c(0,1)
  c <- c(0,1)
  p_col <- 3
  z <- as.matrix(expand.grid(a,b,c)) # orig
  #z <- as.matrix(expand.grid(a,b))
  
  
  #this z_design matrix is the second stage matrix 
  z_design <- cbind(1,z)
  M <- nrow(z_design)
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
  K <- ncol(z_design)
  
  # z <- kronecker(diag(2),z)
  theta_intercept <- runif(M,-3,-1)
  theta_test <- rep(0,K)
  theta_covar <- rep(0.2,K*(p_col-1))
  
  #this theta is the true value
  theta <- c(theta_intercept,theta_test,theta_covar)
  
  #this is the true beta
  beta <- z_all%*%theta
  beta <- matrix(beta,nrow=p_col+1)
  #alpha <- c(0,rep(1,length(beta)-1))
  
  n <- 100000
  x <-  matrix(rnorm(p_col*n),nrow = n)
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
  
  y <- y[,-1] # this is the y matrix for the model
  
  y.case.control <- rowSums(y)
  n.case <- sum(y.case.control)
  
  
  y.tumor <- matrix(0,n,ncol(z))
  for(i in 1:n){
    if(y.case.control[i]==1){
      temp.id <- which(y[i,]==1)
      y.tumor[i,] = z[temp.id,]
    }else{
      y.tumor[i,] <- rep(NA,ncol(z))
    }
  }
  
  y.pheno <- cbind(y.case.control,y.tumor)
  tumor <- K-1
  idx.case <- which(y.case.control==1)
  idx.mis <- matrix(rbinom(n.case*tumor,1,0.2),n.case,tumor)
  
  for(i in 1:ncol(idx.mis)){
    y.tumor[idx.case,i][idx.mis[,i]] <- 888
  }
  
  y.pheno.mis <- as.data.frame(cbind(y.case.control,y.tumor))
  colnames(y.pheno.mis) <- c("behavior","ER","PR","HER2")
  x <- as.data.frame(x)
  colnames(x) <- c("gene","pc1","pc2")
  
  library(bc2)
  result <- EMmvpoly(y.pheno.mis,baselineonly = NULL,additive = x,pairwise.interaction = NULL,saturated = NULL,missingTumorIndicator = 888)
  
  p.value <- c(as.numeric(result[[5]][1,2:3]),result[[4]][,6][1:4])
  p.value.simulation[i.simu,] <- p.value
}

save(p.value.simulation,file=paste0("/spin1/users/zhangh24/breast_cancer_data_analysis/simulation/type_one_error/result/pvalue",i1,".Rdata"))


# #result <- Mvpoly(delta0,y,x,z_design)
# 
# print("Finished mv_poly")
# 
# delta0 <- c(theta_intercept+runif(K),theta_covar+runif(length(theta_covar)))
# #this gives you the score test result
# #p_value <- score_test(delta0,y,x_test,x_covar,z_design)
# 
# result <- score_test(delta0,y,x_test,x_covar,z_design)
# 
# 
# 
# I_M <- diag(M)
# p_covar <- ncol(x_covar)
# 
# z_covar_temp <- NULL
# 
# for(i in 1:M){
#   z_covar_temp <- rbind(z_covar_temp,kronecker(diag(p_covar),t(z_design[i,])))
# }
# 
# z_covar <- matrix(0,nrow = M*(p_covar+1),ncol= M+p_covar*ncol(z_design))
# for(i in 1:M){
#   z_covar[1+(i-1)*(p_covar+1),i] <- 1
# }
# for(i in 1:M){
#   z_covar[(2+(i-1)*(p_covar+1)):(i*(p_covar+1)),(M+1):ncol(z_covar)] <- 
#     z_covar_temp[(1+(i-1)*p_covar):(i*p_covar),]
# }
# x_covar <- cbind(1,x_covar)
# xx_covar <- kronecker(I_M,x_covar)
# xxz_covar <- xx_covar%*%z_covar
# 
# 
# # sof <- "source.so"
# # dyn.load(sof)
# # 
# # 
# # 
# # nparm <- length(delta0)
# # NITER <- 100
# # tol   <- 1e-4
# # NCOV  <- ncol(x_covar)
# # ncat  <- NCOL(z_design)-1
# # DEBUG     <- 2
# # ret_rc    <- as.integer(1)
# # ret_delta <- as.numeric(rep(-9999, nparm))
# # ret_p     <- as.numeric(1)
# # 
# # temp <- .C("score_test", as.numeric(delta0), as.integer(nparm), as.numeric(as.vector(y)),
# #            as.numeric(x_covar), as.numeric(x_test), as.numeric(z_design), as.integer(n), as.integer(M), 
# #            as.integer(ncat), as.integer(NCOV), as.integer(NITER), as.numeric(tol), 
# #            as.integer(DEBUG), ret_rc=ret_rc, ret_delta=ret_delta, ret_p=ret_p)
# # print(paste("rc = ", temp$ret_rc, sep=""))
# # print(paste("pval = ", temp$ret_p, sep=""))
# 
# 
# 
