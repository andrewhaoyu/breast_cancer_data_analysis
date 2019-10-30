###########simulate data with four tumor characteristics containing three binary tumor characteristics and one oridinal tumor characteristics.
###########one genotype with MAF 0.25 is simulated
###########one covariate with rnorm(n) simulated



SimulateData <- function(beta_intercept,beta_covar,x_covar,n){
  a <- c(0,1)
  b <- c(0,1)
  c <- c(0,1)
  d <- c(1,2,3)
  ##number of the other covarites
  p_col <- 1
  z.standard <- as.matrix(expand.grid(a,b,c,d)) # orig
  M <- nrow(z.standard)
  #z <- as.matrix(expand.grid(a,b))
  
  
  #this z_design matrix is the second stage matrix 
  z.design.additive <- cbind(1,z.standard)
  M <- nrow(z.design.additive)
  additive.number <- 2
  additive.second.cat <- ncol(z.design.additive)
  total.covar.number <- 1+ additive.number
  
  
  beta <- c(beta_intercept,beta_covar)
  beta <- matrix(beta,nrow=p_col+1,byrow=T)
  
  #alpha <- c(0,rep(1,length(beta)-1))
  
  #n <- 100000
  #G <-  rbinom(n,2,0.25)
  #x_covar <- rnorm(n)
  x_all <- cbind(1,x_covar) ##adding the intercept into the model
  
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
  y.tumor[idx.control,] <- rep(NA,4)
  idx.case <- which(y.case.control==1)
  rate <- c(0.17,0.25,0.42,0.27)
  for(i in 1:4){
    idx.mis <-  sample(idx.case,size=length(idx.case)*rate[i])
    y.tumor[idx.mis,i] <- 888
  }
  y.pheno.mis <- cbind(y.case.control,y.tumor)
  return(y.pheno.mis)  
}


FixedMixedTwoStageModel <- function(y.pheno.mis,G,x_covar,score.test.support.fixed,beta_intercept,theta_test,theta_covar){
  model1 <- TwoStageModel(y.pheno.mis,
                          additive=cbind(G,x_covar),
                          missingTumorIndicator = 888,
                          delta0 = c(beta_intercept,theta_test,theta_covar))
  z.standard <- model1[[12]]
  M <- nrow(z.standard)
  odds <- model1[[1]][M+(1:5)]
  sigma <-  (model1[[2]][M+(1:5),M+(1:5)])
  fixed.result <- DisplaySecondStageTestResult(odds,sigma)
  p_global <- fixed.result[11]
  p_heter <- fixed.result[12]
  p_indi <- fixed.result[2]
  ##########MTOP global test for association 
  z.design.fixed <- cbind(rep(1,M),z.standard[,1])
  z.design.random <-z.standard[,2:ncol(z.standard)]
  score.test.fixed<- ScoreTestMixedModel(y=y.pheno.mis,
                                         x=as.matrix(G),
                                         z.design=z.design.fixed,
                                         score.test.support=score.test.support.fixed,
                                         missingTumorIndicator=888)
  score.fixed <- score.test.fixed[[1]]
  infor.fixed <- score.test.fixed[[2]]
  
  score.test.support.random <- ScoreTestSupportMixedModelSelfDesign(
    y.pheno.mis,
    x.self.design  = G,
    z.design = z.design.fixed,
    additive =  as.matrix(x_covar),
    pairwise.interaction = NULL,
    saturated = NULL,
    missingTumorIndicator = 888,
    delta0 = c(beta_intercept,theta_test[1:ncol(z.design.fixed)],theta_covar)
  )
  score.test.random<- ScoreTestMixedModel(y=y.pheno.mis,
                                          x=as.matrix(G),
                                          z.design=z.design.random,
                                          score.test.support=score.test.support.random,
                                          missingTumorIndicator=888)
  score.random <- score.test.random[[1]]
  infor.random <- score.test.random[[2]]
  p_mglobal <- DisplayMixedScoreTestResult(score.fixed,
                                           infor.fixed,
                                           score.random,
                                           infor.random)[1]
  ##########MTOP global test for heterogeneity
  z.design.support <- matrix(rep(1,M),M,ncol=1)
  z.design.fixed <- z.standard[,1,drop=F]
  score.test.heter.support <- ScoreTestSupportMixedModelSelfDesign(
    y.pheno.mis,
    x.self.design  = G,
    z.design = z.design.support,
    additive =  as.matrix(x_covar),
    pairwise.interaction = NULL,
    saturated = NULL,
    missingTumorIndicator = 888,
    delta0 = c(beta_intercept,theta_test[1:ncol(z.design.support)],theta_covar)
  )
  score.test.heter.fixed<- ScoreTestMixedModel(y=y.pheno.mis,
                                               x=as.matrix(G),
                                               z.design=z.design.fixed,
                                               score.test.support=score.test.heter.support,
                                               missingTumorIndicator=888)
  score.fixed.heter <- score.test.heter.fixed[[1]]
  infor.fixed.heter <- score.test.heter.fixed[[2]]
  p_mheter <- DisplayMixedScoreTestResult(score.fixed.heter,
                                          infor.fixed.heter,
                                          score.random,
                                          infor.random)[1]
result <- list(p_global,p_heter,p_indi,p_mglobal,p_mheter)  
return(result)
}


# GenerateComplete <- function(y.pheno.mis,x_covar){
#   idx.mis <- which(y.pheno.mis[,2]==888|y.pheno.mis[,3]==888|
#                      y.pheno.mis[,4]==888|y.pheno.mis[,5]==888)
#   y.pheno.complete <- y.pheno.mis[-idx.mis,]
#   x.covar.complete <- x_covar[-idx.mis,]
#   return(list(y.pheno.complete,x.covar.complete))
# }



# args = commandArgs(trailingOnly = T)
# i1 = as.numeric(args[[1]])

setwd('/data/zhangh24/breast_cancer_data_analysis/')
library(bc2)

beta_intercept <- c(-6.51, -3.64, -3.71, -3.93, -4.74, -3.43, -4.45, -2.40, -3.60, -5.85,-1.20,-3.50, -4.51, -2.39, -4.46, -3.53, -5.95,-4.00, -3.62,-2.14,-5.14, -2.65, -3.88,-2.91)
#theta_test <- -c(0.35, 0.15, 0.25, 0.05, 0.20)
beta_covar <- rep(0.05,24)
#SimulateData <- function(beta_intercept,beta_covar,x_covar,n)
#sizes <- 100000
#s_times <- 2
library(foreach)
library(doParallel)
no.cores <- 2
registerDoParallel(no.cores)

theta_test <- rep(0,5)
theta_covar <- c(0.05,0,0,0,0)
delta_ini = c(beta_intercept,theta_test,beta_covar)

result.list <- foreach(job.i = 1:2)%dopar%{
  set.seed(2*i1-job.i)
  s_times <- 2000
  p_global_result <- rep(0,3*s_times)
  p_heter_result <- rep(0,3*s_times)
  p_indi_result <- rep(0,3*s_times)
  p_mglobal_result <- rep(0,3*s_times)
  p_mheter_result <- rep(0,3*s_times)
  sizes <- c(5000,50000,100000)
  
    temp <- 1  
  for(n in sizes){
    x_covar <- rnorm(n)
    y.pheno.mis <- SimulateData(beta_intercept,beta_covar,x_covar,n)
   y <- y.pheno.mis
   missing.data.vec <- GenerateMissingPosition(y,missingTumorIndicator=888)
   y.pheno.complete <- y[-missing.data.vec,]
   freq.subtypes <- GenerateFreqTable(y.pheno.complete)
   idx.del <- which(freq.subtypes[,ncol(freq.subtypes)]<=10)
   if(length(idx.del)!=0){
     beta_intercept_input = beta_intercept[-idx.del]
   }else{
     beta_intercept_input = beta_intercept
   }
    print("simulation")
    score.test.support.fixed <- ScoreTestSupportMixedModel(
      y.pheno.mis,
      baselineonly = NULL,
      additive = as.matrix(x_covar),
      pairwise.interaction = NULL,
      saturated = NULL,
      missingTumorIndicator = 888,
      delta0 = c(beta_intercept_input,theta_covar)
    )
    print("finished")
    for(i in 1:s_times){
      print(i)
      G <- rbinom(n,2,0.25)
      
      ###########TOP model
      model.result <- FixedMixedTwoStageModel(y.pheno.mis,G,x_covar,score.test.support.fixed,beta_intercept_input,theta_test,theta_covar)
      p_global_result[temp] <- as.numeric(model.result[[1]])
      p_heter_result[temp] <- as.numeric(model.result[[2]])
      p_indi_result[temp] <- as.numeric(model.result[[3]])
      p_mglobal_result[temp] <- as.numeric(model.result[[4]])
      p_mheter_result[temp] <- as.numeric(model.result[[5]])
      
      temp = temp+1
      
    }
    
  }
  
  result <- list(p_global_result,p_heter_result,p_indi_result,p_mglobal_result,p_mheter_result)
  return(result)
}
stopImplicitCluster()

  save(result.list,file=paste0("./simulation/type_one_error/result/simu_result",i1,".Rdata"))
