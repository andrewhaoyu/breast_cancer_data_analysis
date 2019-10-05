###########simulate data with four tumor characteristics containing five binary tumor characteristics and one oridinal tumor characteristics.
###########one genotype with MAF 0.25 is simulated
###########one covariate with rnorm(n) simulated

library(nnet)
SimulateDataPower <- function(theta_intercept,theta_test,theta_covar,n){
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




PowerCompare <- function(y.pheno.mis,G,x_covar,theta_intercept,theta_test,theta_covar){
    
    model1 <- TwoStageModel(y.pheno.mis,
                            additive=cbind(G,x_covar),
                            missingTumorIndicator = 888,
                            delta0 = c(theta_intercept,theta_test,theta_covar))
    
      z.standard <- model1[[12]]
      M <- nrow(z.standard)
      K <- ncol(z.standard)+1
      odds <- model1[[1]][M+(1:(K))]
      sigma <-  (model1[[2]][M+(1:K),M+(1:K)])
      fixed.result <- DisplaySecondStageTestResult(odds,sigma)
      p_global <- fixed.result[length(fixed.result)-1]
      
    
    # p_heter <- fixed.result[12]
    # p_indi <- fixed.result[2]
    ##########MTOP global test for association 
    z.design.fixed <- cbind(rep(1,M),z.standard[,1])
    z.design.random <-z.standard[,2:ncol(z.standard)]
    score.test.support.fixed <- ScoreTestSupportMixedModel(
      y.pheno.mis,
      baselineonly = NULL,
      additive = as.matrix(x_covar),
      pairwise.interaction = NULL,
      saturated = NULL,
      missingTumorIndicator = 888,
      c(theta_intercept,
        theta_covar)
    )
    score.test.fixed<- ScoreTestMixedModel(y=y.pheno.mis,
                                           x=as.matrix(G),
                                           z.design=z.design.fixed,
                                           score.test.support=score.test.support.fixed,
                                           missingTumorIndicator=888
    )
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
      delta0 = c(theta_intercept,theta_test[1:ncol(z.design.fixed)],
                 theta_covar)
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
    model.standard <- glm(y.pheno.mis[,1]~G+x_covar,family = binomial(link='logit')) 
    p_standard <- summary(model.standard)$coefficients[2,4]
    
    idx.mis <- GenerateMissingPosition(y.pheno.mis,missingTumorIndicator=888)
    y.pheno.com <- y.pheno.mis[-idx.mis,,drop=F]
    x.covar.com <- x_covar[-idx.mis,drop=F]
    G.com <- G[-idx.mis,drop=F]
    
    model2 <- TwoStageModel(y.pheno.com,additive=cbind(G.com,x.covar.com),missingTumorIndicator =NULL,delta0 = c(theta_intercept,theta_test,theta_covar))
    z.standard <- model2[[12]]
    M <- nrow(z.standard)
    odds <- model2[[1]][M+(1:K)]
    sigma <-  (model2[[2]][M+(1:K),M+(1:K)])
    fixed.result <- DisplaySecondStageTestResult(odds,sigma)
    p_global_complete <- fixed.result[length(fixed.result)-1]
 

 
  result <- list(p_global,p_mglobal,p_standard,p_global_complete)  
  return(result)
}

  # if(length(temp[[1]])<=2){
  #   p_poly = 1
  # }else{
  #   subtypes <- temp[[1]]
  #   idx.remove <- temp[[2]]
  #   if(length(idx.remove)!=0){
  #     x.covar.poly <- x.covar.com[-idx.remove]
  #     G.poly <- G.com[-idx.remove]
  #     
  #   }else{
  #     x.covar.poly <- x.covar.com
  #     G.poly <- G.com
  #     
  #   }
  #   poly.model <- multinom(subtypes~G.poly+x.covar.poly,maxit = 1000)
  #   
  #   if(poly.model$convergence==0){
  #     tryCatch({
  #       poly.model.coef <- coef(poly.model)
  #       M <- nrow(poly.model.coef)
  #       p.covariate <- ncol(poly.model.coef)
  #       snp.cov <- vcov(poly.model)[2+p.covariate*(0:(M-1)),2+p.covariate*(0:(M-1))]
  #       snp.coef <- poly.model.coef[,2]
  #       
  #       result_temp <- DisplaySecondStageTestResult(snp.coef,snp.cov)  
  #       p_poly <- result_temp[length(result_temp)-1]
  #     },
  #     error = function(e){
  #       p_poly<- 1
  #     }
  #     
  #     )
  #     
  #     
  #   }else{
  #     p_poly = 1
  #   }
  #   
  #   
  #   
  # }
  # 
  
 

# GenerateComplete <- function(y.pheno.mis,x_covar){
#   idx.mis <- which(y.pheno.mis[,2]==888|y.pheno.mis[,3]==888|
#                      y.pheno.mis[,4]==888|y.pheno.mis[,5]==888)
#   y.pheno.complete <- y.pheno.mis[-idx.mis,]
#   x.covar.complete <- x_covar[-idx.mis,]
#   return(list(y.pheno.complete,x.covar.complete))
# }



args = commandArgs(trailingOnly = T)
i1 = as.numeric(args[[1]])
print(i1)
#setwd("/dcl01/chatterj/data/hzhang1/breast_cancer_data_analysis/")
setwd('/spin1/users/zhangh24/breast_cancer_data_analysis/')
library(bc2, lib.loc ="/home/zhangh24/R/x86_64-pc-linux-gnu-library/3.6/")

theta_intercept <- c(-6.51, -3.64, -3.71, -3.93, -4.74, -3.43, -4.45, -2.40, -3.60, -5.85,-1.20,-3.50, -4.51, -2.39, -4.46, -3.53, -5.95,-4.00, -3.62,-2.14,-5.14, -2.65, -3.88,-2.91,
                     -6.51, -3.64, -3.71, -3.93, -4.74, -3.43, -4.45, -2.40, -3.60, -5.85,-1.20,-3.50, -4.51, -2.39, -4.46, -3.53, -5.95,-4.00, -3.62,-2.14,-5.14, -2.65, -3.88,-2.91,
                     -6.51, -3.64, -3.71, -3.93, -4.74, -3.43, -4.45, -2.40, -3.60, -5.85,-1.20,-3.50, -4.51, -2.39, -4.46, -3.53, -5.95,-4.00, -3.62,-2.14,-5.14, -2.65, -3.88,-2.91,
                     -6.51, -3.64, -3.71, -3.93, -4.74, -3.43, -4.45, -2.40, -3.60, -5.85,-1.20,-3.50, -4.51, -2.39, -4.46, -3.53, -5.95,-4.00, -3.62,-2.14,-5.14, -2.65, -3.88,-2.91)-1.45

theta_covar <- c(0.05,0,0,0,0,0,0)
#####three scenaiors
sc <- 3

library(foreach)
library(doParallel)
no.cores <- 2
registerDoParallel(no.cores)


result.list <- foreach(job.i = 1:2)%dopar%{
  set.seed(2*i1-job.i)
  s_times <- 1
  #sizes <- c(5000,25000,50000,100000)
 #sizes <- c()
  sizes <- c(5000)
  n.sizes <- length(sizes)
  p_global_result <- rep(0,n.sizes*s_times)
  p_mglobal_result <- rep(0,n.sizes*s_times)
  p_standard <- rep(0,n.sizes*s_times)
  p_global_complete <- rep(0,n.sizes*s_times)
  #p_poly <- rep(0,n.sizes*s_times)
  
  
  
  temp <- 1  
  for(s in 1:sc){
    if(s==1){
      #theta_test <- c(0.08,0,0,0,0,0,0)
      theta_test <- c(0.25,0,0,0,0,0,0)
      }else if(s==2){
      #theta_test <- c(0,0.08,0,0,0,0,0)
        theta_test <- c(0,0.25,0,0,0,0,0)
        }else{
      #theta_test <- c(c(0,0.08),rnorm(5,0,0.02))
      theta_test <- c(c(0,0.25),rnorm(5,0,0.02))
    }
    for(n in sizes){
      for(i in 1:s_times){
        temp.simu <- SimulateDataPower(theta_intercept,theta_test,theta_covar,n)
        y.pheno.mis <- temp.simu[[1]]
        table(y.pheno.mis[,1])
        G <- temp.simu[[2]]
        x_covar <- temp.simu[[3]]
        
        y <- y.pheno.mis
        missing.data.vec <- GenerateMissingPosition(y,missingTumorIndicator=888)
        y.pheno.complete <- y[-missing.data.vec,]
        freq.subtypes <- GenerateFreqTable(y.pheno.complete)
        idx.del <- which(freq.subtypes[,ncol(freq.subtypes)]<=10)
        if(length(idx.del)!=0){
          theta_intercept_input = theta_intercept[-idx.del]
        }else{
          theta_intercept_input = theta_intercept
        }
        print("simulation")
        
        model.result <- PowerCompare(y.pheno.mis,G,x_covar,
                                     theta_intercept_input,
                                     theta_test,
                                     theta_covar)
        p_global_result[temp] <- as.numeric(model.result[[1]])
        p_mglobal_result[temp] <- as.numeric(model.result[[2]])
        p_standard[temp] <- as.numeric(model.result[[3]])
        p_global_complete[temp] <- as.numeric(model.result[[4]])
       # p_poly[temp] <- as.numeric(model.result[[5]])
        temp = temp+1
        
      }
      # x_covar <- rnorm(n)
      
    }
    
  }
  
  result <- list(p_global_result,p_mglobal_result,p_standard,p_mglobal_result,p_global_complete)
  return(result)
}

stopImplicitCluster()
save(result.list,file=paste0("./simulation/power/result/simu_high_0.25_",i1,".Rdata"))

