###########simulate data with four tumor characteristics containing three binary tumor characteristics and one oridinal tumor characteristics.
###########one genotype with MAF 0.25 is simulated
###########one covariate with rnorm(n) simulated
###########simulate results for MTOP, standard log istic regresssion, and TOPO

library(nnet)
#SmulateDataPower function generate phenotype and genotype for power anaysis
#theta_intercept sets the second-stage parameters for the intercepts
#theta_test sets the second-stage parameters for the SNP
#theta_covar sets the second-stage parameters for the other onecovariate, e.g. PC1
#n sets the sample size for the data (sizes = 25000, 50,000 and 100,000)
SimulateDataPower <- function(theta_intercept,theta_test,theta_covar,n,s){
  a <- c(0,1)
  b <- c(0,1)
  c <- c(0,1)
  d <- c(1,2,3)
  ##number of the other covarites
  p_col <- 1
  z.standard <- as.matrix(expand.grid(a,b,c,d)) # orig
  M <- nrow(z.standard)
  z.design.additive <- cbind(1,z.standard)
  z.design.pairwise.interaction <- z.design.additive
  all.combn <- combn(ncol(z.standard),2)
  z.temp <- matrix(0,nrow(z.standard),ncol(all.combn))
  for(i in 1:ncol(all.combn)){
    z.temp[,i] <- z.standard[,all.combn[1,i]]*z.standard[,all.combn[2,i]] 
  }
  z.design.pairwise.interaction <- cbind(z.design.pairwise.interaction,z.temp)
  #z <- as.matrix(expand.grid(a,b))
  
  if(s==1|s==2){
    #this z_design matrix is the second stage matrix 
    
    M <- nrow(z.design.additive)
    additive.number <- 2
    pairwise.interaction.number <- 0
    additive.second.cat <- ncol(z.design.additive)
    total.covar.number <- 1+ additive.number 
    z.all <- matrix(0,nrow=(M*total.covar.number),
                    ncol = (M+                                                        
                              additive.second.cat*additive.number))
  }else if(s==3){
    additive.number <- 1
    pairwise.interaction.number <- 1
    additive.second.cat <- ncol(z.design.additive)
    pairwise.interaction.second.cat <- ncol(z.design.pairwise.interaction)
    total.covar.number <- 1+ additive.number+
      pairwise.interaction.number
    z.all <- matrix(0,nrow=(M*total.covar.number),ncol = (M+
                                    additive.second.cat*additive.number)+
                                    pairwise.interaction.second.cat*pairwise.interaction.number)
  }

  for(i in c("intercept",
             "pairwise.interaction",
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
    }else if(i == "pairwise.interaction"){
      column.start <- M+additive.number*additive.second.cat
      row.start <- 1+additive.number
      if(pairwise.interaction.number!=0){
        for(j in 1:M){
          for(k in 1:pairwise.interaction.number){
            z.all[row.start+k+(j-1)*total.covar.number,
                  (column.start+(k-1)*pairwise.interaction.second.cat+1):
                    (column.start+k*pairwise.interaction.second.cat)] <- as.matrix(z.design.pairwise.interaction[j,])
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
  y.tumor[idx.control,] <- rep(NA,4)
  idx.case <- which(y.case.control==1)
  rate <- c(0.17,0.25,0.42,0.27)
  for(i in 1:4){
    idx.mis <-  sample(idx.case,size=length(idx.case)*rate[i])
    y.tumor[idx.mis,i] <- 888
  }
  y.pheno.mis <- cbind(y.case.control,y.tumor)
  return(list(y.pheno.mis=y.pheno.mis,
              G=G,
              x_covar=x_covar))  
}






Generatesubtypes<- function(ER,PR,HER2,Grade){
  n <- length(ER)
  subtypes <- rep("control",n)
  temp = 1
  for(i in 0:1){
    for(j in 0:1){
      for(k in 0:1){
        for(l in 1:3){
          
          idx <- which(ER==i&PR==j&HER2==k&Grade==l)
          if(length(idx)!=0){
            subtypes[idx] <- temp
            temp = temp+1  
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




#try
# idx.try <- which(y.pheno.mis[,5]==1|y.pheno.mis[,5]==2|
#                    y.pheno.mis[,5]==3)
# 
# y.pheno.mis[idx.try,5] <- y.pheno.mis[idx.try,5]-2



#this power compare function contains the results for existing methods: MTOP, standard logistic regression, and polytomous model
MTOP <- function(y.pheno.mis,G,x_covar, score.test.support.fixed){
  # model1 <- TwoStageModel(y.pheno.mis,
  #                         additive=cbind(G,x_covar),
  #                         missingTumorIndicator = 888,
  #                         delta0 = c(theta_intercept,theta_test,theta_covar))
  # z.standard <- model1[[12]]
  # M <- nrow(z.standard)
  # odds <- model1[[1]][M+(1:5)]
  # sigma <-  (model1[[2]][M+(1:5),M+(1:5)])
  # fixed.result <- DisplaySecondStageTestResult(odds,sigma)
  # p_global <- fixed.result[11]
  # p_heter <- fixed.result[12]
  # p_indi <- fixed.result[2]
  ##########MTOP global test for association################## 
  z.standard = GenerateZstandard(y.pheno.mis)
  M = nrow(z.standard)
  z.design.fixed <- cbind(rep(1,M),z.standard[,1])
  z.design.random <-as.matrix(z.standard[,2:ncol(z.standard)])
  
  score.test.fixed<-  ScoreTest(y=y.pheno.mis,
                                x=G,
                                second.stage.structure="additive",
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
    missingTumorIndicator = 888
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
  
  ######################################################
  
 
    
  
    
    result <- list(p_mglobal)  
    return(result)
  }
  
#TOPO function
TOPO <- function(y.pheno.mis,G, score.test.support.fixed){
  model.result1 <- ScoreTest(y=y.pheno.mis,
                             x=G,
                             second.stage.structure="additive",
                             score.test.support=score.test.support.fixed,
                             missingTumorIndicator=888)
  #the first element is score
  score_result1 <- model.result1[[1]]
  #the second element is the efficient information matrix
  infor_result1 <- model.result1[[2]]
  
  model.result2 <- ScoreTest(y=y.pheno.mis,
                             x=G,
                             second.stage.structure="pairwise.interaction",
                             score.test.support=score.test.support.fixed,
                             missingTumorIndicator=888)
  #the first element is score
  score_result2 <- model.result2[[1]]
  #the second element is the efficient information matrix
  infor_result2 <- model.result2[[2]]
  
  model.result3 <- ScoreTest(y=y.pheno.mis,
                             x=G,
                             second.stage.structure="saturated",
                             score.test.support=score.test.support.fixed,
                             missingTumorIndicator=888)
  #the first element is score
  score_result3 <- model.result3[[1]]
  #the second element is the efficient information matrix
  infor_result3 <- model.result3[[2]]
  
  ## p.value for fixed-effect test using chi-square test
  pvalue1 = ScoreGlobalTestForAssoc(score_result1, infor_result1)
  ## p.value for random-effect test using mixture chi-square test
  pvalue2 = ScoreMixedGlobalTestForHeter(score_result2, infor_result2)
  pvalue3 = ScoreMixedGlobalTestForHeter(score_result3, infor_result3)
  
  pvalue1 <- as.numeric(pvalue1)
  pvalue2 <- as.numeric(pvalue2)
  pvalue3 <- as.numeric(pvalue3)
  
  ##ACAT combination to generate a global pvalue
  p.values = c(pvalue1,pvalue2,pvalue3)
  return(p.values)
  
}










  #i1 controls the simulation seed
  args = commandArgs(trailingOnly = T)
  i1 = as.numeric(args[[1]])
  set.seed(i1)
  print(i1)
  #setwd("/dcl01/chatterj/data/hzhang1/breast_cancer_data_analysis/")
  #change to your own data directory
  setwd('/data/zhangh24/breast_cancer_data_analysis/')
  #library(bc2)
  #library bc2 or TOP package in your directory
  library(bc2, lib = )
  #library(TOP,lib.loc ="/home/zhangh24/R/x86_64-pc-linux-gnu-library/3.6/")
  theta_intercept <- c(-6.51, -3.64, -3.71, -3.93, -4.74, -3.43, -4.45, -2.40, -3.60, -5.85,-1.20,-3.50, -4.51, -2.39, -4.46, -3.53, -5.95,-4.00, -3.62,-2.14,-5.14, -2.65, -3.88,-2.91)
  
  theta_covar <- c(0.05,0,0,0,0)
  
  #number of simulation replicates in each job
  s_times <- 100
  sizes <- c(25000,50000,100000)
  n.sizes <- length(sizes)
  #p_global_result <- rep(0,3*n.sizes*s_times)
  time_mtop <- matrix(0,s_times, n.sizes)
  time_topo <- matrix(0,s_times, n.sizes)
  #put p_topo here
  #p_topo
  
  

    
  #for computation time, we just need one setting
    theta_test <- c(0,0,0,0,0)
     
    temp <- 1  
  #run MTOP  
    for(s_size in 1:length(sizes)){
      n = sizes[s_size]
      G_all = matrix(rbinom(n*s_times, size = 1, prob =0.25), n, s_times)
      temp.simu <- SimulateDataPower(theta_intercept,theta_test,theta_covar,n,1)
      y.pheno.mis <- temp.simu[[1]]
      x_covar <- temp.simu[[3]]
      y <- y.pheno.mis
      #fit the null model
      score.test.support.fixed <- ScoreTestSupport(
        y.pheno.mis,
        baselineonly = NULL,
        additive = as.matrix(x_covar),
        pairwise.interaction = NULL,
        saturated = NULL,
        missingTumorIndicator = 888
      )
      
      time_start = proc.time()
      for(i in 1:s_times){  
        G <- G_all[,i]
        print("simulation")
        
        model.result <- MTOP(y.pheno.mis,G,x_covar, score.test.support.fixed)
    
      }
     
      mtop_run_time = (proc.time() - time_start)[3]
      time_mtop[s] = mtop_run_time
    
    
  
    #run TOPO
      time_start = proc.time()
      for(i in 1:s_times){  
        G <- G_all[,i]
        print("simulation")
        
        model.result <- TOPO(y.pheno.mis,G,score.test.support.fixed)
        
      }
      
      topo_run_time = (proc.time() - time_start)[3]
      time_topo[s] = topo_run_time
    }
 

result <- list(time_mtop,time_topo)

#change the output directory to your own directory
#save(result.list,file=paste0("./simulation/power/result/simu_result_0.25_",i1,".Rdata"))

