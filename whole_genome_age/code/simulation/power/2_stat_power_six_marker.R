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
  e <- c(0,1)
  f <- c(0,1)
  ##number of the other covarites
  p_col <- 1
  z.standard <- as.matrix(expand.grid(a,b,c,d,e,f)) # orig
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




#try
# idx.try <- which(y.pheno.mis[,5]==1|y.pheno.mis[,5]==2|
#                    y.pheno.mis[,5]==3)
# 
# y.pheno.mis[idx.try,5] <- y.pheno.mis[idx.try,5]-2



#this power compare function contains the results for existing methods: standard logistic regression, and MTOP
PowerCompare <- function(y.pheno.mis,G,x_covar){
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
  score.test.support.fixed <- ScoreTestSupportMixedModel(
    y.pheno.mis,
    baselineonly = NULL,
    additive = as.matrix(x_covar),
    pairwise.interaction = NULL,
    saturated = NULL,
    missingTumorIndicator = 888
  )
  score.test.fixed<- ScoreTestMixedModel(y=y.pheno.mis,
                                         x=as.matrix(as.numeric(G)),
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
  
  ###############standard logistic regression#############
  model.standard <- glm(y.pheno.mis[,1]~G+x_covar,family = binomial(link='logit')) 
  p_standard <- summary(model.standard)$coefficients[2,4]
  #########################################################
  
  
  ##############polytomous logistic regression###########
  idx.mis <- GenerateMissingPosition(y.pheno.mis,missingTumorIndicator=888)
  y.pheno.com <- y.pheno.mis[-idx.mis,,drop=F]
  x.covar.com <- x_covar[-idx.mis,drop=F]
  G.com <- G[-idx.mis,drop=F]
  
  
  temp <-  Generatesubtypes(y.pheno.com[,2],y.pheno.com[,3],y.pheno.com[,4],y.pheno.com[,5],y.pheno.com[,6],y.pheno.com[,7])
  if(length(temp[[1]])<=2){
    p_poly = 1
  }else{
    subtypes <- temp[[1]]
    idx.remove <- temp[[2]]
    if(length(idx.remove)!=0){
      x.covar.poly <- x.covar.com[-idx.remove]
      G.poly <- G.com[-idx.remove]
      
    }else{
      x.covar.poly <- x.covar.com
      G.poly <- G.com
      
    }
    poly.model <- multinom(subtypes~G.poly+x.covar.poly,maxit = 1000)
    
    if(poly.model$convergence==0){
      tryCatch({
        poly.model.coef <- coef(poly.model)
        M <- nrow(poly.model.coef)
        p.covariate <- ncol(poly.model.coef)
        snp.cov <- vcov(poly.model)[2+p.covariate*(0:(M-1)),2+p.covariate*(0:(M-1))]
        snp.coef <- poly.model.coef[,2]
        
        result_temp <- DisplaySecondStageTestResult(snp.coef,snp.cov)  
        p_poly <- result_temp[length(result_temp)-1]
      },
      error = function(e){
        p_poly<- 1
      }
      
      )
      
      
    }else{
      p_poly = 1
    }
    
  }
  
  
  # idx.mis <- GenerateMissingPosition(y.pheno.mis,missingTumorIndicator=888)
  # y.pheno.com <- y.pheno.mis[-idx.mis,,drop=F]
  # x.covar.com <- x_covar[-idx.mis,drop=F]
  # G.com <- G[-idx.mis,drop=F]
  # 
  # model2 <- TwoStageModel(y.pheno.com,additive=cbind(G.com,x.covar.com),missingTumorIndicator =NULL,delta0 = c(theta_intercept,theta_test,
  #                                                                                                              theta_covar))
  # z.standard <- model2[[12]]
  # M <- nrow(z.standard)
  # odds <- model2[[1]][M+(1:5)]
  # sigma <-  (model2[[2]][M+(1:5),M+(1:5)])
  # fixed.result <- DisplaySecondStageTestResult(odds,sigma)
  # p_global_complete <- fixed.result[11]
  result <- list(p_mglobal,p_standard,p_poly)  
  return(result)
}


#calculate the p-value for based on score statistics and information matrix
ScoreGlobalTestForAssoc <- function(score,infor){
  infor <- as.matrix(infor)
  df <- length(score)
  GTA.stat <- score%*%solve(infor)%*%t(score)
  p.value.GTA <- pchisq(as.numeric(GTA.stat),df=df,lower.tail = F)
  places = 3
  power.number <- floor(-log10(p.value.GTA))+places
  if(p.value.GTA==0){
    return(1E-300)
  }else{
    ###format the output with three digits in total
    p.value.GTA <- round(p.value.GTA*10^power.number)/(10^power.number)
    
  }
  
  return(p.value.GTA)
  
}

#to perform random-effect test using mixture chi-square test
ScoreMixedGlobalTestForHeter <- function(score.casecase,infor.casecase){
  GTH.stat <- as.numeric(score.casecase%*%t(score.casecase)) #score*t(score)
  lambda <- eigen(infor.casecase)$values #the eigen value of information matrix
  
  p.value.GTH <- Saddle(GTH.stat,lambda) #p-value for linear combination of chi-square
  
  if(p.value.GTH==2){
    acc = 1e-09
    lim = 2000000
    
    #davies method is a numerical algorithm
    #lim is the number of monte carol runs
    #acc is the accuracy
    #you can increase the lim number and decrease the acc level if the function speed is really fast
    result <- davies(GTH.stat,lambda,lim = lim,acc=acc)
    p.value.GTH <- result[[3]]
    
    if(result[[2]]!=0){
      #if result[[2]] is not 0
      #davies method didn't converge. 
      #I noticed this problem happens for 6_102485159_G_A and 1_87277974_G_A;
      #These two SNPs are oncoarray only SNPs with allele frequency around 0.01
      #just put the p-value as 1 for these SNPs
      print("chisq p value accuracy could't be reached")
      p.value.GTH = 1
    }
    
    if(p.value.GTH <0){
      #davies method is a numerical algorithm
      #when p-value is really small, the value can be negative
      #if it happens, just put the p-value as 1e-09 which is the accuracy for the function
      p.value.GTH <- acc
    }
    #places <- 3
    #power.number <- floor(-log10(p.value.GTH))+places
    #p.value.GTH <- round(p.value.GTH*10^power.number)/(10^power.number)
    
    
  }
  
  return(p.value.GTH)
  
}


SubtypesLogistic <- function(y.pheno.mis,G,x_covar){
  p_num = ncol(y.pheno.mis)
  #put a sufficient large vector to store p-value
  p_vec = rep(0, 100)
  temp = 1
  for(p in 2:p_num){
    if(p==5){
      #for grade, we use 1, 2, 3
      for(l in 1:3){
        idx <- which(y[,p] == l| y[,1]==0)  
        model = lm(y[idx,1]~G[idx])
        p_vec[temp] = summary(model)$coefficients[2,4]
        temp = temp + 1
      }
      
    }else{
      #for non-grade, the value is either 0 or 1
      for(l in 0:1){
        idx <- which(y[,p] == l| y[,1]==0)  
        model = lm(y[idx,1]~G[idx])
        p_vec[temp] = summary(model)$coefficients[2,4]
        temp = temp + 1
      }  
    }
    
    
  }
  p_vec = p_vec[1:(temp-1)]
  n_test = length(p_vec)
  #put bonferoni correction for multiple testing problem
  p_sub = ifelse((min(p_vec)*n_test)<1, min(p_vec)*n_test, 1)
  return(p_sub)
}



#i1 controls the simulation seed
args = commandArgs(trailingOnly = T)
i1 = as.numeric(args[[1]])
set.seed(i1)
print(i1)

theta_test3 = c(0.041063835,-0.061340705,-0.117066388,0.134361714,-0.007079230,
                -0.0334, -0.0087,0.085007289,-0.018221958, -0.003921596,
                0.0151, 0.0175,-0.056921875,0.023980953,0.02157052, 0.02969268,
                -0.050138876,-0.003414437, 0.011953769,0.035195921, 0.002280791,
                0.008654017)/8*6


#setwd("/dcl01/chatterj/data/hzhang1/breast_cancer_data_analysis/")
#change to your own data directory
# setwd('/data/zhangh24/breast_cancer_data_analysis/')
setwd('/data/shengf2/simu_mtop/power/')
#library bc2 or TOP package in your directory
library(bc2, lib.loc ="/home/shengf2/R/4.1/library/")
library(CompQuadForm)
library(ACAT, lib.loc = '/home/shengf2/R/4.1/library/')
#library(TOP,lib.loc ="/home/zhangh24/R/x86_64-pc-linux-gnu-library/3.6/")
library(Rcpp)
library(RcppArmadillo)
sourceCpp("Saddle.cpp")

theta_intercept <- c(-6.51, -3.64, -3.71, -3.93, -4.74, -3.43, -4.45, -2.40, -3.60, -5.85,-1.20,-3.50, -4.51, -2.39, -4.46, -3.53, -5.95,-4.00, -3.62,-2.14,-5.14, -2.65, -3.88,-2.91,
                     -6.51, -3.64, -3.71, -3.93, -4.74, -3.43, -4.45, -2.40, -3.60, -5.85,-1.20,-3.50, -4.51, -2.39, -4.46, -3.53, -5.95,-4.00, -3.62,-2.14,-5.14, -2.65, -3.88,-2.91,
                     -6.51, -3.64, -3.71, -3.93, -4.74, -3.43, -4.45, -2.40, -3.60, -5.85,-1.20,-3.50, -4.51, -2.39, -4.46, -3.53, -5.95,-4.00, -3.62,-2.14,-5.14, -2.65, -3.88,-2.91,
                     -6.51, -3.64, -3.71, -3.93, -4.74, -3.43, -4.45, -2.40, -3.60, -5.85,-1.20,-3.50, -4.51, -2.39, -4.46, -3.53, -5.95,-4.00, -3.62,-2.14,-5.14, -2.65, -3.88,-2.91)-1.45


theta_covar <- c(0.05,0,0,0,0,0,0)
#####three scenaiors
#####s==1: no hetergeneity settings; the SNP is associated with all cancer subtypes with the same effect
#####s==2: one marker driven the heterogeneity. One major tumor marker drives the heterogeneity
#####s==3: both additive effect + interaction effect exist. One major tumor marker drives the hetergeneity, the other markers also contribute some effects
sc <- 3


#number of simulation replicates in each job
s_times <- 5
sizes <- c(25000,50000,100000)
n.sizes <- length(sizes)
#p_global_result <- rep(0,3*n.sizes*s_times)
p_mglobal_result <- matrix(0,s_times,sc*n.sizes)
p_standard <- matrix(0,s_times,sc*n.sizes)
p_poly <- matrix(0,s_times,sc*n.sizes)
#p_global_complete <- rep(0,3*n.sizes*s_times)
p_subtypes <- matrix(0, s_times, sc*n.sizes)
#p_poly <- rep(0,9*s_times)
p_topo <- matrix(0,s_times,sc*n.sizes)

#sizes <- c(5000,25000,50000,100000)



for(i in 1:s_times){
  temp <- 1  
  for(s in 1:sc){
    if(s==1){
      #theta_test <- c(0.25,0,0,0,0)
      theta_test <- c(0.08,0,0,0,0,0,0)
      #theta_test <- c(0.25,0,0,0,0)
    }else if(s==2){
      #theta_test <- c(0,0.25,0,0,0)
      theta_test <- c(0,0.08,0,0,0,0,0)
    }else{
      #theta_test <- c(c(0,0.25),rnorm(3,0,0.02))
      #theta_test <- c(c(0,0.08),rnorm(3,0,0.02))
      # theta_test <- c(0,0.06,rnorm(20,0,0.02))
      # theta_test <- c(0,rnorm(6,0,0.03),rnorm(15,0,0.05))
      theta_test <- theta_test3
    }
    for(n in sizes){
      
      temp.simu <- SimulateDataPower(theta_intercept,theta_test,theta_covar,n,s)
      y.pheno.mis <- temp.simu[[1]]
      G <- temp.simu[[2]]
      x_covar <- temp.simu[[3]]
      #data preprocessing for generate the phenotype data
      y <- y.pheno.mis
      # missing.data.vec <- GenerateMissingPosition(y,missingTumorIndicator=888)
      # y.pheno.complete <- y[-missing.data.vec,]
      # freq.subtypes <- GenerateFreqTable(y.pheno.complete)
      # idx.del <- which(freq.subtypes[,ncol(freq.subtypes)]<=10)
      # if(length(idx.del)!=0){
      #   theta_intercept_input = theta_intercept[-idx.del]
      # }else{
      #   theta_intercept_input = theta_intercept
      # }
      # print("simulation")
      
      model.result <- PowerCompare(y.pheno.mis,G,x_covar)
      #p_global_result[temp] <- as.numeric(model.result[[1]])
      
      #the result matrix contains 9 columns (2 scenarios* 3 sample*sizes)
      #the column 1-3 represents scenarior one (s==1), column 1-3 represents sample size 25000, 50000, 100000
      #the column 4-6 represents scenarior two (s==2), column 4-6 represents sample size 25000, 50000, 100000
      #the column 7-9 represents scenarior three (s==3), column 7-9 represents sample size 25000, 50000, 100000
      #each row of the results matrix contains a simulation replicate
      p_mglobal_result[i,temp] <- as.numeric(model.result[[1]])
      p_standard[i,temp] <- as.numeric(model.result[[2]])
      # p_poly[i,temp] <- as.numeric(model.result[[3]])
      tryCatch({
        p_poly[i,temp] <- as.numeric(model.result[[3]])
      },
      error = function(e){
        p_poly[i,temp] <- 1
      }
      )
      # p_poly[temp] <- as.numeric(model.result[[5]])
      
      p_subtypes[i,temp]  = SubtypesLogistic(y.pheno.mis, G, x_covar)
      
      ##################You need to write the function for TOPO in this section using the data############
      ##################the function input contains y.pheno.mis, G and x_covar############################
      ##################the output is the p-value from ACAT result########################################
      score.test.support.fixed <- ScoreTestSupport(
        y.pheno.mis,
        baselineonly = NULL,
        additive = as.matrix(x_covar),
        pairwise.interaction = NULL,
        saturated = NULL,
        missingTumorIndicator = 888
      )
      
      #fit the ScoreTest,  additive
      #change second.stage.structure to second.stage.structure = pairwise.interaction for interaction model
      #change second.stage.structure to second.stage.structure = saturated for saturated
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
      p_topo[i,temp] = ACAT(Pvals=p.values)
      
      
      temp = temp+1
      
    }

    # x_covar <- rnorm(n)
    
  }
  
  
  

}




result <- list(p_mglobal_result,p_standard,p_poly,p_subtypes,p_topo)

#change the output directory to your own directory
save(result,file=paste0("./result_six/simu_result_0.25_",i1,".Rdata"))




