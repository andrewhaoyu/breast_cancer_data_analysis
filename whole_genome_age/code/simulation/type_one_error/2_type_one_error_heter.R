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
  #two covariates are generated under additive structure
  additive.number <- 2
  #four tumor characteristics mimicking ER, PR, HER2 and grade
  additive.second.cat <- ncol(z.design.additive)
  #total covar number need to add one for intercept
  total.covar.number <- 1+ additive.number
  
  #prepare the first stage parameters
  beta <- c(beta_intercept,beta_covar)
  beta <- matrix(beta,nrow=p_col+1,byrow=T)
  
  #alpha <- c(0,rep(1,length(beta)-1))
  
  #n <- 100000
  #G <-  rbinom(n,2,0.25)
  #x_covar <- rnorm(n)
  x_all <- cbind(1,x_covar) ##adding the intercept into the model
  #calcualte X*beta
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
  #set up the missing rate for each of the four tumor characteristics
  rate <- c(0.17,0.25,0.42,0.27)
  for(i in 1:4){
    idx.mis <-  sample(idx.case,size=length(idx.case)*rate[i])
    y.tumor[idx.mis,i] <- 888
  }
  y.pheno.mis <- cbind(y.case.control,y.tumor)
  return(y.pheno.mis)  
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
  
  
  GTH.stat <- as.numeric(score.casecase%*%t(score.casecase))
  lamda <- eigen(infor.casecase)$values
  
  acc = 1e-09
  lim = 2000000
  
  #davies method is a numerical algorithm
  #lim is the number of monte carol runs
  #acc is the accuracy
  #you can increase the lim number and decrease the acc level if the function speed is really fast
  ### some improvement rather than davies
  result <- davies(GTH.stat,lamda,lim = lim,acc=acc)
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
  
  
  return(p.value.GTH)
  
  
}

### 2400 args
args = commandArgs(trailingOnly = T)
i1 = as.numeric(args[[1]])

set.seed(i1)

#setwd("/data/shengf2/simu_mtop/type1error/")
#library(bc2, lib.loc ="/home/shengf2/R/4.1/library/")
library(CompQuadForm)
#library(ACAT, lib.loc = '/home/shengf2/R/4.1/library/')
library(bc2)
library(ACAT)
beta_intercept <- c(-6.51, -3.64, -3.71, -3.93, -4.74, -3.43, -4.45, -2.40, -3.60, -5.85,-1.20,-3.50, -4.51, -2.39, -4.46, -3.53, -5.95,-4.00, -3.62,-2.14,-5.14, -2.65, -3.88,-2.91)
#theta_test <- -c(0.35, 0.15, 0.25, 0.05, 0.20)
beta_covar <- rep(0.05,24)

#SimulateData <- function(beta_intercept,beta_covar,x_covar,n)
#sizes <- 100000
#s_times <- 2
# library(foreach)
# library(doParallel)
# no.cores <- 2
# registerDoParallel(no.cores)


#result.list <- foreach(job.i = 1:2)%dopar%{
# s_times <- 2000
#simulated times as 100,000 
s_times = 10000
#simulate data under three different sample sizes
size <- c(5000,50000,100000)


p_global_result_ftop <- matrix(0, s_times, length(size))
p_global_result_mtop1 <- matrix(0, s_times, length(size))
p_global_result_mtop2 <- matrix(0, s_times, length(size))
p_acat <- matrix(0, s_times, length(size))


for(k in 1:length(size)){
  n = size[k]
  ## should I set a time seed for all three TOP models
  x_covar <- rnorm(n)
  y.pheno.mis <- SimulateData(beta_intercept,beta_covar,x_covar,n)
  y <- y.pheno.mis
  print("simulation")
  
  #get the three different z design matrix
  z.design.list = GenerateZDesignCombination(y.pheno.mis)
  z.additive = z.design.list[[1]]
  z.interaction = z.design.list[[2]]
  z.saturated = z.design.list[[3]]
  M = nrow(z.additive)
  #number of second stage parameters
  #additive model second. num is 5; 11--pairwise; 16--saturated

  #keep the parameter for baseline effect estimated from the model
  z.design.support <- matrix(rep(1,M),M,ncol=1)
  G <- rbinom(n,2,0.25)
  score.test.heter.support <- ScoreTestSupportMixedModelSelfDesign(
    y.pheno.mis,
    x.self.design  = G,
    z.design = z.design.support,
    additive =  as.matrix(x_covar),
    pairwise.interaction = NULL,
    saturated = NULL,
    missingTumorIndicator = 888
  )
  
  #######################################################
  print("finished")
  for(i in 1:s_times){
    print(i)
    #generate a genotype data
    G <- rbinom(n,2,0.25)
    
    ###########rewrite this section using the new model pipeline
    
    #fit the ScoreTest,  additive
    #change second.stage.structure to second.stage.structure = pairwise.interaction for interaction model
    #change second.stage.structure to second.stage.structure = saturated for saturated
    z.additive.heter = z.additive[,-1,drop = F]
    z.interaction.heter = z.interaction[,-1,drop = F]
    z.saturated.heter = z.saturated[,-1, drop = F]
    model.result1 <- ScoreTestMixedModel(y=y.pheno.mis,
                                                 x=as.matrix(G),
                                                 z.design=z.additive.heter,
                                                 score.test.support=score.test.heter.support,
                                                 missingTumorIndicator=888)
    #the first element is score
    score_result1 <- model.result1[[1]]
    #the second element is the efficient information matrix
    infor_result1 <- model.result1[[2]]
    
    model.result2 <- ScoreTestMixedModel(y=y.pheno.mis,
                                         x=as.matrix(G),
                                         z.design=z.interaction.heter,
                                         score.test.support=score.test.heter.support,
                                         missingTumorIndicator=888)
    #the first element is score
    score_result2 <- model.result2[[1]]
    #the second element is the efficient information matrix
    infor_result2 <- model.result2[[2]]
    
    model.result3 <- ScoreTestMixedModel(y=y.pheno.mis,
                                         x=as.matrix(G),
                                         z.design=z.saturated.heter,
                                         score.test.support=score.test.heter.support,
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
    
    p_global_result_ftop[i,k] <- as.numeric(pvalue1)
    p_global_result_mtop1[i,k] <- as.numeric(pvalue2)
    p_global_result_mtop2[i,k] <- as.numeric(pvalue3)
    
    
    ##ACAT combination to generate a global pvalue
    p.values = c(pvalue1,pvalue2,pvalue3)
    p_acat[i,k] = ACAT(Pvals=p.values)
  }
  
  
}

#save the results to your folder
save(p_global_result_ftop,file=paste0("./simulation/result1/ftop_",i1,".Rdata"))
save(p_global_result_mtop1,file=paste0("./simulation/result2/mtop1_",i1,".Rdata"))
save(p_global_result_mtop2,file=paste0("./simulation/result3/mtop2_",i1,".Rdata"))
save(p_acat,file=paste0("./simulation/result/acat",i1,".Rdata"))
