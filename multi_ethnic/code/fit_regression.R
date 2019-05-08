#Goal: Generate phenotypes data for
#simulate phenotypes data for AFR, EUR, LAC
#MAF based on 1000KG
#sample size EUR n =120000
#sample size AFR n = 18000
#sample size LAC n = 18000
#heritability for EUR 0.5
#heritability for AFR 0.5
#heritability for LAC 0.5
#5000 causal SNPs for each population
#4000 shared causal SNPs
#Genetic correlation for the between EUR and AFR is 0.4
#Genetic correlation for the between EUR and LAC is 0.6
#Genetic correlation for the between LAC and AFR is 0.6
setwd('/spin1/users/zhangh24/breast_cancer_data_analysis')
arg <- commandArgs(trailingOnly=T)
i1 <- as.numeric(arg[[1]])

load(paste0("./multi_ethnic/result/pruned_geno/geno_",i1))
for(k in 1:100){
  load(paste0("./multi_ethnic/result/y_",1))
  
  ####regression on the training,testing and validation dataset
  ####the ratio between training,testing and validation is 10:1:1
  ####record the summary level statistics
  
  n.snp <- ncol(genotype[[1]])
  beta_summary_train <- matrix(0,n.snp,9)
  beta_summary_test <- matrix(0,n.snp,9)
  beta_summary_vad <- matrix(0,n.snp,9)
  
  
  colnames(beta_summary_train) <- c("beta_EUR","sd_EUR","p_EUR",
                                    "beta_AFR","sd_AFR","p_AFR", "beta_LAC","sd_LAC","p_LAC")
  colnames(beta_summary_vad) = colnames(beta_summary_test) = colnames(beta_summary_train)
  library(RcppArmadillo)
  
  FitLinearmodel <- function(y,x){
    model <- fastLm(X=cbind(1,x),y=y)
    if(is.na(coef(model)[2])){
      result <- c(0,1,1)
    }else{
      result <- coef(summary(model))[2,c(1,2,4)]  
    }
    return(result)
  }
  Fitmodelall <- function(y,G,ind){
    result.train <- NULL
    result.test <- NULL
    result.vad <- NULL
    #i from 1 to 3 represent EUR, AFR, LAC
    for(i in 1:3){
      
      n.sub <- nrow(y[[i]])
      n.train <- n.sub*10/12
      n.test <- n.sub/12
      n.vad <- n.sub/12
      result.i <- NULL
      result.temp <- FitLinearmodel(y[[i]][1:n.train],
                                    G[[i]][1:n.train,ind])
      result.train <- c(result.train,result.temp)
      result.temp <- FitLinearmodel(y[[i]][(1:n.test)+n.train],
                                    G[[i]][(1:n.test)+n.train,ind])
      result.test <- c(result.test,result.temp)
      result.temp <- FitLinearmodel(y[[i]][(1:n.vad)+n.test+n.train],
                                    G[[i]][(1:n.vad)+n.test+n.train,ind])
      result.vad <- c(result.vad,result.temp)
    }
    return(list(result.train,result.test,result.vad))
    
    
  }
  
  for(i in 1:n.snp){
    if(i%%100==0){
      print(i)
    }
    temp_result <- Fitmodelall(y,genotype,i)
    beta_summary_train[i,] <- temp_result[[1]]
    beta_summary_test[i,] <- temp_result[[2]]
    beta_summary_vad[i,] <- temp_result[[2]]
  }
  
  
  beta_result <- list(beta_summary_train,
                      beta_summary_test,
                      beta_summary_vad)
  save(beta_result,file = paste0("./multi_ethnic/result/pruned_geno/beta_estimate_",k,"_",i1,".Rdata"))
  }











