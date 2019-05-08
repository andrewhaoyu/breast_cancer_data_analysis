#goal: implement the dynamic p-value methods
PostBeta <- function(beta,Sigma,Sigma0){
  n <- length(beta)
  beta_post <- solve(Sigma+solve(Sigma0))%*%(Sigma%*%beta)
  return(beta_post)
}


LDPDyW <- function(y_all,
                  beta.train,
                  sd.train,
                  p.train,
                  p.thr,
                  pop.ind,
                  beta.ref,
                  sd.ref,
                  p.ref,
                  alpha){
  y = y_all[[pop.ind]]
  n.sub <- length(y)
  n.train <- n.sub*10/12
  n.test <- n.sub/12
  n.vad <- n.sub/12
  y.test <- y[n.train+(1:n.test)]
  y.vad <- y[(n.train+n.test)+(1:n.vad)]
  n.cat <- length(alpha)*length(p.thr)
  r2.vad <- r2.test <- rep(0,n.cat)
  
  prop <- rep(0,n.cat)
  n.snp.sec <- rep(0,n.cat)
  alpha.vec <- rep(0,n.cat)
  pthr.vec <- rep(0,n.cat)
  prs.mat <- matrix(0,n.test+n.vad,n.cat)
  
  ######create a p value matrix to contain all different alpha
  p.mat <- matrix(0,length(p.train),length(alpha))
  temp <- 1
  for(k in 1:length(alpha)){
      p.mat[,k] <-   apply(cbind(p.train*alpha[k],p.ref),1,min)
  }
  #####create posterior beta for different calibration point
  beta.train.update<- matrix(0,length(beta.train),
                              n.cat)
  beta.all <- cbind(beta.train,beta.ref)
  sd.all <- cbind(sd.train,sd.ref)
  temp <- 1
  for(k in 1:length(alpha)){
    print(k)
    for(j in 1:length(p.thr)){
      print(j)
      #find the  corresponding updated p
      p.update <- p.mat[,k]
      #find the corresponding snps index passing thres
      idx <- which(p.update<=p.thr[j])
      Sigma0 <- cov(beta.all[idx,])
      for(i in 1:length(idx)){
        Sigma <- diag(sd.all[idx[i],]^2)
        beta.train.update[idx[i],temp]  <- PostBeta(beta.all[idx[i],],
                                                    Sigma,
                                                    Sigma0
        )[1]
      }
      temp <- temp+1
    }
    
  }
  
  
  n.snp <- 0
  
  #p.update <- apply(cbind(p.train*alpha[k],p.ref),1,min)
  for(i in 1:500){
    print(i)
    load(paste0("./multi_ethnic/result/pruned_geno/geno_",i))
    file.snp <- ncol(genotype[[pop.ind]])
    for(k in 1:length(alpha)){
      p.update <- p.mat[,k]
      for(j in 1:file.snp){
        jdx <- which(p.update[n.snp+j]<=p.thr)
        if(length(jdx)!=0){
          #different p.value thr (jdx) get different beta
          #different alpha (k) get differnt beta
          prs.temp <- matrix(0,n.test+n.vad,length(jdx))
          for(l in 1:length(jdx)){
            prs.temp[,l] <- genotype[[pop.ind]][n.train+(1:(n.test+n.vad)),j]*beta.train.update[n.snp+j,jdx+(k-1)*length(p.thr),drop=F][,l]  
          }
          
            prs.mat[,jdx+(k-1)*length(p.thr)] <- prs.mat[,jdx+(k-1)*length(p.thr)]+
            prs.temp
          }
      }
      }
    n.snp <- n.snp+file.snp
  }
  
  temp <- 1
  for(j in 1:length(alpha)){
    p.update <- p.mat[,j]
    for(k in 1:length(p.thr)){
      idx <- which(p.update<=p.thr[k])
      if(length(idx)==0){
        n.snp.sec[temp] <- 0
        prop[temp] <- 0
        r2.test[temp] <- 0
        r2.vad[temp] <- 0
        alpha.vec[temp] <- alpha[j]
        pthr.vec[temp] <- p.thr[k]
        temp <- temp+1
      }else{
        n.snp.sec[temp] <- length(idx)
        cau.snp <- c(1:4000,4000+c(1:1000)+1000*(pop.ind-1))
        prop[temp] <- sum(idx%in%cau.snp)/length(idx)
        prs.test <- prs.mat[(1:n.test),temp]
        prs.vad <- prs.mat[n.test+(1:n.vad),temp]
        model1 <- lm(y.test~prs.test)
        r2.test[temp] <- summary(model1)$adj.r.squared
        model2 <- lm(y.vad~prs.vad)
        r2.vad[temp] <- summary(model2)$adj.r.squared
        alpha.vec[temp] <- alpha[j]
        pthr.vec[temp] <- p.thr[k]
        temp <- temp+1
      }
    }
  }

  result <- data.frame(pthr.vec,
                       alpha.vec,
                       r2.test,
                       r2.vad,
                       n.snp.sec,
                       prop)
  
  return(result)
}  






LDPDy <- function(y_all,
                  beta.train,
                  p.train,
                  p.thr,
                  pop.ind,
                  beta.ref,
                  p.ref,
                  alpha){
  y = y_all[[pop.ind]]
  n.sub <- length(y)
  n.train <- n.sub*10/12
  n.test <- n.sub/12
  n.vad <- n.sub/12
  y.test <- y[n.train+(1:n.test)]
  y.vad <- y[(n.train+n.test)+(1:n.vad)]
  n.cat <- length(alpha)*length(p.thr)
  r2.vad <- r2.test <- rep(0,n.cat)
  
  prop <- rep(0,n.cat)
  n.snp.sec <- rep(0,n.cat)
  alpha.vec <- rep(0,n.cat)
  pthr.vec <- rep(0,n.cat)
  prs.mat <- matrix(0,n.test+n.vad,n.cat)
  
  ######create a p value matrix to contain all different situation
  p.mat <- matrix(0,length(p.train),length(alpha))
  temp <- 1
  for(k in 1:length(alpha)){
    p.mat[,k] <-   apply(cbind(p.train*alpha[k],p.ref),1,min)
  }
  n.snp <- 0
  
  
  #p.update <- apply(cbind(p.train*alpha[k],p.ref),1,min)
  for(i in 1:500){
    print(i)
    load(paste0("./multi_ethnic/result/pruned_geno/geno_",i))
    file.snp <- ncol(genotype[[pop.ind]])
    for(k in 1:length(alpha)){
      p.update <- p.mat[,k]
      for(j in 1:file.snp){
        jdx <- which(p.update[n.snp+j]<=p.thr)
        if(length(jdx)!=0){
          prs.temp <- genotype[[pop.ind]][,j]*beta.train[n.snp+j]
          prs.mat[,jdx+(k-1)*length(p.thr)] <- prs.mat[,jdx+(k-1)*length(p.thr)]+
            prs.temp[n.train+(1:(n.test+n.vad))]
        }
      }
    }
    n.snp <- n.snp+file.snp
  }
  
  temp <- 1
  for(j in 1:length(alpha)){
    p.update <- p.mat[,j]
    for(k in 1:length(p.thr)){
      idx <- which(p.update<=p.thr[k])
      if(length(idx)==0){
        n.snp.sec[temp] <- 0
        prop[temp] <- 0
        r2.test[temp] <- 0
        r2.vad[temp] <- 0
        alpha.vec[temp] <- alpha[j]
        pthr.vec[temp] <- p.thr[k]
        temp <- temp+1
      }else{
        n.snp.sec[temp] <- length(idx)
        cau.snp <- c(1:4000,4000+c(1:1000)+1000*(pop.ind-1))
        prop[temp] <- sum(idx%in%cau.snp)/length(idx)
        prs.test <- prs.mat[(1:n.test),temp]
        prs.vad <- prs.mat[n.test+(1:n.vad),temp]
        model1 <- lm(y.test~prs.test)
        r2.test[temp] <- summary(model1)$adj.r.squared
        model2 <- lm(y.vad~prs.vad)
        r2.vad[temp] <- summary(model2)$adj.r.squared
        alpha.vec[temp] <- alpha[j]
        pthr.vec[temp] <- p.thr[k]
        temp <- temp+1
      }
    }
  }
  
  result <- data.frame(pthr.vec,
                       alpha.vec,
                       r2.test,
                       r2.vad,
                       n.snp.sec,
                       prop)
  
  return(result)
}  


arg <- commandArgs(trailingOnly=T)
i1 <- as.numeric(arg[[1]])
pop.ind <- as.numeric(arg[[2]])
i3 <- as.numeric(arg[[3]])

setwd('/spin1/users/zhangh24/breast_cancer_data_analysis')
load(paste0("./multi_ethnic/result/pruned_geno/beta_all_",i1,".Rdata"))
load(paste0("./multi_ethnic/result/y_",i1))
y_all <- y
beta_train <- beta_result[[1]]
sd_train1 <- beta_train[,2]
p_train1 <- beta_train[,3]
sd_train2 <- beta_train[,5]
p_train2 <- beta_train[,6]
sd_train3 <- beta_train[8]
p_train3 <- beta_train[,9]
beta_train1 <- beta_train[,1]
beta_train2 <- beta_train[,4]
beta_train3 <- beta_train[,7]

p.thr <- c(10^-8,10^-7,10^-6,10^-5,10^-4,10^-3,10^-2,0.1,0.3,0.5)
alpha <- c(1,0.5,0.3,0.15,0.1,0.01,10^-3)
beta.ref <- beta_train1
p.ref <- p_train1
sd.ref <- sd_train1
p.train <- beta_train[,3*pop.ind]
beta.train <- beta_train[,3*pop.ind-2]
sd.train <- beta_train[,3*pop.ind-1]
if(i3==1){
  result <- LDPDy(y_all,
        beta.train,
        p.train,
        p.thr,
        pop.ind,
        beta.ref,
        p.ref,
        alpha)
save(result,file=paste0("./multi_ethnic/result/Dy_result_",i1,"_",pop.ind,"_",i3))  
  
}else{
  result <- LDPDyW(y_all,
                     beta.train,
                     sd.train,
                     p.train,
                     p.thr,
                     pop.ind,
                     beta.ref,
                     sd.ref,
                     p.ref,
                     alpha)
  save(result,file=paste0("./multi_ethnic/result/Dy_result_",i1,"_",pop.ind,"_",i3))  
}            

