#Goal: implement three prediction methods
#method1: LD pruning and threshold
#method2: weighted combination of the three
#method3: E-Bayes
arg <- commandArgs(trailingOnly=T)
pop.ind <- as.numeric(arg[[1]])
setwd('/spin1/users/zhangh24/breast_cancer_data_analysis')
load(paste0("./multi_ethnic/result/y_",1))
y_all = y


#find all genotypes data files cutpoint
# all.cut <- rep(0,500)
# n.snp <- 0
# for(i1 in 1:500){
#   if(i1%%50==0){
#     print(i1)  
#   }
#   load(paste0("./multi_ethnic/result/pruned_geno/beta_estimate_",i1,".Rdata"))
#  temp <- nrow(beta_result[[1]])
#   n.snp <- n.snp +temp
#   all.cut[i1] <- n.snp
# }
load(paste0("./multi_ethnic/result/pruned_geno/beta_all_",1,".Rdata"))




#pop.ind population indicator 1 EUR, 2 AFR, 3 LAC
LDP <- function(y_all,beta.train,beta.test,beta.vad,p.train,
                p.thr,pop.ind){
  y = y_all[[pop.ind]]
  n.sub <- length(y)
  n.train <- n.sub*10/12
  n.test <- n.sub/12
  n.vad <- n.sub/12
  y.test <- y[n.train+(1:n.test)]
  y.vad <- y[(n.train+n.test)+(1:n.vad)]
  r2.vad <- r2.test <- rep(0,length(p.thr))
  prop <- rep(0,length(p.thr))
  n.snp.sec <- rep(0,length(p.thr))
  prs.mat <- matrix(0,n.test+n.vad,length(p.thr))
  
  
  n.snp <- 0
  for(i in 1:500){
    print(i)
    load(paste0("./multi_ethnic/result/pruned_geno/geno_",i))
    file.snp <- ncol(genotype[[pop.ind]])
    for(j in 1:file.snp){
      jdx <- which(p.train[n.snp+j]<=p.thr)
        if(length(jdx)!=0){
          prs.temp <- genotype[[pop.ind]][,j]*beta.train[n.snp+j]
          prs.mat[,jdx] <- prs.mat[,jdx]+
          prs.temp[n.train+(1:(n.test+n.vad))]
        }
    }
    n.snp <- n.snp+file.snp
  }

  for(k in 1:length(p.thr)){
    idx <- which(p.train<=p.thr[k])
    if(length(idx)==0){
      n.snp.sec[k] <- 0
      prop[k] <- 0
      r2.test[k] <- 0
      r2.vad[k] <- 0
    }else{
      n.snp.sec[k] <- length(idx)
      cau.snp <- c(1:4000,4000+c(1:1000)+1000*(pop.ind-1))
      prop[k] <- sum(idx%in%cau.snp)/length(idx)
      prs.test <- prs.mat[(1:n.test),k]
      prs.vad <- prs.mat[n.test+(1:n.vad),k]
      model1 <- lm(y.test~prs.test)
      r2.test[k] <- summary(model1)$adj.r.squared
      model2 <- lm(y.vad~prs.vad)
      r2.vad[k] <- summary(model2)$adj.r.squared
    }
  }

  result <- list(n.snp.sec,prop,
                 r2.test,r2.vad,
                 prs.mat)
  return(result)
}  
p.thr <- c(10^-8,10^-7,10^-6,10^-5,10^-4,10^-3,10^-2,0.1,0.3,0.5)
beta.train <- beta_result[[1]][,3*pop.ind-2]
beta.test <- beta_result[[2]][,3*pop.ind-2]
beta.vad <- beta_result[[3]][,3*pop.ind-2]
p.train <- beta_result[[1]][,3*pop.ind]

LDP.result <-  LDP(y,beta.train,beta.test,beta.vad,p.train,p.thr,pop.ind)
save(LDP.result,file = paste0("./multi_ethnic/result/LDP.result_",pop.ind))
  
# try <- genotype[[pop.ind]][,idx]%*%beta.train[idx]
# try2 <- try[n.train+(1:(n.test+n.vad))]
# 
#   
# LDP <- function(y_all,beta.train,beta.test,beta.vad,p.train,
#                 p.thr,pop.ind){
#   y = y_all[[pop.ind]]
#   n.sub <- length(y)
#   n.train <- n.sub*10/12
#   n.test <- n.sub/12
#   n.vad <- n.sub/12
#   y.test <- y[n.train+(1:n.test)]
#   y.vad <- y[(n.train+n.test)+(1:n.vad)]
#   r2.vad <- r2.test <- rep(0,length(p.thr))
#   prop <- rep(0,length(p.thr))
#   n.snp.sec <- rep(0,length(p.thr))
#   prs.mat <- matrix(0,n.test+n.vad,length(p.thr))
#   
#   for(k in 1:length(p.thr)){
#     #print(k)
#     idx <- which(p.train<=p.thr[k])
#     if(length(idx)==0){
#       n.snp.sec[k] <- 0
#       prop[k] <- 0
#       r2.test[k] <- 0
#       r2.vad[k] <- 0
#     }else{
#       n.snp.sec[k] <- length(idx)
#       cau.snp <- c(1:4000,4000+c(1:1000)+1000*(pop.ind-1))
#       prop[k] <- sum(idx%in%cau.snp)/length(idx)
#       prs <- rep(0,n.sub)
#       temp.result <- Findfilenum(idx,all.cut)
#       filenum <- temp.result[[1]]
#       colnum <- temp.result[[2]]
#       if(is.null(filenum)){
#         r2.test[k] <- NULL
#         r2.vad[k] <- NULL
#       }else{
#         
#         for(i in 1:length(filenum)){
#           print(i)
#           if(i==1){
#             #####load the first genotype file
#             file.temp <- filenum[i]
#             load(paste0("./multi_ethnic/result/pruned_geno/geno_",file.temp))
#           }else if(file.temp!=filenum[i]){
#             #avoid reload genotype data
#             load(paste0("./multi_ethnic/result/pruned_geno/geno_",filenum[i]))
#             file.temp <- filenum[i]
#           }
#           #take the corresponding column number for selected SNP
#           geno <- genotype[[pop.ind]][,colnum[i]]
#           prs.temp <- beta.train[idx[i]]*geno
#           prs <- prs+prs.temp  
#         }
#       }
#       prs.test <- prs[n.train+(1:n.test)]
#       prs.vad <- prs[n.train+n.test+(1:n.vad)]
#       model1 <- lm(y.test~prs.test)
#       r2.test[k] <- summary(model1)$adj.r.squared
#       model2 <- lm(y.vad~prs.vad)
#       r2.vad[k] <- summary(model2)$adj.r.squared
#       prs.mat[,k] <- c(prs.test,prs.vad)
#     }
#     
#   }
#   result <- list(n.snp.sec,prop,
#                  r2.test,r2.vad,
#                  prs.mat)
#   return(result)
#   
# }
# 

  


# #find the file number and column number for significant SNP
# Findfilenum <- function(idx, all.cut){
#   n.idx <- length(idx)
#   if(n.idx==0){
#     filenum = NULL
#     colnum = NULL
#   }else{
#     filenum <- rep(0,n.idx)
#     colnum <- rep(0,n.idx)
#     for(j in 1:length(idx)){
#       for(i in 1:500){
#         
#         if(idx[j] <=all.cut[i]){
#           filenum[j] <- i
#           if(i==1){
#             colnum[j] <- idx[j]
#           }else{
#             colnum[j] <- idx[j]-all.cut[i-1]
#           }
#           break;}
#         }
#     }
#   }
#   return(list(filenum,colnum))
# }


