mine <- c(1:6)
table.combos <- matrix(data = 1:18, nrow = 10, ncol = 6, byrow=T)
apply(table.combos,1,function(x){identical(mine,x)})


ProbFitting <- function(delta0,y,x.all,z.standard,z.all,missingTumorIndicator=NULL){
  idx.drop = NULL
  if(is.null(missingTumorIndicator)==1){
    n <- nrow(y)
    M <- nrow(z.standard)
    y_em <- matrix(0,nrow=n,ncol = M)
    beta <- matrix(z.all%*%delta0,ncol = M)
    ##add intercept to x.all
    x.all.inter <- as.matrix(cbind(1,x.all))
    index = 1
    for(i in 1:nrow(y)){
      if(y[i,1]==1) {
        ###find out which tumor characteristic is observed
        ###-3.14 is just a random number to make the algorithm run
        ###since there is no missing, all of the 2:ncol(y) will be chose
        idx <- 1:(ncol(y)-1)
        ###jdx is the potential subtype this missing person could be
        jdx <- apply(z.standard,1,function(t){all(t[idx]==y[i,idx+1])})
        if(sum(jdx)==0){
          idx.drop = c(idx.drop,i)
        }else{
          jdx <- which(jdx==T)
          ####get the conditional probability
          y_em[i,jdx] <- 1  
        }
        
      }
    }
    return(list(y_em=y_em,missing.vec = NULL , missing.mat = NULL,complete.vec = NULL,
                idx.drop = idx.drop))
  }else{
    n <- nrow(y)
    M <- nrow(z.standard)
    y_em <- matrix(0,nrow=n,ncol = M)
    missing.vec = rep(0,n)
    missing.mat = matrix(0,nrow=n,ncol=M)
    
    beta <- matrix(z.all%*%delta0,ncol = M)
    ##add intercept to x.all
    x.all.inter <- as.matrix(cbind(1,x.all))
    index = 1
    
    
    for(i in 1:nrow(y)){
      if(y[i,1]==1) {
        ###find out which tumor characteristic is observed
        idx <- which(y[i,2:ncol(y)]!=missingTumorIndicator)
        ###jdx is the potential subtype this missing person could be
        jdx <- apply(z.standard,1,function(t){all(t[idx]==y[i,idx+1])})
        if(sum(jdx)>=1){
          missing.vec[index] <- i
          missing.mat[index,] <- jdx
          index <- index + 1
          jdx <- which(jdx==T)
          ####get the conditional probability
          temp <- exp(x.all.inter[i,]%*%beta[,jdx])
          y_em[i,jdx] <- temp/sum(temp)
        }else if(sum(jdx)==0){
          idx.drop = c(idx.drop,i)
        }else{
          y_em[i,jdx] <- 1
        }
        
      }
    }
    
    missing.vec <- missing.vec[1:(index-1)]
    missing.mat <- missing.mat[1:(index-1),]
    complete.vec <- c(1:n)
    complete.vec <- complete.vec[!(complete.vec%in%missing.vec)]
    
    return(list(y_em=y_em,missing.vec = missing.vec , missing.mat = missing.mat,complete.vec = complete.vec,
                idx.drop = idx.drop))
  }
  
}






#remove all the subtypes below the cutoff threshold
idx.sub <- which(freq.subtypes[,ncol(freq.subtypes)]<=cutoff)

if(length(idx.sub)!=0){

  
  idx.cut <- NULL
  for(i in 1:length(idx.sub)){
    if(freq.subtypes[idx.sub[i],1:(ncol(freq.subtypes)-1)]==0){
      idx.temp <- NULL
    }else{
      idx.temp <-which(
        apply(y[,2:ncol(y),drop=F],1,function(x){
          identical(x,freq.subtypes[idx.sub[i],1:(ncol(freq.subtypes)-1)])
        })
        )
    }
    idx.cut <- c(idx.cut,idx.temp)
  }
}
if(length(idx.cut)!=0){
  y = y[-idx.cut,drop=F]
  if(is.null(baselineonly)!=0){
    baselineonly = baselineonly[-idx.cut,drop=F]
  }
  if(is.null(additive)!=0){
    additive = additive[-idx.cut,drop=F]
  }
  if(is.null(pairwise.interaction)!=0){
    pairwise.interaction = pairwise.interaction[-idx.cut,drop=F]
  }
  if(is.null(saturated)!=0){
    saturated = saturated[-idx.cut,drop=F]
  }
}



y <- data[,1:5]
SNP <- data[,6,drop=F]
PC1 <- data[,7,drop=F]

baselineonly=NULL;
pairwise.interaction=NULL;
additive = PC1;
saturated=NULL;
missingTumorIndicator <- 888
missing.data.vec <- GenerateMissingPosition(y,missingTumorIndicator)
y.pheno.complete <- y[-missing.data.vec,,drop=F]
initial.set <- InitialSetup(y.pheno.complete,
                            baselineonly,
                            additive,
                            pairwise.interaction,
                            saturated
)
###z standard matrix means the additive model z design matrix without baseline effect
###z standard matrix is used to match the missing tumor characteristics to the complete subtypes
if(is.null(delta0)==T){
  delta0 = initial.set$delta0
}else{
  delta0 = delta0
}

z.all = initial.set$z.all
z.standard = initial.set$z.standard
z.deisign.baselineonly = initial.set$z.design.baseline.only
z.design.additive = initial.set$z.design.additive
z.design.pairwise.interaction = initial.set$z.design.pairwise.interaction
z.design.saturated = initial.set$z.design.saturated
x.all <- as.matrix(GenerateXAll(y,baselineonly,additive,pairwise.interaction,saturated))
covar.names <- initial.set$covar.names
tumor.names <- initial.set$tumor.names
tol <- as.numeric(1e-04)
tolMaxstep <- as.numeric(1e-03)
#delta_old <- rep(0,length(delta0))
delta_old <- delta0
##EM algorithm
##first E step
print(paste0("Begin EM algorithm"))
print(paste0("EM round: 1"))
prob.fit.result <- ProbFitting(delta_old,y,x.all,
                               z.standard,z.all,missingTumorIndicator)
y_em <- prob.fit.result[[1]]
idx <- which(y[,2]==0&y[,3]==0&y[,4]==0&y[,5]==0)
y_em[idx,]
