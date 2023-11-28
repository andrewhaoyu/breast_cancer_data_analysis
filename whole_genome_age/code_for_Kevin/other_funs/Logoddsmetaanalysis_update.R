#two studies (e.g. icogs and oncoarray)
LogoddsMetaAnalysis <- function(logodds1,sigma1,logodds2,sigma2){
  sigma1.inv <- solve(sigma1)
  sigma2.inv <- solve(sigma2)
  sigma.meta <- solve(sigma1.inv+sigma2.inv)
  logodds.meta <- sigma.meta%*%(sigma1.inv%*%logodds1+sigma2.inv%*%logodds2)
  
  return(list(logodds.meta = logodds.meta,
              sigma.meta = sigma.meta))
  
}


LogoddsMetaAnalysis <- function(logodds1,sigma1,logodds2,sigma2,logodds3,sigma3){
  sigma1.inv <- solve(sigma1)
  sigma2.inv <- solve(sigma2)
  sigma3.inv <- solve(sigma3)
  sigma.meta <- solve(sigma1.inv+sigma2.inv+sigma3.inv)
  logodds.meta <- sigma.meta%*%(sigma1.inv%*%logodds1+sigma2.inv%*%logodds2+
                                  sigma3.inv%*%logodds3)
  
  return(list(logodds.meta = logodds.meta,
              sigma.meta = sigma.meta))
  
}


#four ancestries (e.g. EUR, AFR, AMR, EAS)
LogoddsMetaAnalysis <- function(logodds1,sigma1,logodds2,sigma2,
                                logodds3,sigma3,logodds4,sigma4){
  sigma1.inv <- solve(sigma1)
  sigma2.inv <- solve(sigma2)
  sigma3.inv <- solve(sigma3)
  sigma4.inv <- solve(sigma4)
  sigma.meta <- solve(sigma1.inv+sigma2.inv+sigma3.inv+sigma4.inv)
  logodds.meta <- sigma.meta%*%(sigma1.inv%*%logodds1+sigma2.inv%*%logodds2+
                                  sigma3.inv%*%logodds3+sigma4.inv%*%logodds4 )
  
  return(list(logodds.meta = logodds.meta,
              sigma.meta = sigma.meta))
  
}

#Global test for associaiton using logodds and sigma
GlobalTestForAssoc <- function(logodds,sigma){
  sigma <- as.matrix(sigma)
  df <- length(logodds)
  GTA.stat <- t(logodds)%*%solve(sigma)%*%logodds
  p.value.GTA <- pchisq(as.numeric(GTA.stat),df=df,lower.tail = F)
  places = 3
  power.number <- floor(-log10(p.value.GTA))+places
  ###format the output with three digits in total
  p.value.GTA <- round(p.value.GTA*10^power.number)/(10^power.number)
  
  return(p.value.GTA)
  
}

