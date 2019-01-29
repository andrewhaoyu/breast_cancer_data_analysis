setwd("/spin1/users/zhangh24/breast_cancer_data_analysis/risk_prediction/two_stage_model/result/")
log.odds.meta.triple <- rep(0,205)
log.odds.meta.two.stage.all <- matrix(0,205,5)
sigma.log.odds.two.stage <- matrix(0,205,25)
log.odds.meta.la.all <- matrix(0,205,5)
p.heter.intrinsic <- rep(0,205)
heter.sigma <- rep(0,205)
#heter.sigma2 <- rep(0,205)



load("/spin1/users/zhangh24/breast_cancer_data_analysis/risk_prediction/standard_analysis/result/log.odds.meta.Rdata")

# solve_variance <- function(x,k,lamda){
#   result <- sum(1/(lamda+x))-sum(k^2/(lamda+x)^2)
#   return(result)
# }
# findvariance <- uniroot(solve_variance,
#                       lower=0,upper = 0.5,k=k,lamda=lamda,a=0)
# 
# log.odds <- meta.result[[1]]
# sigma <- meta.result[[2]]
# beta0 <- log.odds.meta[i]
# solve_variance(0,k,lamda)
# solve_variance(0.3,k,lamda)
# heter.variance.estimate2 <- function(log.odds,sigma,beta0){
#   M <- length(log.odds)
#   sigma.eigen <- eigen(sigma)
#   lamda <- sigma.eigen$values
#   P <- t(sigma.eigen$vectors)
#   k <- P%*%(log.odds-beta0)
#   if(solve_variance(0,k,lamda)*solve_variance(0.3,k,lamda)>=0){
#     result=0
#   }else{
#     findvariance <- uniroot(solve_variance,lower=0,upper = 0.3,k=k,lamda=lamda)
#     result <- as.numeric(findvariance$root)
#     
#   }
#   
#   if(result <= 0){
#     result <- 0
#   }
#   return(result)
# }
# 
# 






heter.variance.estimate <- function(log.odds,sigma){
  M <- length(log.odds)
  result <- (sum((log.odds-mean(log.odds))^2)-sum(diag(sigma))+sum(sigma)/M)/(M-1)
  if(result <= 0){
    result <- 0
  }
  return(result)
}

library(bc2)

for(i1 in 1:205){
  print(i1)
  load(paste0("/spin1/users/zhangh24/breast_cancer_data_analysis/risk_prediction/two_stage_model/result/meta.result",i1,".Rdata"))
  load(paste0("/spin1/users/zhangh24/breast_cancer_data_analysis/risk_prediction/two_stage_model/result/log.odds.meta.la",i1,".Rdata"))
  log.odds.meta.triple[i1] <- meta.result[[1]][5]
  log.odds.meta.two.stage.all[i1,] <- meta.result[[1]]
  log.odds.meta.la.all[i1,] <- log.odds.meta.la
  sigma.log.odds.two.stage[i1,] <- as.vector(meta.result[[2]])
  p.heter.intrinsic[i1] <- GlobalTestForHeter(meta.result[[1]],meta.result[[2]],self.design =T)
  heter.sigma[i1] <-     heter.variance.estimate(meta.result[[1]],meta.result[[2]])
  #heter.sigma2[i1] <- heter.variance.estimate2(meta.result[[1]],meta.result[[2]],log.odds.meta[i])
  
}
j <- 1

log.odds.meta <- log.odds.meta.two.stage.all[,j]
var.log.odds <- rep(0,205)
for(i in 1:205){
  var.log.odds[i]<-  diag(matrix(sigma.log.odds.two.stage[i,],5,5))[j]  
}

#var.odds.meta <- diag(meta.result[[2]])
var.odds.meta <- var.log.odds
idx.control <- which(y.pheno.mis2.train[,1]==0)
snp.control <- x.snp.all.train2[idx.control,]
p <- apply(snp.control,2,function(x){sum(x)/(2*length(x))})
2*sum(p*(1-p)*(log.odds.meta^2-var.odds.meta))
new <- c(178:205)
sum(p[new]*(1-p[new])*(log.odds.meta[new]^2-var.odds.meta[new]))/log(2)








save(log.odds.meta.triple,file="/spin1/users/zhangh24/breast_cancer_data_analysis/risk_prediction/two_stage_model/result/log.odds.meta.triple.Rdata")
save(log.odds.meta.two.stage.all,file="/spin1/users/zhangh24/breast_cancer_data_analysis/risk_prediction/two_stage_model/result/log.odds.meta.two.stage.all.Rdata")
save(sigma.log.odds.two.stage,file="/spin1/users/zhangh24/breast_cancer_data_analysis/risk_prediction/two_stage_model/result/sigma.log.odds.two.stage.Rdata")
save(p.heter.intrinsic,file="/spin1/users/zhangh24/breast_cancer_data_analysis/risk_prediction/two_stage_model/result/p.heter.intrinsic.Rdata")
save(heter.sigma,file = "/spin1/users/zhangh24/breast_cancer_data_analysis/risk_prediction/two_stage_model/result/heter.sigma.Rdata")
save(log.odds.meta.la.all,file="/spin1/users/zhangh24/breast_cancer_data_analysis/risk_prediction/two_stage_model/result/log.odds.meta.la.all.Rdata")
#save(heter.sigma2,file= "/spin1/users/zhangh24/breast_cancer_data_analysis/risk_prediction/two_stage_model/result/heter.sigma2.Rdata")
true.false.calculate <- function(prs,test.data){
  idx.true <- which(test.data==1)
  idx.false <- which(test.data==0)
  n <- length(test.data)
  min.prs <- range(prs)[1]
  max.prs <- range(prs)[2]
  cut.point <- seq(from=min.prs,to=max.prs,by=(max.prs-min.prs)/1000)
  true.pos <- rep(0,length(cut.point))
  false.pos <- rep(0,length(cut.point))
  true.pos[length(cut.point)] <- 1
  false.pos[length(cut.point)] <- 1
  for(i in 2:(length(cut.point)-1)){
    
    predict.result <- ifelse(prs<cut.point[i],1,0)
    
    temp <- table(predict.result,test.data)
    true.pos[i] <- temp[2,2]/colSums(temp)[2]
    false.pos[i] <- (temp[2,1])/colSums(temp)[1]
    
  }
  return(cbind(true.pos,false.pos))
}

auc_cal <- function(roc){
  n <- nrow(roc.result)
  auc <- 0
  for(i in 1:(n-1)){
    temp <- (roc[i+1,1]-roc[i,1])*(roc[i+1,2]+roc[i,2])/2
    auc <- temp+auc
  }
  return(auc)
}


library(data.table)
icog.data <- as.data.frame(fread("/spin1/users/zhangh24/breast_cancer_data_analysis/data/sig_snps_icog.csv",header=T))
onco.data <- as.data.frame(fread("/spin1/users/zhangh24/breast_cancer_data_analysis/data/sig_snps_onco.csv",header=T))
library(tidyverse)
y.pheno.mis1 <- select(icog.data,Behaviour1,ER_status1,PR_status1,HER2_status1,Grade1)
x.covar1 <- select(icog.data,5:14)
x.snp.all1 <- select(icog.data,26:230)
colnames(y.pheno.mis1)
y.pheno.mis2 <- select(onco.data,Behaviour1,ER_status1,PR_status1,HER2_status1,Grade1)
x.covar2 <- select(onco.data,5:14)
x.snp.all2 <- select(onco.data,26:230)
colnames(y.pheno.mis2)

idx.control1 <- which(icog.data$Behaviour1==0)
length(idx.control1)
idx.triple1 <- which(y.pheno.mis1[,2]==0&y.pheno.mis1[,3]==0&y.pheno.mis1[,4]==0)
length(idx.triple1)
idx.control2 <- which(onco.data$Behaviour1==0)
length(idx.control2)
idx.triple2 <- which(y.pheno.mis2[,2]==0&y.pheno.mis2[,3]==0&y.pheno.mis2[,4]==0)
length(idx.triple2)
#############random sample 200 cases & 200 controls from icog
#############random sample 500 cases & 500 controls from onco

n.test.control.icog <- 200
n.test.cases.icog <- 200
n.test.control.onco <- 500
n.test.cases.onco <- 500
set.seed(1)
idx.test.control.icog <- idx.control1[sample(length(idx.control1),n.test.control.icog)]
idx.test.triple.icog <- idx.triple1[sample(length(idx.triple1),n.test.cases.icog)]
idx.test1 <- c(idx.test.control.icog,idx.test.triple.icog)
y.pheno.mis1.test <- y.pheno.mis1[idx.test1,]
x.covar.test1 <- x.covar1[idx.test1,]
x.snp.all.test1 <- x.snp.all1[idx.test1,]


y.pheno.mis1.train <- y.pheno.mis1[-idx.test1,]
x.covar.train1 <- x.covar1[-idx.test1,]
x.snp.all.train1 <- x.snp.all1[-idx.test1,]




idx.test.control.onco <- idx.control2[sample(length(idx.control2),n.test.control.onco)]
idx.test.triple.onco <- idx.triple2[sample(length(idx.triple2),n.test.cases.onco)]
idx.test2 <- c(idx.test.control.onco,idx.test.triple.onco)
y.pheno.mis2.test <- y.pheno.mis2[idx.test2,]
x.covar.test2 <- x.covar2[idx.test2,]
x.snp.all.test2 <- x.snp.all2[idx.test2,]


y.test <- rbind(y.pheno.mis1.test,y.pheno.mis2.test)
x.snp.all.test <- rbind(as.matrix(x.snp.all.test1),as.matrix(x.snp.all.test2))

y.pheno.mis2.train <- y.pheno.mis2[-idx.test2,]
x.covar.train2 <- x.covar2[-idx.test2,]
x.snp.all.train2 <- x.snp.all2[-idx.test2,]


load("/spin1/users/zhangh24/breast_cancer_data_analysis/risk_prediction/two_stage_model/result/log.odds.meta.two.stage.Rdata")
prs <- x.snp.all.test%*%log.odds.meta.two.stage
min.prs <- range(prs)[1]
max.prs <- range(prs)[2]
cut.data <- seq(from=min.prs,to=max.prs,by=(max.prs-min.prs)/100)

roc.result <- true.false.calculate(prs,y.test[,1])
plot(roc.result[,1],roc.result[,2],xlab="false_p",ylab="true_p")
abline(a=0,b=1,col="red")
auc <- auc_cal(roc.result)

