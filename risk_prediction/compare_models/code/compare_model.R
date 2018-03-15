auc_cal <- function(roc){
  n <- nrow(roc)
  auc <- 0
  for(i in 1:(n-1)){
    temp <- (roc[i+1,1]-roc[i,1])*(roc[i+1,2]+roc[i,2])/2
    auc <- temp+auc
  }
  return(auc)
}


# prior_sigma <- function(log.odds,sigma){
#   p <- ncol(sigma)
#   n <- nrow(log.odds)
#   mean.log.odds <- apply(log.odds,1,mean)
#   n.eff <- 1/sigma
#   n.eff.sum <- apply(n.eff,1,sum)
#   prior.sigma.result <- rep(0,205)
#   
#   for(i in 1:n){
#     temp <- 0
#     for(j in 1:p){
#       temp <- (n.eff[i,j]/n.eff.sum[i])*(log.odds[i,j]-mean.log.odds[i])^2
#       prior.sigma.result[i] <- temp+prior.sigma.result[i]
#     }
#   }
#   prior.sigma.result <- prior.sigma.result/(p-1)
#   return(prior.sigma.result)
# }
true.false.calculate <- function(prs,test.data){
  idx.true <- which(test.data==1)
  idx.false <- which(test.data==0)
  n <- length(test.data)
  min.prs <- range(prs)[1]
  max.prs <- range(prs)[2]
  cut.point <- seq(from=min.prs,to=max.prs,by=(max.prs-min.prs)/100)
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

#############################compare different model
######standard analysis result
load("/spin1/users/zhangh24/breast_cancer_data_analysis/risk_prediction/standard_analysis/result/log.odds.meta.Rdata")
#####intrinsic subtype result
load("/spin1/users/zhangh24/breast_cancer_data_analysis/risk_prediction/two_stage_model/result/log.odds.meta.triple.Rdata")
load("/spin1/users/zhangh24/breast_cancer_data_analysis/risk_prediction/two_stage_model/result/p.heter.intrinsic.Rdata")
load("/spin1/users/zhangh24/breast_cancer_data_analysis/risk_prediction/two_stage_model/result/log.odds.meta.two.stage.all.Rdata")
load("/spin1/users/zhangh24/breast_cancer_data_analysis/risk_prediction/two_stage_model/result/sigma.log.odds.two.stage.Rdata")
load("/spin1/users/zhangh24/breast_cancer_data_analysis/risk_prediction/two_stage_model/result/log.odds.meta.triple.Rdata")
load("/spin1/users/zhangh24/breast_cancer_data_analysis/risk_prediction/two_stage_model/result/heter.sigma.Rdata")
####additive two-stage model
load("/spin1/users/zhangh24/breast_cancer_data_analysis/risk_prediction/two_stage_model/result/p.heter.add.Rdata")
load("/spin1/users/zhangh24/breast_cancer_data_analysis/risk_prediction/two_stage_model/result/log.odds.add.triple.Rdata")
load("/spin1/users/zhangh24/breast_cancer_data_analysis/risk_prediction/two_stage_model/result/log.odds.add.triple.eb.Rdata")






##############standard model
calibration <- function(y,prs){
  n <- length(y)
  n.qun <- 11
  idx.control <- which(y==0)
  idx.case <- which(y==1)
  prs.control <- prs[idx.control]
  prs.case <- prs[idx.case]
  control.qun <- quantile(prs.control,probs=seq(0,1,1/(n.qun-1)))
  odds <- rep(0,n.qun-2)
  for(i in 2:(n.qun-1)){
    if(i==2){
      idx <- which(prs.control<=control.qun[i])
      n.control <- length(idx)
      idx <- which(prs.case<= control.qun[i])
      n.case <- length(idx)
      odds[i-1] <- n.case/n.control
    }else if(i==(n.qun-1)){
      
        idx <- which(prs.control>=control.qun[i])
        n.control <- length(idx)
        idx <- which(prs.case>= control.qun[i])
        n.case <- length(idx)
        odds[n.qun-2] <- n.case/n.control
      }else{
        idx <- which(prs.control>=control.qun[i-1]&
                       prs.control<=control.qun[i])
        n.control <- length(idx)
        idx <- which(prs.case>= control.qun[i-1]&
                       prs.case<= control.qun[i])
        n.case <- length(idx)
        odds[i-1] <- n.case/n.control
      }
    
  }
  return(odds)
}




library(pRoc)
prs.standard <- x.snp.all.test%*%log.odds.meta
cal.standard <- calibration(y.test[,1],prs.standard)
roc.try <- roc(y.test[,1],prs.standard)
roc.result.standard <- true.false.calculate(prs.standard,y.test[,1])
plot(roc.result.standard[,1],roc.result.standard[,2],xlab="false_p",ylab="true_p")
abline(a=0,b=1,col="red")
n <- nrow(roc.result.standard)
lab = rep("standard",n)
auc.standard <- auc_cal(roc.result.standard)
roc.result.standard.mer <- data.frame(roc.result.standard,lab)
#############intrinsic subtypes
prs.intrinsic <- x.snp.all.test%*%log.odds.meta.triple
cal.intrinsic <- calibration(y.test[,1],prs.intrinsic)
roc.result.intrinsic <- true.false.calculate(prs.intrinsic,y.test[,1])
plot(roc.result.intrinsic[,1],roc.result.intrinsic[,2],xlab="false_p",ylab="true_p")
abline(a=0,b=1,col="red")
lab = rep("intrinsic subtypes",n)
auc.intrinsic <- auc_cal(roc.result.intrinsic)
roc.reuslt.intrinsic.mer <- data.frame(roc.result.intrinsic,lab)
###########intrinsic subtypes dic mixed
log.odds.intrinsic.dic <- log.odds.meta
idx <- which(p.heter.intrinsic<0.05)
log.odds.intrinsic.dic[idx] <- log.odds.meta.triple[idx]

prs.intrinsic.dic <- x.snp.all.test%*%log.odds.intrinsic.dic
cal.intrinsic.dic<- calibration(y.test[,1],prs.intrinsic.dic)
roc.result.intrinsic.dic <- true.false.calculate(prs.intrinsic.dic,y.test[,1])
plot(roc.result.intrinsic.dic[,1],roc.result.intrinsic.dic[,2],xlab="false_p",ylab="true_p")
abline(a=0,b=1,col="red")
auc.intrinsic.dic <- auc_cal(roc.result.intrinsic.dic)
lab = rep("intrinsic dicho",n)
roc.result.intrinsic.dic.mer <- data.frame(roc.result.intrinsic.dic,lab)
#########intrinsic subtypes eb 
ebestimate <- function(logodds.subtype,
                       sigma.subtype,
                       logodds.standard,
                       prior.sigma
){
  M <- length(logodds.subtype)
  if(prior.sigma==0){
    return(rep(logodds.standard,M))
  }else{
    result <- solve(solve(sigma.subtype)+(1/prior.sigma)*diag(M))%*%(solve(sigma.subtype)%*%logodds.subtype+
                                                           (1/prior.sigma)*rep(logodds.standard,M))
    return(result)
  }
}

log.odds.intrinsic.eb <- rep(0,205)
for(i in 1:205){
 
   logodds.subtype <- log.odds.meta.two.stage.all[i,]
 M <- length(logodds.subtype)
    sigma.subtype <- matrix(sigma.log.odds.two.stage[i,],5,5)
    log.odds.intrinsic.eb[i] <- ebestimate(logodds.subtype,sigma.subtype,
             as.numeric(log.odds.meta[i]),
             as.numeric(heter.sigma[i]))[5]
}

prs.intrinsic.eb <- x.snp.all.test%*%log.odds.intrinsic.eb
cal.intrinsic.eb <- calibration(y.test[,1],prs.intrinsic.eb)
roc.result.intrinsic.eb <- true.false.calculate(prs.intrinsic.eb ,y.test[,1])
plot(roc.result.intrinsic.eb[,1],roc.result[,2],xlab="false_p",ylab="true_p")
abline(a=0,b=1,col="red")
auc.intrinsic.eb <- auc_cal(roc.result.intrinsic.eb)
lab = rep("intrinsic eb",n)
roc.result.intrinsic.eb <- data.frame(roc.result.intrinsic.eb,lab)






#######additive two-stage model
prs.add <- x.snp.all.test%*%log.odds.add.triple
cal.add <- calibration(y.test[,1],prs.add)
roc.result.add <- true.false.calculate(prs.add,y.test[,1])
plot(roc.result.add[,1],roc.result.add[,2],xlab="false_p",ylab="true_p")
abline(a=0,b=1,col="red")
auc.add <- auc_cal(roc.result.add)
lab = rep("add",n)
roc.result.add.mer <- data.frame(roc.result.add,lab)
#######additive two-stage model dichotomized
log.odds.add.dic <- log.odds.meta

idx <- which(p.heter.add<0.05)
log.odds.add.dic[idx] <- log.odds.add.triple[idx]

prs.add.dic <- x.snp.all.test%*%log.odds.add.dic
cal.add.dic <- calibration(y.test[,1],prs.add.dic)
roc.result.add.dic <- true.false.calculate(prs.add.dic,y.test[,1])
plot(roc.result.add.dic[,1],roc.result.add.dic[,2],xlab="false_p",ylab="true_p")
abline(a=0,b=1,col="red")
auc.add.dic <- auc_cal(roc.result.add.dic)
lab = rep("add dic",n)
roc.result.add.dic <- data.frame(roc.result.add.dic,lab)
#######additive two-stage model eb

prs.add.eb <- x.snp.all.test%*%log.odds.add.triple.eb
cal.add.eb <- calibration(y.test[,1],prs.add.eb)
roc.result.add.eb <- true.false.calculate(prs.add.eb,y.test[,1])
plot(roc.result.add.eb[,1],roc.result.add.eb[,2],xlab="false_p",ylab="true_p")
abline(a=0,b=1,col="red")
auc.add.eb <- auc_cal(roc.result.add.eb)
lab = rep("add.eb",n)
roc.result.add.eb.mer <- data.frame(roc.result.add.eb,lab)


all.data <- rbind(roc.result.standard.mer,roc.reuslt.intrinsic.mer,roc.result.intrinsic.dic.mer,roc.result.intrinsic.eb.pin.mer,roc.result.add.mer,roc.result.add.dic,roc.result.add.eb.mer)



cal.result <- rbind(cal.standard,cal.intrinsic,cal.intrinsic.dic,cal.intrinsic.eb,cal.add,cal.add.dic,cal.add.eb)
row.names(cal.result) <- c(
  "standard analysis",
  "intrinsic subtype",
  "intrinsic subtype dichotomized",
  "intrinsic subtype empirical bayesian",
  "additive two-stage model",
  "additive two-stage model dichotomized",
  "additive two-stage model empirical bayesian"
)
write.csv(cal.result,file="/spin1/users/zhangh24/breast_cancer_data_analysis/risk_prediction/two_stage_model/result/cal.result.csv")










sigma.heter <- apply(log.odds.meta.two.stage.all,1,var)
log.odds.meta.mean <- apply(log.odds.meta.two.stage,1,mean)
sigma.heter.improve <- rep(0,205)








log.odds.mix <- log.odds.meta
idx <- which(p.heter.add<0.05)
log.odds.mix[idx] <- log.odds.meta.two.stage[idx]

prs <- x.snp.all.test%*%log.odds.mix
min.prs <- range(prs)[1]
max.prs <- range(prs)[2]
cut.data <- seq(from=min.prs,to=max.prs,by=(max.prs-min.prs)/100)

roc.result <- true.false.calculate(prs,y.test[,1])
plot(roc.result[,1],roc.result[,2],xlab="false_p",ylab="true_p")
abline(a=0,b=1,col="red")
auc <- auc_cal(roc.result)


prior.sigma.result <- prior.sigma(log.odds.meta.two.stage.all,sigma.log.odds.two.stage)





rho <- rep(0,205)
for(i in 1:205){
  rho[i] <- prior.sigma.result[i]/(prior.sigma.result[i]+sigma.log.odds.two.stage[i,5])
}

log.odds.bayes <- rho*log.odds.meta.two.stage+(1-rho)*log.odds.meta
prs <- x.snp.all.test%*%log.odds.bayes
min.prs <- range(prs)[1]
max.prs <- range(prs)[2]
cut.data <- seq(from=min.prs,to=max.prs,by=(max.prs-min.prs)/100)

roc.result <- true.false.calculate(prs,y.test[,1])
plot(roc.result[,1],roc.result[,2],xlab="false_p",ylab="true_p")
abline(a=0,b=1,col="red")
auc <- auc_cal(roc.result)



rho.addprior.tri <- rep(0,205)
for(i in 1:205){
  rho[i] <- prior.sigma.add[i]/(prior.sigma.add[i]+sigma.log.odds.two.stage[i,5])
}

log.odds.bayes.addprio <- rho*log.odds.meta.two.stage+(1-rho)*log.odds.meta
prs <- x.snp.all.test%*%log.odds.bayes.addprio 
min.prs <- range(prs)[1]
max.prs <- range(prs)[2]
cut.data <- seq(from=min.prs,to=max.prs,by=(max.prs-min.prs)/100)

roc.result <- true.false.calculate(prs,y.test[,1])
plot(roc.result[,1],roc.result[,2],xlab="false_p",ylab="true_p")
abline(a=0,b=1,col="red")
auc <- auc_cal(roc.result)



