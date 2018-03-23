setwd('/spin1/users/zhangh24/breast_cancer_data_analysis/risk_prediction/two_stage_model/result/')
auc_cal <- function(roc){
  n <- nrow(roc)
  auc <- 0
  for(i in 1:(n-1)){
    temp <- (roc[i+1,1]-roc[i,1])*(roc[i+1,2]+roc[i,2])/2
    auc <- temp+auc
  }
  return(auc)
}


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
########subtypes risk prediction
library(bcutility)
library(bc2)
library(data.table)
icog.data <- as.data.frame(fread("/spin1/users/zhangh24/breast_cancer_data_analysis/data/sig_snps_icog.csv",header=T))
onco.data <- as.data.frame(fread("/spin1/users/zhangh24/breast_cancer_data_analysis/data/sig_snps_onco.csv",header=T))
library(tidyverse)
y.pheno.mis1 <- select(icog.data,Behaviour1,ER_status1,PR_status1,HER2_status1,Grade1)

subtypes.icog <- GenerateIntrinsicmis(y.pheno.mis1[,2],y.pheno.mis1[,3],
                                      y.pheno.mis1[,4],y.pheno.mis1[,5])

#table(subtypes.icog)+table(subtypes.onco)

x.covar1 <- select(icog.data,5:14)
x.snp.all1 <- select(icog.data,26:230)
colnames(y.pheno.mis1)

icog.test.id <- Generatetestid(subtypes.icog)

y.pheno.mis2 <- select(onco.data,Behaviour1,ER_status1,PR_status1,HER2_status1,Grade1)
subtypes.onco <- GenerateIntrinsicmis(y.pheno.mis2[,2],
                                      y.pheno.mis2[,3],
                                      y.pheno.mis2[,4],
                                      y.pheno.mis2[,5])
x.covar2 <- select(onco.data,5:14)
x.snp.all2 <- select(onco.data,26:230)
colnames(y.pheno.mis2)
subtypes.onco <- GenerateIntrinsicmis(y.pheno.mis2[,2],y.pheno.mis2[,3],
                                      y.pheno.mis2[,4],y.pheno.mis2[,5])
onco.test.id <- Generatetestid(subtypes.onco)



n.snp <- 205
M <- 5
log.odds.standard.all <-matrix(0,n.snp,M)
log.odds.intrinsic.all <- matrix(0,n.snp,M)
log.odds.intrinsic.dic.all <- matrix(0,n.snp,M)
log.odds.intrinsic.eb.all <- matrix(0,n.snp,M)
log.odds.intrinsic.la.all <- matrix(0,n.snp,M)

for(i1 in 1:n.snp){
  load(paste0("/spin1/users/zhangh24/breast_cancer_data_analysis/risk_prediction/two_stage_model/result/all.model.result",i1,".Rdata"))
  
  log.odds.standard.all[i1,] <- all.model.result[[1]][1:5]
  log.odds.intrinsic.all[i1,] <- diag(all.model.result[[2]])
  log.odds.intrinsic.dic.all[i1,] <- diag(all.model.result[[4]])
  log.odds.intrinsic.eb.all[i1,] <- diag(all.model.result[[5]])
  log.odds.intrinsic.la.all[i1,] <- diag(all.model.result[[6]])
}

subtypes.names <- c("Luminal A","Luminal B",
                    "Luminal B HER2Neg",
                    "HER2 Enriched",
                    "Triple Negative")


auc.summary <- matrix("c",5,5)
auc.com <- matrix(0,5,5)
for(i in 1:5){
  print(i)
  idx.test.case <-  icog.test.id[[1]][[i]]
  idx.test.control <- icog.test.id[[2]][[i]]
  idx.test1 <- c(idx.test.case,idx.test.control)
  y.pheno.mis1.test <- y.pheno.mis1[idx.test1,1]
  x.snp.all.test1 <- x.snp.all1[idx.test1,]
  idx.test.case <-  onco.test.id[[1]][[i]]
  idx.test.control <- onco.test.id[[2]][[i]]
  idx.test2 <- c(idx.test.case,idx.test.control)
  y.pheno.mis2.test <- y.pheno.mis2[idx.test2,1]
  x.snp.all.test2 <- x.snp.all2[idx.test2,]
  y.test <- c(y.pheno.mis1.test,y.pheno.mis2.test)
  x.test <- rbind(as.matrix(x.snp.all.test1),as.matrix(x.snp.all.test2))
  log.odds.standard <- log.odds.standard.all[,i]
  log.odds.intrinsic <- log.odds.intrinsic.all[,i]
  log.odds.intrinsic.dic <-  log.odds.intrinsic.dic.all[,i]
  log.odds.intrinsic.eb <- log.odds.intrinsic.eb.all[,i]
  log.odds.intrinsic.la <- log.odds.intrinsic.la.all[,i]
  auc.cal.result <- GenerateAuc_Cal(
    log.odds.standard,
    log.odds.intrinsic,
    log.odds.dic,
    log.odds.intrinsic.eb,
    log.odds.intrinsic.la,
    x.test,
    y.test
  )    
  auc.result <-   auc.cal.result [[1]]
  auc.95 <- auc.cal.result[[2]]
  auc.com[i,] <- auc.cal.result [[1]]
  auc.summary[i,] <- paste0(auc.result," (",auc.95,")")
  cal.result <- auc.cal.result[[3]]
  sensitivities <-   as.vector(auc.cal.result[[4]])
  specificities <- as.vector(auc.cal.result[[5]])
  n <- length(sensitivities)/5
  method <- c(rep("standard analysis",n),
           rep("intrinsic subtypes",n),
           rep("dichotomized analysis",n),
           rep("Empirical Bayesian (Normal Prior)",n),
           rep("Empirical Bayesian (Laplace Prior)",n))
  n <- 10
  quantile <- rep(c(1:10),5)
  method.cal <- c(rep("Standard analysis",n),
                  rep("Intrinsic subtypes",n),
                  rep("Dichotomized analysis",n),
                  rep("Empirical Bayesian (Normal Prior)",n),
                  rep("Empirical Bayesian (Laplace Prior)",n))
  data.auc <- data.frame(sensitivities,specificities,method)
  png(filename=paste0(subtypes.names[i]," risk prediction.png"),
      width=8,height=6,units="in",res=300)
      print({
        p <-  ggplot(data= data.auc,aes(x=1-specificities,y=sensitivities,col=method))+geom_line()+style_roc()+ggtitle(paste0(subtypes.names[i]," AUC plot"))
                      p       
      }
        
      )
  dev.off()
  data.cal <- data.frame(calresult = as.vector(t(cal.result)),
                         quantile = quantile,
                         method= method.cal)
  png(filename=paste0(subtypes.names[i]," calibration.png"),
      width=8,height=6,units="in",res=300)
  print({
  p <-   ggplot(data= data.cal,aes(x= quantile,y=calresult,col=method))+geom_line()+ggtitle(paste0(subtypes.names[i]," calibration"))+ylab("Odds ratio")+xlab("risk quantile")+
      scale_x_continuous(breaks=c(1:10))
  p
    
  })
  
  dev.off()
}
rownames(auc.summary) <- c("Luminal A","Luminal B",
                           "Luminal B HER2Neg",
                           "HER2 Enriched",
                           "Triple Negative")
colnames(auc.summary) <- c("standard analysis",
                           "intrinsic subtypes",
                           "dichotomized analysis",
                           "Empirical Bayesian (Normal Prior)",
                           "Empirical Bayesian (Laplace Prior)")
write.csv(auc.summary,file="auc.summary.csv")
auc.com 
n <- 5
method <- c(rep("Standard analysis",n),
            rep("Intrinsic subtypes",n),
            rep("Dichotomized analysis",n),
            rep("Empirical Bayesian (Normal Prior)",n),
            rep("Empirical Bayesian (Laplace Prior)",n))
subtypes <- rep(c("Luminal A","Luminal B",
              "Luminal B HER2Neg",
              "HER2 Enriched",
              "Triple Negative"),5)
subtypes <- factor(subtypes,levels=c("Luminal A",
                                           "Luminal B",
                                           "Luminal B HER2Neg",
                                           "HER2 Enriched",
                                           "Triple Negative"))
data <- data.frame(auc = as.vector(auc.com),method=method,subtypes=subtypes)
png(filename=paste0(subtypes.names[i]," auc_summary.png"),
    width=9,height=6,units="in",res=300)
ggplot(data,aes(x= subtypes,y=auc,group=method,col=method))+geom_line()+
  ggtitle(paste0("AUC summary"))
dev.off()


















# +
#   ylab("AUC")
# library(ggplot2)
# library(plotROC)
# 
# 
# library(pROC)
# 
# 
# 


#######standard breast cancer risk prediction

n.snp <- 205
M <- 5
log.odds.standard.all <-rep(0,n.snp,M)
log.odds.intrinsic.all <- matrix(0,n.snp,M)
log.odds.intrinsic.dic.all <- matrix(0,n.snp,M)
log.odds.intrinsic.eb.all <- matrix(0,n.snp,M)
log.odds.intrinsic.la.all <- matrix(0,n.snp,M)

for(i1 in 1:n.snp){
  load(paste0("/spin1/users/zhangh24/breast_cancer_data_analysis/risk_prediction/two_stage_model/result/all.model.result",i1,".Rdata"))
  
  log.odds.standard.all[i1] <- all.model.result[[1]][6]
  log.odds.intrinsic.all[i1,] <- all.model.result[[2]][6,]
  log.odds.intrinsic.dic.all[i1,] <- all.model.result[[4]][6,]
  log.odds.intrinsic.eb.all[i1,] <- all.model.result[[5]][6,]
  log.odds.intrinsic.la.all[i1,] <- all.model.result[[6]][6,]
}

subtypes.names <- c("Luminal A","Luminal B",
                    "Luminal B HER2Neg",
                    "HER2 Enriched",
                    "Triple Negative")
i <- 6
print(i)
idx.test.case <-  icog.test.id[[1]][[i]]
idx.test.control <- icog.test.id[[2]][[i]]
idx.test1 <- c(idx.test.case,idx.test.control)
y.pheno.mis1.test <- y.pheno.mis1[idx.test1,1]
y.pheno.mi1.train <- y.pheno.mis1[-idx.test1,]
x.snp.all.test1 <- x.snp.all1[idx.test1,]
x.snp.all.train1 <- x.snp.all1[-idx.test1,]
idx.test.case <-  onco.test.id[[1]][[i]]
idx.test.control <- onco.test.id[[2]][[i]]
idx.test2 <- c(idx.test.case,idx.test.control)
y.pheno.mis2.test <- y.pheno.mis2[idx.test2,1]
y.pheno.mis2.train <- y.pheno.mis2[-idx.test2,]
x.snp.all.test2 <- x.snp.all2[idx.test2,]
x.snp.all.train2 <- x.snp.all2[-idx.test2,]
y.test <- c(y.pheno.mis1.test,y.pheno.mis2.test)
y.pheno.mis.test <- rbind(y.pheno.mis1[idx.test1,],y.pheno.mis2[idx.test2,])
y.pheno.mis.train <- rbind(y.pheno.mi1.train,y.pheno.mis2.train)
x.snp.train <- rbind(as.matrix(x.snp.all.train1) ,as.matrix(x.snp.all.train2))
x.test <- rbind(as.matrix(x.snp.all.test1),as.matrix(x.snp.all.test2))
#log.odds.standard <- log.odds.standard.all[,i]
prs.standard <- x.test%*%log.odds.standard

cal.standard <- calibration(y.test,prs.standard)
roc.standard <- roc(y.test,as.vector(prs.standard),ci=T,plot=F)


subtypes.prs.intrinsic.train <- x.snp.train%*%log.odds.intrinsic.all
subtypes.train <- as.character(GenerateIntrinsicmis(y.pheno.mis.train[,2],y.pheno.mis.train[,3],
                                                   y.pheno.mis.train[,4],y.pheno.mis.train[,5]))
data.train <- data.frame(Luminal_A_PRS=subtypes.prs.intrinsic.train[,1],
                   Triple_Neg_PRS = subtypes.prs.intrinsic.train[,5],
                   subtypes=subtypes.train)
png(filename=paste0("Luminal_A_vs_TP.png"),
    width=8,height=6,units="in",res=300)
ggplot(data.train,aes(x= Triple_Neg_PRS,y=Luminal_A_PRS,col=subtypes))+geom_point()+ggtitle("PRS Luminal A vs Triple Negative")+
  xlab("Triple Negative PRS")+
  ylab("Luminal A PRS")


cor(subtypes.prs.intrinsic[,1],subtypes.prs.intrinsic[,5])

subtypes.prs.intrinsic <- x.test%*% log.odds.intrinsic.all
try.prs <- apply(subtypes.prs.intrinsic,1,max)
scale.prs <- apply(subtypes.prs.intrinsic,2,function(x){x-mean(x)/sd(x)})
try.prs <- apply(scale.prs,1,max)
try <- roc(y.test,try.prs)
try
try <- roc(y.test~subtypes.prs.intrinsic[,1]+subtypes.prs.intrinsic[,2])

TrueFalseCalculateMult <- function(prs,test.data){
  p <- ncol(prs)
  idx.true <- which(test.data==1)
  idx.false <- which(test.data==0)
  n <- length(test.data)
  n <- 11
  cut.point <- matrix(0,n,p)
  for(i in 1:p){
    min.prs <- range(prs[,i])[1]
    max.prs <- range(prs[,i])[2]
    cut.point[,i] <- seq(from=min.prs,to=max.prs,by=(max.prs-min.prs)/10)
  }
  true.pos <- rep(0,n^p)
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






subtypes.test <- as.character(GenerateIntrinsicmis(y.pheno.mis.test[,2],y.pheno.mis.test[,3],
                                      y.pheno.mis.test[,4],y.pheno.mis.test[,5]))
idx.other <- which(subtypes.test!="Luminal_A"&subtypes.test!="TripleNeg"&subtypes.test!="control")
subtypes.test[idx.other] <- "other subtypes"
data <- data.frame(Luminal_A_PRS=subtypes.prs.intrinsic[,1],
                   Triple_Neg_PRS = subtypes.prs.intrinsic[,5],
                   subtypes=subtypes.test)

png(filename=paste0("Luminal_A_vs_TP.png"),
    width=8,height=6,units="in",res=300)
ggplot(data,aes(x= Triple_Neg_PRS,y=Luminal_A_PRS,col=subtypes))+geom_point()+ggtitle("PRS Luminal A vs Triple Negative")+
  xlab("Triple Negative PRS")+
  ylab("Luminal A PRS")
dev.off()
  

data.clean <- data[!(data$subtypes=="control"|data$subtypes=="other subtypes"),]
cor(data.clean[,1],data.clean[,2])
png(filename=paste0("Luminal_A_vs_TP_clean.png"),
    width=8,height=6,units="in",res=300)

ggplot(data.clean,aes(x= Triple_Neg_PRS,y=Luminal_A_PRS,col=subtypes))+geom_point()+ggtitle("PRS Luminal A vs Triple Negative")+
  xlab("Triple Negative PRS")+
  ylab("Luminal A PRS")
dev.off()







log.odds.intrinsic <- log.odds.intrinsic.all[,i]
log.odds.intrinsic.dic <-  log.odds.intrinsic.dic.all[,i]
log.odds.intrinsic.eb <- log.odds.intrinsic.eb.all[,i]
log.odds.intrinsic.la <- log.odds.intrinsic.la.all[,i]





















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
load("/spin1/users/zhangh24/breast_cancer_data_analysis/risk_prediction/two_stage_model/result/log.odds.meta.la.all.Rdata")
##triple vs nontriple
load("/spin1/users/zhangh24/breast_cancer_data_analysis/risk_prediction/two_stage_model/result/log.odds.meta.tvn.triple.Rdata")
load("/spin1/users/zhangh24/breast_cancer_data_analysis/risk_prediction/two_stage_model/result/log.odds.meta.two.stage.tvn.all.Rdata")
load("/spin1/users/zhangh24/breast_cancer_data_analysis/risk_prediction/two_stage_model/result/sigma.log.odds.two.stage.tvn.Rdata")
load("/spin1/users/zhangh24/breast_cancer_data_analysis/risk_prediction/two_stage_model/result/heter.sigma.tvn.Rdata")

####additive two-stage model
load("/spin1/users/zhangh24/breast_cancer_data_analysis/risk_prediction/two_stage_model/result/p.heter.add.Rdata")
load("/spin1/users/zhangh24/breast_cancer_data_analysis/risk_prediction/two_stage_model/result/log.odds.add.triple.Rdata")
load("/spin1/users/zhangh24/breast_cancer_data_analysis/risk_prediction/two_stage_model/result/log.odds.add.triple.eb.Rdata")






##############standard model


library(pROC)
library(ggplot2)
prs.standard <- x.snp.all.test%*%log.odds.meta
cal.standard <- calibration(y.test[,1],prs.standard)
plot(cal.standard,type="l",xlab="quantile",ylab="odds ratio",xlim=c(1,10),ylim=c(0,3),main="calibration")
roc.standard <- roc(y.test[,1],as.vector(prs.standard),ci=T,plot=T)
auc.standard <- roc.standard$auc
roc.result.standard <- true.false.calculate(prs.standard,y.test[,1])
n <- nrow(roc.result.standard)
lab = rep("standard",n)
# plot(roc.result.standard[,1],roc.result.standard[,2],xlab="false_p",ylab="true_p")
# abline(a=0,b=1,col="red")
 

# auc.standard <- auc_cal(roc.result.standard)
roc.result.standard.mer <- data.frame(roc.result.standard,lab)
#############intrinsic subtypes
prs.intrinsic <- x.snp.all.test%*%log.odds.meta.triple
cal.intrinsic <- calibration(y.test[,1],prs.intrinsic)
plot(cal.intrinsic,type="l",xlab="quantile",ylab="odds ratio",xlim=c(1,10),ylim=c(0,3),main="calibration")
roc.intrinsic <- roc(y.test[,1],as.vector(prs.intrinsic),plot=T,ci=T)
auc.intrinsic <- roc.intrinsic$auc
roc.result.intrinsic <- true.false.calculate(prs.intrinsic,y.test[,1])
# plot(roc.result.intrinsic[,1],roc.result.intrinsic[,2],xlab="false_p",ylab="true_p")
# abline(a=0,b=1,col="red")
#lab = rep("intrinsic subtypes",n)
#auc.intrinsic <- auc_cal(roc.result.intrinsic)
roc.reuslt.intrinsic.mer <- data.frame(roc.result.intrinsic,lab)
###########intrinsic subtypes dic mixed
log.odds.intrinsic.dic <- log.odds.meta
idx <- which(p.heter.intrinsic<0.05)
log.odds.intrinsic.dic[idx] <- log.odds.meta.triple[idx]

prs.intrinsic.dic <- x.snp.all.test%*%log.odds.intrinsic.dic
cal.intrinsic.dic<- calibration(y.test[,1],prs.intrinsic.dic)
plot(cal.intrinsic.dic,type="l",xlab="quantile",ylab="odds ratio",xlim=c(1,10),ylim=c(0,3),main="calibration")
roc.result.intrinsic.dic <- true.false.calculate(prs.intrinsic.dic,y.test[,1])
roc.intrinsic.dic <- roc(y.test[,1],as.vector(prs.intrinsic.dic),ci=T,plot=T)
auc.intrinsic.dic <- roc.intrinsic.dic$auc
# plot(roc.result.intrinsic.dic[,1],roc.result.intrinsic.dic[,2],xlab="false_p",ylab="true_p")
# abline(a=0,b=1,col="red")

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
roc.intrinsic.eb <- roc(y.test[,1],as.vector(prs.intrinsic.eb),ci=T,plot=T)
plot(cal.intrinsic.eb,type="l",xlab="quantile",ylab="odds ratio",xlim=c(1,10),ylim=c(0,3),main="calibration")
auc.intrinsic.eb <- roc.intrinsic.eb$auc
#plot(roc.result.intrinsic.eb[,1],roc.result[,2],xlab="false_p",ylab="true_p")
#abline(a=0,b=1,col="red")

lab = rep("intrinsic eb",n)
roc.result.intrinsic.eb <- data.frame(roc.result.intrinsic.eb,lab)

###########intrinsic laplace
log.odds.intrinsic.la <- log.odds.meta.la.all[,5]

prs.intrinsic.la <- x.snp.all.test%*%log.odds.intrinsic.la
cal.intrinsic.la <- calibration(y.test[,1],prs.intrinsic.la)
plot(cal.intrinsic.la,type="l",xlab="quantile",ylab="odds ratio",xlim=c(1,10),ylim=c(0,3),main="calibration")
roc.result.intrinsic.la <- true.false.calculate(prs.intrinsic.la ,y.test[,1])
roc.intrinsic.la <- roc(y.test[,1],as.vector(prs.intrinsic.la),ci=T,plot=T)
auc.intrinsic.la <- roc.intrinsic.la$auc
#plot(roc.result.intrinsic.la[,1],roc.result[,2],xlab="false_p",ylab="true_p")
#abline(a=0,b=1,col="red")
#auc_cal(roc.result.intrinsic.la)
lab = rep("intrinsic la",n)
roc.result.intrinsic.la <- data.frame(roc.result.intrinsic.la,lab)



##triple vs nontriple 

prs.tvn <- x.snp.all.test%*%log.odds.meta.tvn.triple
cal.tvn <- calibration(y.test[,1],prs.tvn)
roc.result.tvn <- true.false.calculate(prs.tvn,y.test[,1])
plot(roc.result.tvn[,1],roc.result.tvn[,2],xlab="false_p",ylab="true_p")
abline(a=0,b=1,col="red")
lab = rep("tvn subtypes",n)
auc.tvn <- auc_cal(roc.result.tvn)
roc.reuslt.tvn.mer <- data.frame(roc.result.tvn,lab)

#########triple vs n triple eb

log.odds.tvn.eb <- rep(0,205)
for(i in 1:205){
  
  logodds.subtype <- log.odds.meta.two.stage.all[i,]
  M <- length(logodds.subtype)
  sigma.subtype <- matrix(sigma.log.odds.two.stage[i,],5,5)
  log.odds.tvn.eb[i] <- ebestimate(logodds.subtype,sigma.subtype,
                                         as.numeric(log.odds.meta[i]),
                                         as.numeric(heter.sigma.tvn[i]))[5]
}

prs.tvn.eb <- x.snp.all.test%*%log.odds.tvn.eb
cal.tvn.eb <- calibration(y.test[,1],prs.tvn.eb)
roc.result.tvn.eb <- true.false.calculate(prs.tvn.eb ,y.test[,1])
plot(roc.result.tvn.eb[,1],roc.result[,2],xlab="false_p",ylab="true_p")
abline(a=0,b=1,col="red")
auc.tvn.eb <- auc_cal(roc.result.tvn.eb)
lab = rep("tvn eb",n)
roc.result.tvn.eb <- data.frame(roc.result.tvn.eb,lab)





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





