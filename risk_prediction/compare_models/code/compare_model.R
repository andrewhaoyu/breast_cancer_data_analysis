setwd('/data/zhangh24/breast_cancer_data_analysis/risk_prediction/two_stage_model/result/')
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
library(bcutility,lib.loc ="/home/zhangh24/R/x86_64-pc-linux-gnu-library/3.6/")
library(bc2,lib.loc ="/home/zhangh24/R/x86_64-pc-linux-gnu-library/3.6/")
library(pROC)
library(data.table)
library(plotROC)
icog.data <- as.data.frame(fread("/data/zhangh24/breast_cancer_data_analysis/data/sig_snp_icog_prs.csv",header=T))
icog.data <- icog.data[,-1]
onco.data <- as.data.frame(fread("/data/zhangh24/breast_cancer_data_analysis/data/sig_snp_onco_prs.csv",header=T))
onco.data <- onco.data[,-1]
library(tidyverse)
y.pheno.mis1 <- select(icog.data,Behavior,ER,PR,HER2,Grade)
############put the people with unknown case status as 1
idx.insi.unknown.icog <- which(y.pheno.mis1[,1]==888|
                                 y.pheno.mis1[,1]==2) 
y.pheno.mis1[idx.insi.unknown.icog,1] = 1
subtypes.icog <- GenerateIntrinsicmis(y.pheno.mis1[,2],y.pheno.mis1[,3],
                                      y.pheno.mis1[,4],y.pheno.mis1[,5])

#table(subtypes.icog)+table(subtypes.onco)

x.covar1 <- select(icog.data,7:16)
x.snp.all1 <- select(icog.data,17:228)
colnames(y.pheno.mis1)
# split.id <- list(id1.train,
#                  id2.train,
#                  id2.test,
#                  id.cohort.clean1,
#                  id.cohort.clean2)
#load the training and testing data row number
genetic_correlation <- as.matrix(read.csv("/data/zhangh24/breast_cancer_data_analysis/risk_prediction/two_stage_model/result/genetic_correlation.csv",header=T))


setwd('/data/zhangh24/breast_cancer_data_analysis/')
load(paste0("./risk_prediction/result/split.id.rdata"))
#icog.test.id <- Generatetestid(subtypes.icog)
icog.train.id <- split.id[[1]]
onco.train.id <- split.id[[2]]
onco.test.id <- split.id[[3]]
icog.cohort.id <- split.id[[4]]
onco.cohort.id <- split.id[[5]]

y.pheno.mis2 <- select(onco.data,Behavior,ER,PR,HER2,Grade)
############put onco array unknown cases as 1
idx.insi.unknown.onco <- which(y.pheno.mis2[,1]==888|
                                 y.pheno.mis2[,1]==2) 
y.pheno.mis2[idx.insi.unknown.onco,1] <- 1
subtypes.onco <- GenerateIntrinsicmis(y.pheno.mis2[,2],
                                      y.pheno.mis2[,3],
                                      y.pheno.mis2[,4],
                                      y.pheno.mis2[,5])
x.covar2 <- select(onco.data,7:16)
x.snp.all2 <- select(onco.data,17:228)
colnames(y.pheno.mis2)
subtypes.onco <- GenerateIntrinsicmis(y.pheno.mis2[,2],y.pheno.mis2[,3],
                                      y.pheno.mis2[,4],y.pheno.mis2[,5])
n.snp <- ncol(x.snp.all1)

M <- 5
log.odds.standard.all <-matrix(0,n.snp,M)
log.odds.poly.all <-matrix(0,n.snp,M)
log.odds.intrinsic.all <- matrix(0,n.snp,M)
log.odds.intrinsic.eb.all <- matrix(0,n.snp,M)
log.odds.intrinsic.ge.all <- matrix(0,n.snp,M)
log.odds.intrinsic.ep.all <- matrix(0,n.snp,M)
# log.odds.intrinsic.la.all <- matrix(0,n.snp,M)
# log.odds.intrinsic.tree.all <- matrix(0,n.snp,M)

for(i1 in 1:n.snp){
  load(paste0("/data/zhangh24/breast_cancer_data_analysis/risk_prediction/two_stage_model/result/all.model.result",i1,".Rdata"))
  
  log.odds.standard.all[i1,] <- rep(all.model.result[[1]],M)
  log.odds.poly.all[i1,] <- all.model.result[[4]]
  log.odds.intrinsic.all[i1,] <- all.model.result[[2]]
  log.odds.intrinsic.eb.all[i1,] <- all.model.result[[5]]
  log.odds.intrinsic.ge.all[i1,] <- all.model.result[[6]]
  # log.odds.intrinsic.la.all[i1,] <- diag(all.model.result[[6]])
  # log.odds.intrinsic.tree.all[i1,] <- diag(all.model.result[[8]])
}
#prior.sigma <- cov(log.odds.intrinsic.all)


for(i1 in 1:n.snp){
  load(paste0("/data/zhangh24/breast_cancer_data_analysis/risk_prediction/two_stage_model/result/all.model.result",i1,".Rdata"))
  log.odds.intrinsic <- all.model.result[[2]]
  sigma.log.odds.intrinsic <- matrix(all.model.result[[3]],M,M)
  #######empirical bayesian with whole genome genetic correlaiton
  heter.variance.estiamte.R <- HeterVarianceEstimateNew(
    log.odds.intrinsic,
    sigma.log.odds.intrinsic,
    genetic_correlation
  )
  prior.sigma <-    heter.variance.estiamte.R*as.matrix(genetic_correlation)
  log.odds.intrinsic.ge.all[i1,] <- EbestimateNew(
    log.odds.intrinsic,
    sigma.log.odds.intrinsic,
    all.model.result[[1]],
    prior.sigma
  )
  #######empirical bayesian with empirical genetic correlation
  heter.variance.estiamte.R <- HeterVarianceEstimateNew(
    log.odds.intrinsic,
    sigma.log.odds.intrinsic,
    cor(log.odds.intrinsic.all)
  )
  prior.sigma <-    heter.variance.estiamte.R*as.matrix(genetic_correlation)
  
  log.odds.intrinsic.ep.all[i1,] <- EbestimateNew(
    log.odds.intrinsic,
    sigma.log.odds.intrinsic,
    all.model.result[[1]],
    prior.sigma
  )
  # log.odds.standard.all[i1,] <- rep(all.model.result[[1]],M)
  # log.odds.poly.all <- all.model.result[[4]]
  # log.odds.intrinsic.all[i1,] <- all.model.result[[2]]
  # log.odds.intrinsic.eb.all[i1,] <- all.model.result[[5]]
  # log.odds.intrinsic.ge.all[i1,] <- all.model.result[[6]]
  # log.odds.intrinsic.la.all[i1,] <- diag(all.model.result[[6]])
  # log.odds.intrinsic.tree.all[i1,] <- diag(all.model.result[[8]])
}

#save(temp,file="/data/zhangh24/breast_cancer_data_analysis/risk_prediction/two_stage_model/result/temp.Rdata")
subtypes.names <- c("Luminal_A","Luminal_B",
                    "Luminal_B_HER2Neg",
                    "HER2Enriched",
                    "TripleNeg")
plot(log.odds.standard.all[,5],log.odds.intrinsic.all[,5])
plot(log.odds.intrinsic.all[,5],log.odds.intrinsic.eb.all[,5])
plot(log.odds.intrinsic.all[,5],log.odds.intrinsic.ge.all[,5])

#########
n.method <- 6
n.subtype <- 5
auc.summary <- matrix("c",n.subtype,n.method)
auc.com <- matrix(0,n.subtype,n.method)
coeff <- matrix(0,n.subtype,n.method)
coeff_95 <- matrix("c",n.subtype,n.method)

##########split the data
onco.test <- which(onco.data[,1]%in%onco.test.id)
y.pheno.mis2.test <- y.pheno.mis2[onco.test,]
x.covar.test <- x.covar2[onco.test,]
x.snp.all.test <- x.snp.all2[onco.test,]
#sample.id.order <- onco.data[onco.test,1]
#sample.id.order.j <- sample.id.order[onco.test.sub.j]
insub.onco.test <- GenerateIntrinsicmis(y.pheno.mis2.test[,2],
                                        y.pheno.mis2.test[,3],
                                        y.pheno.mis2.test[,4],
                                        y.pheno.mis2.test[,5])



for(j in 1:M){

  onco.test.sub.j <- which(insub.onco.test==subtypes.names[j]|
                             insub.onco.test=="control")
  y.pheno.mis2.test.casecon.j <- as.vector(y.pheno.mis2.test[onco.test.sub.j,1])
  x.snp.j <- as.matrix(x.snp.all.test[onco.test.sub.j,])
  #sample.id.order.j <- sample.id.order[onco.test.sub.j]
  #prs.intrinsic = x.snp.j%*%log.odds.intrinsic.all[,j]
  #true_result <- data.frame(sample.id.order.j,y.pheno.mis2.test.casecon.j,prs.intrinsic,stringsAsFactors=F)
  #colnames(true_result) <- c("IID","case_status_t","prs_true")
  # for(i in 1:nrow(true_result)){
  # if(as.numeric(true_result[i,1])==100000){
  #   true_result[i,1]= paste0("sample_1e+05")
  # }else{
  #   true_result[i,1]=paste0("sample_",true_result[i,1])
  # }
  # }
  # idx.match <- match(test.sample_new[,1],true_result[,1])
  # x.snp.j.order <- x.snp.j[idx.match,]
  

  
  auc.cal.result <- GenerateAuc_Cal(
    log.odds.standard=
      as.vector(log.odds.standard.all[,j]),
    log.odds.intrinsic=
      log.odds.poly.all[,j],
    log.odds.dic=
      log.odds.intrinsic.all[,j],
    log.odds.intrinsic.eb=
      log.odds.intrinsic.eb.all[,j],
    log.odds.intrinsic.la=
      log.odds.intrinsic.ge.all[,j],
    log.odds.intrinsic.tree=
      log.odds.intrinsic.ep.all[,j],
    x.test=
      x.snp.j,
    y.test=
      y.pheno.mis2.test.casecon.j
  )   
  #snp =  names(select(onco.data,17:228))
  #logodds_result = data.frame(snp,log.odds.intrinsic.all,stringsAsFactors=F)
  #colnames(logodds_result) <- c("snp",subtypes.names)
  #write.csv(logodds_result,file = "/data/zhangh24/breast_cancer_data_analysis/risk_prediction/two_stage_model/result/temp_logodds_result.csv")
  #colnames(temp) <- c("y_test","score")
  #temp[,1] = as.numeric(temp[,1])
  #temp = as.data.frame(temp)
  #temp_new = temp %>% group_by(y_test) %>% mutate(means=mean(score))
  coeff.result <- GenerateCoeff(
    log.odds.standard.all[,j],
    log.odds.poly.all[,j],
    log.odds.intrinsic.all[,j],
    log.odds.intrinsic.eb.all[,j],
    log.odds.intrinsic.ge.all[,j],
    log.odds.intrinsic.ep.all[,j],
    x.snp.j,
    y.pheno.mis2.test.casecon.j
  )
  coeff[j,] <- coeff.result[[1]]
  coeff_95[j,] <- coeff.result[[2]]
  auc.result <-   auc.cal.result [[1]]
  auc.95 <- auc.cal.result[[2]]
  auc.com[j,] <- auc.cal.result [[1]]
  auc.summary[j,] <- paste0(auc.result," (",auc.95,")")
  cal.result <- auc.cal.result[[3]]
  sensitivities <-   as.vector(auc.cal.result[[4]])
  specificities <- as.vector(auc.cal.result[[5]])
  
}









rownames(auc.summary) <- c("Luminal A","Luminal B",
                           "Luminal B HER2Neg",
                           "HER2 Enriched",
                           "Triple Negative")
colnames(auc.summary) <- c("standard analysis",
                           "poly subtypes",
                           "two-stage analysis",
                           "Empirical Bayesian (Normal Prior)",
                           "Empirical Bayesian (Genetic correlation)",
                           "Empirical Bayesian (Empirical prior)")

# colnames(auc.summary) <- c("standard analysis",
#                            "intrinsic subtypes",
#                            "dichotomized analysis",
#                            "Empirical Bayesian (Normal Prior)",
#                            "Empirical Bayesian (Laplace Prior)",
#                            "Classification tree")
write.csv(auc.summary,file="/data/zhangh24/breast_cancer_data_analysis/risk_prediction/two_stage_model/result/auc.summary.csv")
auc.com 
n <- 5
method <- c(rep("Standard analysis",n),
            rep("Intrinsic subtypes",n),
            rep("Dichotomized analysis",n),
            rep("Empirical Bayesian (Normal Prior)",n),
            rep("Empirical Bayesian (Laplace Prior)",n),
            rep("Classification tree"),n)
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


















library(ggplot2)
library(plotROC)


library(pROC)





#######standard breast cancer risk prediction

n.snp <- 205
M <- 5
log.odds.standard.all <-rep(0,n.snp,M)
log.odds.intrinsic.all <- matrix(0,n.snp,M)
log.odds.intrinsic.dic.all <- matrix(0,n.snp,M)
log.odds.intrinsic.eb.all <- matrix(0,n.snp,M)
log.odds.intrinsic.la.all <- matrix(0,n.snp,M)

for(i1 in 1:n.snp){
  load(paste0("/data/zhangh24/breast_cancer_data_analysis/risk_prediction/two_stage_model/result/all.model.result",i1,".Rdata"))

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
prs.standard <- x.test%*%log.odds.standard.all
library(pROC)
cal.standard <- calibration(y.test,prs.standard)
roc.standard <- roc(y.test,as.vector(prs.standard),ci=T,plot=F)


#######Generate the PRS using empirical Bayesian with normal proior

subtypes.prs.intrinsic.eb.test <- x.test%*%log.odds.intrinsic.eb.all
subtypes.test <- as.character(GenerateIntrinsicmis(y.pheno.mis.test[,2],y.pheno.mis.test[,3],
                                                   y.pheno.mis.test[,4],y.pheno.mis.test[,5]))
data.test <- data.frame(Luminal_A_PRS=subtypes.prs.intrinsic.eb.test[,1],
                        Luminal_B_PRS = subtypes.prs.intrinsic.eb.test[,2],
                        Luminal_B_HER2_en_PRS = subtypes.prs.intrinsic.eb.test[,3],
                        HER2_en_PRS = subtypes.prs.intrinsic.eb.test[,4],
                        
                         
                         Triple_Neg_PRS = subtypes.prs.intrinsic.eb.test[,5],
                         subtypes=subtypes.test)

all.subtypes <- c("Luminal_A",
                  "Luminal_B",
                  "Luminal_B_HER2Neg",
                  "HER2Enriched",
                  "TripleNeg")

coefficients <- rep(0,5)
control.mean <- rep(0,5)
control.sd <- rep(0,5)
for(i in 1:5){
  control.mean[i] <- mean(data.test[,i])
  control.sd[i] <- sd(data.test[,i])
  subtypes.names.temp <- all.subtypes[i]
  subtype <- all.subtypes[i]
  idx <- which(subtypes.test==subtype|subtypes.test=="control")
  subtype.y <- factor(subtypes.test[idx],levels=c("control",subtype))
  prs.x <- data.test[idx,i]
  model <- glm(subtype.y~prs.x,family=binomial(link='logit'))
  coefficients[i] <- coef(model)[2]
  
}

control.prs.infor <- data.frame(subtypes=all.subtypes,
                                coefficients=coefficients,
                                control.mean = control.mean,
                                control.sd = control.sd
)
write.csv(control.prs.infor,file="control.prs.infor.csv")




data.test.clean <- data.test[data.test$subtypes=="control",]




#############generate the control PRS using empirical bayesian on training dataset
subtypes.prs.intrinsic.train <- x.snp.train%*%log.odds.intrinsic.eb.all
standard.prs <- as.numeric(x.snp.train%*%log.odds.standard.all)
subtypes.train <- as.character(GenerateIntrinsicmis(y.pheno.mis.train[,2],y.pheno.mis.train[,3],
                                                   y.pheno.mis.train[,4],y.pheno.mis.train[,5]))
data.train <- data.frame(Luminal_A_PRS=subtypes.prs.intrinsic.train[,1],
                         
                   Triple_Neg_PRS = subtypes.prs.intrinsic.train[,5],
                   standard_prs = standard.prs,
                   subtypes=subtypes.train)
png(filename=paste0("Luminal_A_vs_TP.png"),
    width=8,height=6,units="in",res=300)
ggplot(data.train,aes(x= Triple_Neg_PRS,y=Luminal_A_PRS,col=subtypes))+geom_point()+ggtitle("PRS Luminal A vs Triple Negative")+
  xlab("Triple Negative PRS")+
  ylab("Luminal A PRS")
dev.off()

data.train.clean <- data.train[data.train$subtypes=="control",]

cor(data.train.clean[,1],data.train.clean[,2])
cor(data.train.clean[,1],data.train.clean[,3])
cor(data.train.clean[,2],data.train.clean[,3])
prs.LA <- data.train.clean[,1]
prs.TN <- data.train.clean[,2]
prs.sd <- data.train.clean[,3]






library(dplyr)
##########plot for LA






# 
# try <-  data.train%>%
#   group_by(group1,group2) %>%
#   count(HD)

plot.subtype.name <- c("Luminal A","Triple Negative")
for(i in 1:2){
  names.col <- colnames(data.train)
  emprical.pro <- crosstaub(data.train.clean[,i],prs.sd)
  high.emprical.pro <- emprical.pro[9:11,]/sum(emprical.pro[9:11,])
  high.emprical.pro.m <- melt(high.emprical.pro)
  conditional.pro <- apply(emprical.pro,2,function(x){x/sum(x)})
  conditional.m <- melt(conditional.pro)
  
  emprical.pro.new <- crosstaub_new(data.train.clean[,i],prs.sd)
  emprical.pro.new.m <- melt(emprical.pro.new)
  upper <- 0.04
  emprical.pro.new.m$value[emprical.pro.new.m$value>=upper] <- upper
  #write.csv(emprical.pro,file=paste0(plot.subtype.name[i],"vs standard PRS probability.csv"))
  
  
  
  png(filename=paste0(plot.subtype.name[i],"_vs_standard.png"),
      width=8,height=6,units="in",res=300)
  print(
    {p <- ggplot(data.train.clean,aes(x= standard_prs,y=data.train.clean[,i]))+geom_point()+ggtitle(paste0("PRS" ,plot.subtype.name[i], " vs Standard analysis"))+
      xlab("Standard PRS")+
      ylab(paste0(plot.subtype.name[i]," PRS"))
    p
    }
  )
  dev.off()
  
  png(filename=paste0(plot.subtype.name[i],"_vs_standard_heatmap_absolute.png"),
      width=8,height=6,units="in",res=300)
  print({
    ggplot(emprical.pro.new.m, aes(group2, group1)) + 
      geom_tile(aes(fill = value),colour = "white") + 
      scale_fill_gradient(limits=c(0,0.04),low = "white",  high = "dodgerblue4")+
      xlab("Standard PRS percentile")+
      ylab( paste0(plot.subtype.name[i]," PRS percentile"))+
      ggtitle( paste0(plot.subtype.name[i]," vs Standard Analysis"))
    p
  })
  dev.off()
  
  
  
  
  png(filename=paste0(plot.subtype.name[i],"_vs_standard_heatmap.png"),
      width=8,height=6,units="in",res=300)
  print({
  p <-  ggplot(conditional.m, aes(group2, group1)) + 
      geom_tile(aes(fill = value),colour = "white") + 
      scale_fill_gradient(limits=c(0,0.7),low = "white",  high = "dodgerblue4")+
      xlab("Standard PRS percentile")+
      ylab( paste0(plot.subtype.name[i]," PRS percentile"))+
      ggtitle( paste0(plot.subtype.name[i]," vs Standard"))
    dev.off()
    p
  })
  
  
  png(filename=paste0(plot.subtype.name[i],"_vs_standard_heatmap_high.png"),
      width=8,height=6,units="in",res=300)
  print({
    p <- ggplot(high.emprical.pro.m, aes(group2, group1)) + 
      geom_tile(aes(fill = value),colour = "white") + 
      scale_fill_gradient(low = "white",  high = "dodgerblue4")+
      xlab("Standard PRS percentile")+
      ylab( paste0(plot.subtype.name[i]," PRS percentile"))+
      ggtitle( paste0(plot.subtype.name[i]," vs Standard"))
    dev.off()
  p  
  })
  
  
  
}


# write.csv(result,file="./LuAvsTN.csv")
# heatmap.2(result,
#           tracecol=NA)
library(plyr)
library(scales)


  
png(filename=paste0("Luminal_A_vs_TP_heatmap_absolute.png"),
    width=8,height=6,units="in",res=300)
ggplot(emprical.m,aes(group1,group2))+
  geom_tile(aes(fill = value),colour = "white") + 
  scale_fill_gradient(low = "white",  high = "steelblue")+
  xlab("Luminal A PRS percentile")+
  ylab("Triple Negative PRS percentile")+
  ggtitle("Luminal A vs Triple Negative")
dev.off()
# ,cexRow=1,cexCol=1,margins=c(10,12),col = col,breaks=pal.breaks,key.ylab="",key.title = "",
#           main=" Genetic Correlation Heatmap",dendrogram="row",density.info="none",lwid = c(1.5,4))


qun1 <- quantile(data.train.clean[,1],probs=
                   c(0.01,0.05,0.10,0.20,
                     0.40,0.60,0.80,
                     0.9,
                     0.95,0.99))








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
  

data.clean <- data[data$subtypes=="control",]
cor(data.clean[,1],data.clean[,2])
png(filename=paste0("Luminal_A_vs_TP_clean.png"),
    width=8,height=6,units="in",res=300)

ggplot(data.clean,aes(x= Triple_Neg_PRS,y=Luminal_A_PRS))+geom_point(col="red")+ggtitle("PRS Luminal A vs Triple Negative")+
  xlab("Triple Negative PRS")+
  ylab("Luminal A PRS")
dev.off()







log.odds.intrinsic <- log.odds.intrinsic.all[,i]
log.odds.intrinsic.dic <-  log.odds.intrinsic.dic.all[,i]
log.odds.intrinsic.eb <- log.odds.intrinsic.eb.all[,i]
log.odds.intrinsic.la <- log.odds.intrinsic.la.all[,i]





















library(data.table)
icog.data <- as.data.frame(fread("/data/zhangh24/breast_cancer_data_analysis/data/sig_snps_icog.csv",header=T))
onco.data <- as.data.frame(fread("/data/zhangh24/breast_cancer_data_analysis/data/sig_snps_onco.csv",header=T))
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
load("/data/zhangh24/breast_cancer_data_analysis/risk_prediction/standard_analysis/result/log.odds.meta.Rdata")
#####intrinsic subtype result
load("/data/zhangh24/breast_cancer_data_analysis/risk_prediction/two_stage_model/result/log.odds.meta.triple.Rdata")
load("/data/zhangh24/breast_cancer_data_analysis/risk_prediction/two_stage_model/result/p.heter.intrinsic.Rdata")
load("/data/zhangh24/breast_cancer_data_analysis/risk_prediction/two_stage_model/result/log.odds.meta.two.stage.all.Rdata")
load("/data/zhangh24/breast_cancer_data_analysis/risk_prediction/two_stage_model/result/sigma.log.odds.two.stage.Rdata")
load("/data/zhangh24/breast_cancer_data_analysis/risk_prediction/two_stage_model/result/log.odds.meta.triple.Rdata")
load("/data/zhangh24/breast_cancer_data_analysis/risk_prediction/two_stage_model/result/heter.sigma.Rdata")
load("/data/zhangh24/breast_cancer_data_analysis/risk_prediction/two_stage_model/result/log.odds.meta.la.all.Rdata")
##triple vs nontriple
load("/data/zhangh24/breast_cancer_data_analysis/risk_prediction/two_stage_model/result/log.odds.meta.tvn.triple.Rdata")
load("/data/zhangh24/breast_cancer_data_analysis/risk_prediction/two_stage_model/result/log.odds.meta.two.stage.tvn.all.Rdata")
load("/data/zhangh24/breast_cancer_data_analysis/risk_prediction/two_stage_model/result/sigma.log.odds.two.stage.tvn.Rdata")
load("/data/zhangh24/breast_cancer_data_analysis/risk_prediction/two_stage_model/result/heter.sigma.tvn.Rdata")

####additive two-stage model
load("/data/zhangh24/breast_cancer_data_analysis/risk_prediction/two_stage_model/result/p.heter.add.Rdata")
load("/data/zhangh24/breast_cancer_data_analysis/risk_prediction/two_stage_model/result/log.odds.add.triple.Rdata")
load("/data/zhangh24/breast_cancer_data_analysis/risk_prediction/two_stage_model/result/log.odds.add.triple.eb.Rdata")






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
write.csv(cal.result,file="/data/zhangh24/breast_cancer_data_analysis/risk_prediction/two_stage_model/result/cal.result.csv")










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





