#install_github("andrewhaoyu/bc2",ref='development', args = c('--library="/home/zhangh24/R/x86_64-pc-linux-gnu-library/3.4"'))
###1 represent Icog
###2 represent Onco

rm(list=ls())
#commandarg <- commandArgs(trailingOnly=F)
#myarg <- commandarg[length(commandarg)]
#myarg <- sub("-","",myarg)
#i1 <- as.numeric(myarg)
args = commandArgs(trailingOnly = T)
i1 = as.numeric(args[[1]])

setwd("/spin1/users/zhangh24/breast_cancer_data_analysis/")
library(readr)
library(devtools)
library(CompQuadForm)
library(bc2)
library(data.table)
z.design <- matrix(c(
  c(0,1,1,1,0,0,0,0),
  c(0,0,0,0,0,1,1,1),
  c(0,0,0,0,1,0,0,0),
  c(1,0,0,0,0,0,0,0)
),ncol=4)

colnames(z.design) <- c("Luminial A",
                        "Luminal B ",
                        "HER2 Enriched ",
                        "Triple Negative")

generate_pairwise_differ <- function(n){
  temp <- combn(n,2)
  n.col <- ncol(temp)
  n.row <- n
  result <- matrix(0,n.row,n.col)
  for(i in 1:n.col){
    result[temp[1,i],i] <- 1
    result[temp[2,i],i] <- -1
  }
  return(t(result))
}

distance.fun <- function(logodds,sigma,dis.trans){

n.dis <- nrow(dis.trans)
result <- rep(0,n.dis)
for(i in 1:n.dis){
  dis.temp <- dis.trans[i,,drop=F]
  result[i] <- dis.temp%*%logodds%*%
    solve(dis.temp%*%sigma%*%t(dis.temp))%*%
    t(dis.temp%*%logodds)
  
}
  return(result)
  
}


dis.trans <- generate_pairwise_differ(4)

if(i1<=177){
  ##analysis for Icog
  data1 <- as.data.frame(fread("./data/iCOGS_euro_v10_10232017.csv",header=T))
  y.pheno.mis1 <- cbind(data1$Behaviour1,data1$PR_status1,data1$ER_status1,data1$HER2_status1)
  colnames(y.pheno.mis1) = c("Behavior","PR","ER","HER2")
  # Grade1.fake <- data1$Grade1
  # Grade1.fake[data1$Grade1==2|data1$Grade1==3] <- 1
  # Grade1.fake[data1$Grade1==1] <- 0
  #y.pheno.mis1 <- cbind(data1$Behaviour1,data1$PR_status1,data1$ER_status1,data1$HER2_status1,Grade1.fake)
  # y.pheno.mis1 <- cbind(data1$Behaviour1,data1$PR_status1,data1$ER_status1,data1$HER2_status1)
  
  x.test.all.mis1 <- data1[,c(27:204)]
  
  x.covar.mis1 <- data1[,5:14]
  x.all.mis1 <- as.matrix(cbind(x.test.all.mis1[,i1],x.covar.mis1))
  colnames(x.all.mis1)[1] <- "gene"
  
  Heter.result.Icog = EMmvpolySelfDesign(y.pheno.mis1,x.self.design = x.all.mis1[,1,drop=F],z.design=z.design,baselineonly = NULL,additive = x.all.mis1[,2:11],pairwise.interaction = NULL,saturated = NULL,missingTumorIndicator = 888)
  z.standard <- Heter.result.Icog[[12]]
  z.additive.design <- as.matrix(cbind(1,z.standard))
  M <- nrow(z.standard)
  number.of.tumor <- ncol(z.standard)
  log.odds.icog <- Heter.result.Icog[[1]][(M+1):(M+1+number.of.tumor)]
  
  sigma.log.odds.icog <- Heter.result.Icog[[2]][(M+1):(M+1+number.of.tumor),(M+1):(M+1+number.of.tumor)]

  
  
  #analysis for Onco Array
  #data2 <- read.csv("./V10/Onco_euro_v10_05242017.csv",header=T)
  data2 <- as.data.frame(fread("./data/Onco_euro_v10_10232017.csv",header=T))
  
  
  y.pheno.mis2 <- cbind(data2$Behaviour1,data2$PR_status1,data2$ER_status1,data2$HER2_status1)
  #y.pheno.mis2 <- cbind(data2$Behaviour1,data2$PR_status1,data2$ER_status1,data2$HER2_status1)
  colnames(y.pheno.mis2) = c("Behaviour","PR",
                             "ER","HER2")
  
  x.test.all.mis2 <- data2[,c(27:203)]
  x.covar.mis2 <- data2[,5:14]
  x.all.mis2 <- as.matrix(cbind(x.test.all.mis2[,i1],x.covar.mis2))
  colnames(x.all.mis2)[1] = "gene"
  
  
  Heter.result.Onco = EMmvpolySelfDesign(y.pheno.mis2,x.self.design = x.all.mis2[,1,drop=F],z.design = z.design,baselineonly = NULL,additive = x.all.mis2[,2:11],pairwise.interaction = NULL,saturated = NULL,missingTumorIndicator = 888)
  z.standard <- Heter.result.Onco[[12]]
  M <- nrow(z.standard)
  number.of.tumor <- ncol(z.standard)
  log.odds.onco <- Heter.result.Onco[[1]][(M+1):(M+1+number.of.tumor)]
  sigma.log.odds.onco <- Heter.result.Onco[[2]][(M+1):(M+1+number.of.tumor),(M+1):(M+1+number.of.tumor)]
  
  
  meta.result <- LogoddsMetaAnalysis(log.odds.icog,
                                     sigma.log.odds.icog,
                                     log.odds.onco,
                                     sigma.log.odds.onco)
  
  second.stage.logodds.meta <- meta.result[[1]]
  second.stage.sigma.meta <- meta.result[[2]]
  
  dis.result <- distance.fun(second.stage.logodds.meta,
                             second.stage.sigma.meta,
                             dis.trans)
  
  save(dis.result,file=paste0("/spin1/users/zhangh24/breast_cancer_data_analysis/known_SNPs/cluster_analysis/result",i1,".Rdata"))
  
  
}else{
  data2 <- as.data.frame(fread("./data/Onco_euro_v10_10232017.csv",header=T))
  #data2 <- read.csv("/data/zhangh20/breast_cancer/known_SNPs_anlysis/Onco_euro_v10_rs554219.csv",1)
  
  # names1 = colnames(data1)[27:206]
  #rm(data1)
  names2 = colnames(data2)
  
  idxi1 = which(names2=="rs554219")
  y.pheno.mis2 <- cbind(data2$Behaviour1,data2$PR_status1,data2$ER_status1,data2$HER2_status1)
  x.test.all.mis2 <- data2
  
  x.covar.mis2 <- data2[,5:14]
  x.all.mis2 <- as.matrix(cbind(x.test.all.mis2[,idxi1],x.covar.mis2))
  
  
  
  Heter.result.Onco = EMmvpolySelfDesign(y.pheno.mis2,x.self.design = x.all.mis2[,1,drop=F],z.design = z.design,baselineonly = NULL,additive = x.all.mis2[,2:11],pairwise.interaction = NULL,saturated = NULL,missingTumorIndicator = 888)
  z.standard <- Heter.result.Onco[[12]]
  M <- nrow(z.standard)
  number.of.tumor <- ncol(z.standard)
  log.odds.onco <- Heter.result.Onco[[1]][(M+1):(M+1+number.of.tumor)]
  sigma.log.odds.onco <- Heter.result.Onco[[2]][(M+1):(M+1+number.of.tumor),(M+1):(M+1+number.of.tumor)]
  
  
  
  
  second.stage.logodds.meta <-log.odds.onco
  second.stage.sigma.meta <- sigma.log.odds.onco
  
  
  
  
  
  dis.result <- distance.fun(second.stage.logodds.meta,
                             second.stage.sigma.meta,
                             dis.trans)
  
  save(dis.result,file=paste0("/spin1/users/zhangh24/breast_cancer_data_analysis/known_SNPs/cluster_analysis/result",i1,".Rdata"))
  
  
  
  
  
  
}





