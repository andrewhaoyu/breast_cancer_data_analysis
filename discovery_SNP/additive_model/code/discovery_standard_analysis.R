#install_github("andrewhaoyu/bc2",ref='development', args = c('--library="/home/zhangh24/R/x86_64-pc-linux-gnu-library/3.4"'))
#install_github("andrewhaoyu/bc2", args = c('--library="/home/zhangh24/R/x86_64-pc-linux-gnu-library/3.4"'))
#install_github("andrewhaoyu/bc2", args = c('--library="/home/zhangh24/R/x86_64-pc-linux-gnu-library/3.4"'))
###1 represent Icog
###2 represent Onco

rm(list=ls())
#commandarg <- commandArgs(trailingOnly=F)
#myarg <- commandarg[length(commandarg)]
#myarg <- sub("-","",myarg)
#i1 <- as.numeric(myarg)
# args = commandArgs(trailingOnly = T)
# i1 = as.numeric(args[[1]])

setwd("/data/zhangh24/breast_cancer_data_analysis/")

library(readr)
library(devtools)
library(CompQuadForm)
library(bc2)
library(data.table)
#icog.julie <- fread("/data/zhangh24/breast_cancer_data_analysis/data/Julie_snp_icog.csv")
#icog.julie <- icog.julie[,-1]
discovery.snp.icog <- fread("/data/zhangh24/breast_cancer_data_analysis/data/discovery_icog_data.csv",header=T)
#onco.julie <- fread("/data/zhangh24/breast_cancer_data_analysis/data/Julie_snp_onco.csv")
#onco.julie <- onco.julie[,-1]
discovery.snp.onco <- fread("/data/zhangh24/breast_cancer_data_analysis/data/discovery_onco_data.csv")
x.test.all.mis1 <- as.data.frame(discovery.snp.icog)
x.test.all.mis2 <- as.data.frame(discovery.snp.onco)

p.value <- rep(0,19)
data1 <- fread("./data/iCOGS_euro_v10_10232017.csv",header=T)
data1 <- as.data.frame(data1)
data2 <- fread("./data/Onco_euro_v10_10232017.csv",header=T)
data2 <- as.data.frame(data2)


for(i1 in 1:19){
  if(i1<=18){
 
    y.pheno.mis1 <- cbind(data1$Behaviour1,data1$ER_status1,data1$PR_status1,data1$HER2_status1,data1$Grade1)
    colnames(y.pheno.mis1) = c("Behavior","ER","PR","HER2","Grade")
    x.covar.mis1 <- cbind(data1[,c(5:14)],as.factor(data1[,3]))
    pc <- as.matrix(data1[,c(5:14)])
    country <- as.factor(data1[,3])
    gene <-   x.test.all.mis1[,i1]                     
    
    x.all.mis1 <- cbind(x.test.all.mis1[,i1],x.covar.mis1)
    
    colnames(x.all.mis1)[1] <- "gene"
    
    
    model1 <- glm(y.pheno.mis1[,1]~gene+pc+country,family = binomial(link = 'logit'))
    coeff.icog <- coef(model1)[2]
    var.icog <- (summary(model1)$coefficient[2,2])^2
    
    
    y.pheno.mis2 <- cbind(data2$Behaviour1,data2$ER_status1,data2$PR_status1,data2$HER2_status1,data2$Grade1)
    colnames(y.pheno.mis2) = c("Behaviour","ER",
                               "PR","HER2","Grade")
    
   
    x.covar.mis2 <- data2[,c(5:14)]
    
    
    
   
    pc <- as.matrix(data2[,c(5:14)])
    country <- as.factor(data2[,4])
    gene <-   x.test.all.mis2[,i1]    
      model2 <- glm(y.pheno.mis2[,1]~gene+pc+country,family = binomial(link = 'logit'))
    coeff.onco <- coef(model2)[2]
    var.onco <- (summary(model2)$coefficient[2,2])^2
    
    meta.result <- LogoddsMetaAnalysis(  coeff.icog,
                                         var.icog,
                                         coeff.onco,
                                         var.onco)
    
    coeff.meta <- meta.result[[1]]
    var.meta <- meta.result[[2]]
    p.value[i1] <- as.numeric(DisplaySecondStageTestResult(coeff.meta,var.meta)[2])  
    
    
    
    
    
    

    
    
  }else{
    data2 <- fread("./data/Onco_euro_v10_10232017.csv",header=T)
    data2 <- as.data.frame(data2)
    names2 = colnames(data2)
    y.pheno.mis2 <- cbind(data2$Behaviour1,data2$ER_status1,data2$PR_status1,data2$HER2_status1,data2$Grade1)
    colnames(y.pheno.mis2) = c("Behaviour","ER",
                               "PR","HER2","Grade")
    idxi1 = 29
    x.covar.mis2 <- data2[,c(5:14)]
    pc <- as.matrix(data2[,c(5:14)])
    country <- as.factor(data2[,4])
    gene <-   x.test.all.mis2[,i1]    
    model2 <- glm(y.pheno.mis2[,1]~gene+pc+country,family = binomial(link = 'logit'))
    coeff.onco <- coef(model2)[2]
    var.onco <- (summary(model2)$coefficient[2,2])^2
    
    
    
    coeff.meta <- coeff.onco
    var.meta <- var.onco
    p.value[i1] <- as.numeric(DisplaySecondStageTestResult(coeff.meta,var.meta)[2]  )
    
  
    
    
  }
  
  
  
  
}

write.csv(p.value,file= "./discovery_SNP/result/standard_pvalue.csv")
