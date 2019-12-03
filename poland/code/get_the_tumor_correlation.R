#goal: get the tumor correlation in PBCS

setwd("/data/zhangh24/breast_cancer_data_analysis/")
library(readr)
library(devtools)
library(CompQuadForm)
library(bc2)
library(data.table)
library(nnet)
Generatesubtypes<- function(ER,PR,HER2,Grade){
  n <- length(ER)
  subtypes <- rep("control",n)
  temp = 1
  for(i in 0:1){
    for(j in 0:1){
      for(k in 0:1){
        for(l in 1:3){
          
          idx <- which(ER==i&PR==j&HER2==k&Grade==l)
          if(length(idx)!=0){
            subtypes[idx] <- temp
            temp = temp+1  
          }
          
        }
      }
    }
  }
  
  subtypes <- factor(subtypes,levels=c("control",
                                       c(1:temp)))
  sum <- table(subtypes)
  idx.cat <- which(sum<=10)
  idx.remove <- which((subtypes%in%(unique(idx.cat)-1))==T)
  subtypes <- subtypes[-idx.remove]
  return(list(subtypes,idx.remove))
}




  data2 <- fread("./data/Onco_euro_v10_10232017.csv",header=T)
  data2 <- as.data.frame(data2)
   y.pheno.mis2 <- cbind(data2$Behaviour1,data2$ER_status1,data2$PR_status1,data2$HER2_status1,data2$Grade1)
  
  #y.pheno.mis2 <- cbind(data2$Behaviour1,data2$PR_status1,data2$ER_status1,data2$HER2_status1)
  colnames(y.pheno.mis2) = c("Behaviour","ER",
                             "PR","HER2","Grade")
  y.pheno.mis2 <- data.frame(y.pheno.mis2,
                             stringsAsFactors = F)
library(dplyr)
  y.pheno.complete <- y.pheno.mis2 %>% 
    filter(Behaviour==1&ER!=888&
             PR!=888&HER2!=888&
             Grade!=888)
  cor(y.pheno.complete[,2:5])
  