#Goal: estimate the correlation between the tumor characteristics

setwd("/data/zhangh24/breast_cancer_data_analysis/")
library(readr)
library(devtools)
library(CompQuadForm)
library(bc2)
library(data.table)


  CountFun <- function(ER,PR,HE2,grade,idx){
    result1 <- c(table(ER[idx]),table(PR[idx]),table(HER2[idx]),
      table(grade[idx]))
    n1 <- sum(HER2[idx]==0&(ER[idx]==1|PR[idx]==1)&(grade[idx]==1|grade[idx]==2))
    n2 <- sum(HER2[idx]==1&(ER[idx]==1|PR[idx]==1))
    n3 <- sum(HER2[idx]==0&(ER[idx]==1|PR[idx]==1)&grade[idx]==3)
    n4 <- sum(HER2[idx]==1&ER[idx]==0&PR[idx]==0)
    n5 <- sum(HER2[idx]==0&ER[idx]==0&PR[idx]==0)
    return(c(result1,n1,n2,n3,n4,n5))
  }
  
  TNcount <- function(ER,PR,HER2,grade,age,study,study.names){
    study.include <- NULL
    for(i in 1:length(study.names)){
      idx <- which(study==study.names[i])
      if(checkcondition(ER,PR,HER2,grade,idx)==T){
        study.include <- c(study.include,study.names[i])
      }
    }
    idx <- which(study%in%study.include)
    ER = ER[idx]
    PR = PR[idx]
    HER2 = HER2[idx]
    grade = grade[idx]
    age = age[idx]
    
    idx.complete <- which(ER!=888&
                          PR!=888&
                            HER2!=888&
                            grade!=888)
    new.data <- data.frame(ER[idx.complete],
                           PR[idx.complete],
                           HER2[idx.complete],
                           grade[idx.complete])
    correlation.result <- cor(new.data)
    
    
    
    
    #missing age people
    idx.mis <- which(age==888)
    
    result.mis.age <- CountFun(ER,PR,HER2,grade,idx.mis)
    result.age <- matrix(0,18,6)
    age.cut.point <- c(16,40,50,60,70,100)
    for(i in 1:5){
    idx <- which(age>=age.cut.point[i]&
                   age<age.cut.point[i+1])  
    result.age[,i] <- CountFun(ER,PR,HER2,grade,idx)
    }
    result.age[,6] <- rowSums(result.age[,1:5])
    
    result.portion <- matrix(0,18,6)
    for(i in 1:6){
      result.portion[,i] <-   round(result.age[,i]/result.age[,6],2)
    }
    result.portion <- result.portion*100
    result.portion <- matrix(paste0("(",result.portion,"%",")"),ncol=6)
    
    result.age.portion <- matrix(paste0(result.age," ",result.portion),ncol=6)
    
    
    
    return(list(result.mis.age,result.age,result.age.portion,correlation.result))
  }
  
  
checkcondition <- function(ER,PR,HER2,grade,idx){
  if(all(names(table(ER[idx]))=="888")|
     all(names(table(PR[idx]))=="888")|
     all(names(table(HER2[idx]))=="888")|
     all(names(table(grade[idx]))=="888")){
    return(F)
  }else{
    return(T)
  }
}
  
  
 
data1 <- fread("./data/iCOGS_euro_v10_10232017.csv",header=T)
data1 <- as.data.frame(data1)
y.pheno.mis1 <- cbind(data1$Behaviour1,data1$ER_status1,data1$PR_status1,data1$HER2_status1,data1$Grade1)
  colnames(y.pheno.mis1) = c("Behavior","ER","PR","HER2","Grade")
  age1 <- data1$age
  study1 <- data1$study
  data2 <- fread("./data/Onco_euro_v10_10232017.csv",header=T)
  data2 <- as.data.frame(data2)
  y.pheno.mis2 <- cbind(data2$Behaviour1,data2$ER_status1,data2$PR_status1,data2$HER2_status1,data2$Grade1)
  colnames(y.pheno.mis2) <- c("Behavior","ER","PR","HER2","Grade")
  age2 <- data2$age
  study2 <- data2$study
  y.pheno.mis <- rbind(y.pheno.mis1,y.pheno.mis2)
  casecontrol <- y.pheno.mis[,1]
  idx.case <- which(casecontrol==1)
  age <- c(age1,age2)[idx.case]
  study <- c(study1,study2)[idx.case]
  study.names <- names(table(study[idx.case]))
  ER <- y.pheno.mis[idx.case,2]
  PR <- y.pheno.mis[idx.case,3]
  HER2 <- y.pheno.mis[idx.case,4]
  grade <- y.pheno.mis[idx.case,5]
  temp.result <- TNcount(ER,PR,HER2,grade,age,study,study.names)
  result.age <- temp.result[[3]]
  correlation.result <- temp.result[[4]]
  write.csv(result.age,file = "./known_SNPs/known_SNPs_analysis_G_revised/tumor_correlation/result/result_age.csv")
  write.csv(correlation.result,file= "./known_SNPs/known_SNPs_analysis_G_revised/tumor_correlation/result/tumor_correlation.csv")
  
  
  
  
     # idx.complete2 <- which(y.pheno.mis2[,1]==1&
  #                          y.pheno.mis2[,2]!=888&
  #                          y.pheno.mis2[,3]!=888&
  #                          y.pheno.mis2[,4]!=888&
  #                          y.pheno.mis2[,5]!=888)
  # y.pheno.complete2 <- y.pheno.mis2[idx.complete2,]
  
  
  # colnames(sample.size.result2) <- c("study","sample size",
  #                                    "Luminal A-like",
  #                                    "Luminal B-like",
  #                                    "Luminal HER2-negative-like",
  #                                    "HER2 enriched-like",
  #                                    "TN")
  # sample.size.reuslt.all <- rbind(sample.size.result,
  #                                 sample.size.result2)
  # write.csv(sample.size.reuslt.all,
  #           file = "./known_SNPs/known_SNPs_analysis_G_revised/tumor_correlation/result/sample_size_by_study.csv")  
  # 
  
  # table(y.pheno.complete2[,1])
  # table(y.pheno.complete2[,2])
  # table(y.pheno.complete2[,3])
  # table(y.pheno.complete2[,4])
  # table(y.pheno.complete2[,5])
  # 
  # 
  
  
  
  
  
  # y.pheno.complete.all <- rbind(y.pheno.complete1,
  #                               y.pheno.complete2)  
  # correlation.result <- cor(y.pheno.complete.all[,2:5])  
 
  
