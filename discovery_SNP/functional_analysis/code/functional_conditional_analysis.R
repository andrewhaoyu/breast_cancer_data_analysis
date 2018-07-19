#install_github("andrewhaoyu/bc2", args = c('--library="/home/zhangh24/R/x86_64-pc-linux-gnu-library/3.5"'))
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
#############get the conditional analysis results for the SNPs with nearby known SNPs
load(paste0("./discovery_SNP/functional_analysis/result/functional_snp_conditional.Rdata"))
load("./discovery_SNP/functional_analysis/result/ICOG/functional_conditional_icog_snpvalue.Rdata")
load("./discovery_SNP/functional_analysis/result/ONCO/functional_conditional_onco_snpvalue.Rdata")
new_filter <- fread("/spin1/users/zhangh24/breast_cancer_data_analysis/data/Filter_based_on_function.csv",header=T)
idx.temp <- which(is.na(functional_snp_conditional$SNP.ICOG))
num = nrow(functional_snp_conditional)
temp = diff(functional_snp_conditional$position)
size = 1000
start.end <- startend(num,size,i1)
start = start.end[1]
end  = start.end[2]

CheckKnow <- function(new_filter,position,chr){
  n <- nrow(new_filter)
  for(i in 1:n){
    chr.new.filter <- new_filter[i,3]
    position.new.filter <- new_filter[i,2]
    if(chr.new.filter==chr&position>=(position.new.filter-500000)&
       position<=(position.new.filter+500000)){
      return(i)
    }
  }
}
DisKnown <- function(dis.idx){
  if(dis.idx==1){
    return(5)
  }else if(dis.idx==2){
    return(c(9,10))
  }else if(dis.idx==3){
    return(157)
  }else{
    return(c(170,171))
  }
}

MixedEffectModel <- function(y.pheno.mis,
                             G,
                             x_covar
){
  z.standard <- GenerateZstandard(y.pheno.mis)
  M <- nrow(z.standard)
  z.design.fixed <- cbind(rep(1,M),z.standard[,1])
  z.design.random <-z.standard[,2:ncol(z.standard)]
  score.test.support.fixed <- ScoreTestSupportMixedModel(
    y.pheno.mis,
    baselineonly = NULL,
    additive = as.matrix(x_covar),
    pairwise.interaction = NULL,
    saturated = NULL,
    missingTumorIndicator = 888
  )
  score.test.fixed<- ScoreTestMixedModel(y=y.pheno.mis,
                                         x=as.matrix(G),
                                         z.design=z.design.fixed,
                                         score.test.support=score.test.support.fixed,
                                         missingTumorIndicator=888
  )
  score.fixed <- score.test.fixed[[1]]
  infor.fixed <- score.test.fixed[[2]]
  
  score.test.support.random <- ScoreTestSupportMixedModelSelfDesign(
    y.pheno.mis,
    x.self.design  = G,
    z.design = z.design.fixed,
    additive =  as.matrix(x_covar),
    pairwise.interaction = NULL,
    saturated = NULL,
    missingTumorIndicator = 888
  )
  score.test.random<- ScoreTestMixedModel(y=y.pheno.mis,
                                          x=as.matrix(G),
                                          z.design=z.design.random,
                                          score.test.support=score.test.support.random,
                                          missingTumorIndicator=888)
  score.random <- score.test.random[[1]]
  infor.random <- score.test.random[[2]]
  # p_mglobal <- DisplayMixedScoreTestResult(score.fixed,
  #                                          infor.fixed,
  #                                          score.random,
  #                                          infor.random)[1]
  return(list(score.fixed,
              infor.fixed,
              score.random,
              infor.random))
}








p.value.result <- rep(0,(end-start+1))
temp = 1


for(i in start:end){
  ##########check out which discovery SNP is nearby
  position = functional_snp_conditional$position[i]
  chr = functional_snp_conditional$CHR[i]
  dis.idx <- CheckKnow(new_filter,position,chr)
  known.idx <- DisKnown(dis.idx)
  SNP.ICOG <- functional_snp_conditional$SNP.ICOG[i]
  if(is.na(SNP.ICOG)!=T){
    data1 <- fread("./data/iCOGS_euro_v10_10232017.csv",header=T)
    data1 <- as.data.frame(data1)
    y.pheno.mis1 <- cbind(data1$Behaviour1,data1$ER_status1,data1$PR_status1,data1$HER2_status1,data1$Grade1)
    colnames(y.pheno.mis1) = c("Behavior","ER","PR","HER2","Grade")
    # Grade1.fake <- data1$Grade1
    # Grade1.fake[data1$Grade1==2|data1$Grade1==3] <- 1
    # Grade1.fake[data1$Grade1==1] <- 0
    #y.pheno.mis1 <- cbind(data1$Behaviour1,data1$PR_status1,data1$ER_status1,data1$HER2_status1,Grade1.fake)
    # y.pheno.mis1 <- cbind(data1$Behaviour1,data1$PR_status1,data1$ER_status1,data1$HER2_status1)
    
    x.known1 <- data1[,c(27:203)]
    
    ###pc1-10 and age
    x.covar.mis1 <- data1[,c(5:14,204)]
    age <- data1[,204]
    idx.complete1 <- which(age!=888)
    y.pheno.mis1 <- y.pheno.mis1[idx.complete1,]
    x.covar.mis1 <- x.covar.mis1[idx.complete1,]
    x.known1 <- x.known1[idx.complete1,]
    x.covar.con1 <- cbind(x.covar.mis1,x.known1[,known.idx])
    x.test1 <- snpvalue.result.icog[idx.complete1,i]
    Mixed.model1 <- MixedEffectModel(y.pheno.mis1,
                                     x.test1,
                                     x.covar.con1)
    score.fixed1 <- Mixed.model1[[1]]
    infor.fixed1 <- Mixed.model1[[2]]
    score.random1 <- Mixed.model1[[3]]
    infor.random1 <- Mixed.model1[[4]]
    
    
    data2 <- fread("./data/Onco_euro_v10_10232017.csv",header=T)
    data2 <- as.data.frame(data2)
    y.pheno.mis2 <- cbind(data2$Behaviour1,data2$ER_status1,data2$PR_status1,data2$HER2_status1,data2$Grade1)
    #y.pheno.mis2 <- cbind(data2$Behaviour1,data2$PR_status1,data2$ER_status1,data2$HER2_status1)
    colnames(y.pheno.mis2) = c("Behaviour","ER",
                               "PR","HER2","Grade")
    
    x.known2 <- data2[,c(27:203,205)]
    colnames(x.known2)[known.idx]
    x.covar.mis2 <- data2[,c(5:8,204)]
    
    ages <- data2[,204]
    idx.complete2 <- which(ages!=888)
    y.pheno.mis2 <- y.pheno.mis2[idx.complete2,]
    x.covar.mis2 <- x.covar.mis2[idx.complete2,]
    x.known2 <- x.known2[idx.complete2,]
    x.covar.con2 <- cbind(x.covar.mis2,x.known2[,known.idx])
    x.test2 <- snpvalue.result.onco[idx.complete2,i]
    Mixed.model2 <- MixedEffectModel(y.pheno.mis2,
                                     x.test2,
                                     x.covar.con2)
    score.fixed2 <- Mixed.model2[[1]]
    infor.fixed2 <- Mixed.model2[[2]]
    score.random2 <- Mixed.model2[[3]]
    infor.random2 <- Mixed.model2[[4]]
    
    fixed.meta <- ScoreMetaAnalysis(score.fixed1,
                      infor.fixed1,
                      score.fixed2,
                      infor.fixed2)
    score.fixed.meta <- fixed.meta[[1]]
    infor.fixed.meta <- fixed.meta[[2]]
    
    random.meta <- ScoreMetaAnalysis(score.random1,
                                    infor.random1,
                                    score.random2,
                                    infor.random2)
    score.random.meta <- random.meta[[1]]
    infor.random.meta <- random.meta[[2]]
    
    p.value.result[temp] = 
      DisplayMixedScoreTestResult(score.fixed.meta,
                                  infor.fixed.meta,
                                  score.random.meta,
                                  infor.random.meta)[1]
    temp = temp+1
    
  }else{
    
    data2 <- fread("./data/Onco_euro_v10_10232017.csv",header=T)
    data2 <- as.data.frame(data2)
    y.pheno.mis2 <- cbind(data2$Behaviour1,data2$ER_status1,data2$PR_status1,data2$HER2_status1,data2$Grade1)
    #y.pheno.mis2 <- cbind(data2$Behaviour1,data2$PR_status1,data2$ER_status1,data2$HER2_status1)
    colnames(y.pheno.mis2) = c("Behaviour","ER",
                               "PR","HER2","Grade")
    
    x.known2 <- data2[,c(27:203,205)]
    colnames(x.known2)[known.idx]
    x.covar.mis2 <- data2[,c(5:8,204)]
    
    ages <- data2[,204]
    idx.complete2 <- which(ages!=888)
    y.pheno.mis2 <- y.pheno.mis2[idx.complete2,]
    x.covar.mis2 <- x.covar.mis2[idx.complete2,]
    x.known2 <- x.known2[idx.complete2,]
    x.covar.con2 <- cbind(x.covar.mis2,x.known2[,known.idx])
    x.test2 <- snpvalue.result.onco[idx.complete2,i]
    Mixed.model2 <- MixedEffectModel(y.pheno.mis2,
                                     x.test2,
                                     x.covar.con2)
    score.fixed2 <- Mixed.model2[[1]]
    infor.fixed2 <- Mixed.model2[[2]]
    score.random2 <- Mixed.model2[[3]]
    infor.random2 <- Mixed.model2[[4]]
    
   
    p.value.result[temp] = 
      DisplayMixedScoreTestResult(score.fixed2,
                                  infor.fixed2,
                                  score.random2,
                                  infor.random2)[1]
    temp = temp+1
  }
}

save(p.value.result,
     file = paste0("./discovery_SNP/functional_analysis/result/ICOG/p.value.result",i1,".Rdata"))


