generate_second_stage_parameter_names = function(tumor_characteristics){
  result = c("baseline effect (95%CI)",
             "P_value for baseline effect")
  for(i in 1:length(tumor_characteristics)){
    result = c(result,paste0(tumor_characteristics[i]," main effect(95%CI)"),
               paste0(tumor_characteristics[i]," main effect P_Value"))
  }
  result = c(result,"Wald global test p value",
             "Wald global heterogneity test p value",
             "Score global test p value",
             "Mixed Model global test p value",
             "Mixed Model global heterogeneity test p value",
             "loglikelihood",
             "AIC")
  return(result)
}
change_binary_to_negative_positive = function(x){
  if(x==0){
    return("-")
  }else{
    return("+")
  }
}
generate_first_stage_parameter_names = function(tumor_characteristics,z_standard){
  max.z_standard = apply(z_standard,2,max)
  idx.not.binary = which(max.z_standard!=1)
  idx.binary = which(max.z_standard==1)
  result= NULL
  for(i in 1:nrow(z_standard)){
    names_each_row = NULL
    for(j in 1:ncol(z_standard)){
      if(j%in%idx.binary){
        temp = paste0(tumor_characteristics[j],change_binary_to_negative_positive(z_standard[i,j]))
      }else{
        temp = paste0(tumor_characteristics[j],z_standard[i,j])
      }
      names_each_row = paste0(names_each_row,temp)
    }
    result = c(result,paste0(names_each_row," OR (95%CI)"),paste("P_value for OR of ",names_each_row))
  }
  return(result)
}

freq <- function(x){
  sum(x)/(length(x)*2)
}
library(bc2)
##tp53
setwd("/data/zhangh24/breast_cancer_data_analysis/")
library(data.table)
data1 <- as.data.frame(fread("./data/iCOGS_euro_v10_05242017.csv",header=T))
data2 <- as.data.frame(fread("./data/Onco_euro_v10_05242017.csv",header=T))
x.test.all.mis.icog = as.data.frame(fread("./data/LD.position.knownSNP.pruning.SNP.icog.value.csv",header=T))
x.test.all.mis.onco = as.data.frame(fread("./data/LD.position.knownSNP.pruning.SNP.onco.value.csv",header=T))
y.pheno.mis1 <- cbind(data1$Behaviour1,data1$ER_status1,data1$PR_status1,data1$HER2_status1,data1$Grade1)
y.pheno.mis.2 <- cbind(data2$Behaviour1,data2$ER_status1,data2$PR_status1,data2$HER2_status1,
                          data2$Grade1)
snpname <- "chr17_7571752_G_T"
SNP.id.icog <- colnames(x.test.all.mis.icog)
idx.icog <- grep(snpname,SNP.id.icog)
i1 <- idx.icog
x.covar.mis1 <- data1[,5:14]
snpvalue <- x.test.all.mis.icog[,i1]
x.all.mis1 <- cbind(snpvalue,x.covar.mis1)

Heter.result.Icog = EMmvpoly(y.pheno.mis1,baselineonly = NULL,additive = x.all.mis1,pairwise.interaction = NULL,saturated = NULL,missingTumorIndicator = 888)
z.standard <- Heter.result.Icog[[12]]
z.additive.design <- as.matrix(cbind(1,z.standard))
M <- nrow(z.standard)
number.of.tumor <- ncol(z.standard)
log.odds.icog <- Heter.result.Icog[[1]][(M+1):(M+1+number.of.tumor)]

sigma.log.odds.icog <- Heter.result.Icog[[2]][(M+1):(M+1+number.of.tumor),(M+1):(M+1+number.of.tumor)]
beta.icog <- z.additive.design%*%log.odds.icog
beta.sigma.icog <- z.additive.design%*%sigma.log.odds.icog%*%t(z.additive.design)
loglikelihood.icog <- Heter.result.Icog[[8]]
rm(Heter.result.Icog)

score.test.support.icog <- ScoreTestSupport(
  y.pheno.mis1,
  baselineonly = NULL,
  additive = x.all.mis1[,2:ncol(x.all.mis1)],
  pairwise.interaction = NULL,
  saturated = NULL,
  missingTumorIndicator = 888
)
score.test.icog<- ScoreTest(y=y.pheno.mis1,
                            x=x.all.mis1[,1,drop=F],
                            second.stage.structure="additive",
                            score.test.support=score.test.support.icog,
                            missingTumorIndicator=888)
z.design.score.baseline <- matrix(rep(1,M),ncol=1)
z.design.score.casecase <-z.standard
z.design.score.baseline.ER <- cbind(z.design.score.baseline,z.standard[,1])
z.design.score.casecase.ER <- z.standard[,2:ncol(z.standard)]

score.icog <- score.test.icog[[1]]
infor.icog <- score.test.icog[[2]]
score.test.support.icog.baseline <- score.test.support.icog
#rm(score.test.support.icog)
rm(score.test.icog)  

score.test.icog.baseline<- ScoreTestSelfDesign(y=y.pheno.mis1,
                                               x=x.all.mis1[,1,drop=F],
                                               z.design=z.design.score.baseline,
                                               score.test.support=score.test.support.icog.baseline,
                                               missingTumorIndicator=888)

score.icog.baseline <- score.test.icog.baseline[[1]]
infor.icog.baseline <- score.test.icog.baseline[[2]]
rm(score.test.support.icog.baseline)  
rm(score.test.icog.baseline)
gc()  
score.test.support.icog.casecase <- ScoreTestSupport(
  y.pheno.mis1,
  baselineonly = x.all.mis1[,1,drop=F],
  additive = x.all.mis1[,2:ncol(x.all.mis1)],
  pairwise.interaction = NULL,
  saturated = NULL,
  missingTumorIndicator = 888
)
score.test.icog.casecase<- ScoreTestSelfDesign(y=y.pheno.mis1,
                                               x=x.all.mis1[,1,drop=F],
                                               z.design=z.design.score.casecase,
                                               score.test.support=score.test.support.icog.casecase,
                                               missingTumorIndicator=888)

score.icog.casecase <- score.test.icog.casecase[[1]]
infor.icog.casecase <- score.test.icog.casecase[[2]]

z.design.score.baseline.heterER <- z.standard[,1,drop=F]


score.icog.baseline.heterER <- score.icog.casecase[1]
infor.icog.baseline.heterER <- infor.icog.casecase[1,1]



rm(score.test.support.icog.casecase)
rm(score.test.icog.casecase)  
gc()  


score.test.icog.baseline.ER<- ScoreTestSelfDesign(y=y.pheno.mis1,
                                                  x=x.all.mis1[,1,drop=F],
                                                  z.design=z.design.score.baseline.ER,
                                                  score.test.support=score.test.support.icog,
                                                  missingTumorIndicator=888)

score.icog.baseline.ER <- score.test.icog.baseline.ER[[1]]
infor.icog.baseline.ER <- score.test.icog.baseline.ER[[2]]

score.test.support.icog.casecase.ER <- ScoreTestSupportSelfDesign(
  y.pheno.mis1,
  x.self.design  = x.all.mis1[,1,drop=F],
  z.design = z.design.score.baseline.ER,
  additive = x.all.mis1[,2:ncol(x.all.mis1)],
  pairwise.interaction = NULL,
  saturated = NULL,
  missingTumorIndicator = 888
)

score.test.icog.casecase.ER<- ScoreTestSelfDesign(y=y.pheno.mis1,
                                                  x=x.all.mis1[,1,drop=F],
                                                  z.design=z.design.score.casecase.ER,
                                                  score.test.support=score.test.support.icog.casecase.ER,
                                                  missingTumorIndicator=888)

score.icog.casecase.ER <- score.test.icog.casecase.ER[[1]]
infor.icog.casecase.ER <- score.test.icog.casecase.ER[[2]]
rm(score.test.support.icog.casecase.ER)
rm(score.test.icog.casecase.ER)
gc()




















idxi1 = i1

y.pheno.mis2 <- cbind(data2$Behaviour1,data2$ER_status1,data2$PR_status1,data2$HER2_status1,data2$Grade1)
#y.pheno.mis2 <- cbind(data2$Behaviour1,data2$PR_status1,data2$ER_status1,data2$HER2_status1)
colnames(y.pheno.mis2) = c("Behaviour","ER",
                           "PR","HER2","Grade")

#snpvalue <- x.test.all.mis.onco[,i1]
x.covar.mis2 <- data2[,5:14]
#x.all.mis2 <- as.matrix(cbind(x.test.all.mis.onco[,idxi1],x.covar.mis2))
snpvalue <- x.test.all.mis.onco[,idxi1]
x.all.mis2 <- cbind(snpvalue,x.covar.mis2)
colnames(x.all.mis2)[1] = "gene"




Heter.result.Onco = EMmvpoly(y.pheno.mis2,baselineonly = NULL,additive = x.all.mis2,pairwise.interaction = NULL,saturated = NULL,missingTumorIndicator = 888)
z.standard <- Heter.result.Onco[[12]]
M <- nrow(z.standard)
number.of.tumor <- ncol(z.standard)
log.odds.onco <- Heter.result.Onco[[1]][(M+1):(M+1+number.of.tumor)]
nparm <- length(Heter.result.Onco)
sigma.log.odds.onco <- Heter.result.Onco[[2]][(M+1):(M+1+number.of.tumor),(M+1):(M+1+number.of.tumor)]
beta.onco <- z.additive.design%*%log.odds.onco
beta.sigma.onco <- z.additive.design%*%sigma.log.odds.onco%*%t(z.additive.design)
loglikelihood.onco <- Heter.result.Onco[[8]]
#rm(Heter.result.Onco)  

score.test.support.onco <- ScoreTestSupport(
  y.pheno.mis2,
  baselineonly = NULL,
  additive = x.all.mis2[,2:ncol(x.all.mis2)],
  pairwise.interaction = NULL,
  saturated = NULL,
  missingTumorIndicator = 888
)
score.test.onco<- ScoreTest(y=y.pheno.mis2,
                            x=x.all.mis2[,1,drop=F],
                            second.stage.structure="additive",
                            score.test.support=score.test.support.onco,
                            missingTumorIndicator=888)
z.design.score.baseline <- matrix(rep(1,M),ncol=1)
z.design.score.casecase <-z.standard
z.design.score.baseline.ER <- cbind(z.design.score.baseline,z.standard[,1])
z.design.score.casecase.ER <- z.standard[,2:ncol(z.standard)]

score.onco <- score.test.onco[[1]]
infor.onco <- score.test.onco[[2]]
score.test.support.onco.baseline <- score.test.support.onco
#rm(score.test.support.onco)
rm(score.test.onco)

score.test.onco.baseline<- ScoreTestSelfDesign(y=y.pheno.mis2,
                                               x=x.all.mis2[,1,drop=F],
                                               z.design=z.design.score.baseline,
                                               score.test.support=score.test.support.onco.baseline,
                                               missingTumorIndicator=888)

score.onco.baseline <- score.test.onco.baseline[[1]]
infor.onco.baseline <- score.test.onco.baseline[[2]]
rm(score.test.support.onco.baseline)
rm(score.test.onco.baseline)
gc()
score.test.support.onco.casecase <- ScoreTestSupport(
  y.pheno.mis2,
  baselineonly = x.all.mis2[,1,drop=F],
  additive = x.all.mis2[,2:ncol(x.all.mis2)],
  pairwise.interaction = NULL,
  saturated = NULL,
  missingTumorIndicator = 888
)
score.test.onco.casecase<- ScoreTestSelfDesign(y=y.pheno.mis2,
                                               x=x.all.mis2[,1,drop=F],
                                               z.design=z.design.score.casecase,
                                               score.test.support=score.test.support.onco.casecase,
                                               missingTumorIndicator=888)

score.onco.casecase <- score.test.onco.casecase[[1]]
infor.onco.casecase <- score.test.onco.casecase[[2]]
z.design.score.baseline.heterER <- z.standard[,1,drop=F]


score.onco.baseline.heterER <-   score.onco.casecase[1]
infor.onco.baseline.heterER <- infor.onco.casecase[1,1]




rm(score.test.support.onco.casecase) 







score.test.onco.baseline.ER<- ScoreTestSelfDesign(y=y.pheno.mis2,
                                                  x=x.all.mis2[,1,drop=F],
                                                  z.design=z.design.score.baseline.ER,
                                                  score.test.support=score.test.support.onco,
                                                  missingTumorIndicator=888)

score.onco.baseline.ER <- score.test.onco.baseline.ER[[1]]
infor.onco.baseline.ER <- score.test.onco.baseline.ER[[2]]

score.test.support.onco.casecase.ER <- ScoreTestSupportSelfDesign(
  y.pheno.mis2,
  x.self.design  = x.all.mis2[,1,drop=F],
  z.design = z.design.score.baseline.ER,
  additive = x.all.mis2[,2:ncol(x.all.mis2)],
  pairwise.interaction = NULL,
  saturated = NULL,
  missingTumorIndicator = 888
)

score.test.onco.casecase.ER<- ScoreTestSelfDesign(y=y.pheno.mis2,
                                                  x=x.all.mis2[,1,drop=F],
                                                  z.design=z.design.score.casecase.ER,
                                                  score.test.support=score.test.support.onco.casecase.ER,
                                                  missingTumorIndicator=888)

score.onco.casecase.ER <- score.test.onco.casecase.ER[[1]]
infor.onco.casecase.ER <- score.test.onco.casecase.ER[[2]]

rm(score.test.support.onco.casecase.ER)
rm(score.test.onco.casecase.ER)
gc()


meta.result <- LogoddsMetaAnalysis(log.odds.icog,
                                   sigma.log.odds.icog,
                                   log.odds.onco,
                                   sigma.log.odds.onco)

second.stage.logodds.meta <- meta.result[[1]]
second.stage.sigma.meta <- meta.result[[2]]



test.result.second.wald <- DisplaySecondStageTestResult(second.stage.logodds.meta,second.stage.sigma.meta)


beta.meta <- z.additive.design%*%second.stage.logodds.meta
beta.sigma.meta <- z.additive.design%*%second.stage.sigma.meta%*%t(z.additive.design)

test.result.first.wald <- DisplayFirstStageTestResult(beta.meta,beta.sigma.meta)

meta.result.score <- ScoreMetaAnalysis(score.icog,infor.icog,
                                       score.onco,infor.onco)
score.meta <- meta.result.score[[1]]
infor.meta <- meta.result.score[[2]]

test.result.second.score <- DisplayFixedScoreTestResult(score.meta,infor.meta)

meta.result.score.baseline <- ScoreMetaAnalysis(score.icog.baseline,
                                                infor.icog.baseline,
                                                score.onco.baseline,
                                                infor.onco.baseline)
score.meta.baseline <- meta.result.score.baseline[[1]]
infor.meta.baseline <- meta.result.score.baseline[[2]]

meta.result.score.casecase <- ScoreMetaAnalysis(score.icog.casecase,
                                                infor.icog.casecase,
                                                score.onco.casecase,
                                                infor.onco.casecase)
score.meta.casecase <- meta.result.score.casecase[[1]]
infor.meta.casecase <- meta.result.score.casecase[[2]]

meta.result.score.baseline.heterER <- 
  ScoreMetaAnalysis(score.icog.baseline.heterER,
                    infor.icog.baseline.heterER,
                    score.onco.baseline.heterER,
                    infor.onco.baseline.heterER)
score.meta.baseline.heterER <- meta.result.score.baseline.heterER[[1]]  
infor.meta.baseline.heterER <- meta.result.score.baseline.heterER[[2]]


meta.result.score.baseline.ER <- ScoreMetaAnalysis(score.icog.baseline.ER,
                                                   infor.icog.baseline.ER,
                                                   score.onco.baseline.ER,
                                                   infor.onco.baseline.ER)
score.meta.baseline.ER <-   meta.result.score.baseline.ER[[1]]
infor.meta.baseline.ER <- meta.result.score.baseline.ER[[2]]

meta.result.score.casecase.ER <- ScoreMetaAnalysis(score.icog.casecase.ER,
                                                   infor.icog.casecase.ER,
                                                   score.onco.casecase.ER,
                                                   infor.onco.casecase.ER)
score.meta.casecase.ER <-   meta.result.score.casecase.ER[[1]]
infor.meta.casecase.ER <- meta.result.score.casecase.ER[[2]]




test.result.second.mixed <- DisplayMixedScoreTestResult(score.meta.baseline,
                                                        infor.meta.baseline,
                                                        score.meta.casecase,
                                                        infor.meta.casecase)  
test.result.second.mixed <- data.frame(t(test.result.second.mixed))



###########DisplayMixedScoreTestResult is set up for one fixed effect,
###########we need to adjust for the heterogeneity test
test.result.second.mixed.ER <- DisplayMixedScoreTestResult(score.meta.baseline.ER,
                                                           infor.meta.baseline.ER,
                                                           score.meta.casecase.ER,
                                                           infor.meta.casecase.ER)  
test.result.second.mixed.ER <- data.frame(t(test.result.second.mixed.ER))

test.result.second.mixed.ER.2 <- DisplayMixedScoreTestResult(score.meta.baseline.heterER,
                                                             infor.meta.baseline.heterER,
                                                             score.meta.casecase.ER,
                                                             infor.meta.casecase.ER)  
test.result.second.mixed.ER.2 <- data.frame(t(test.result.second.mixed.ER.2))

###Both ER and RANDOM EFFECT should be counted into heterogeneity
test.result.second.mixed.ER[1,2] <-  test.result.second.mixed.ER.2[1,1]


colnames(test.result.second.mixed) <- c("mixed model global test for association(baseline fixed)","mixed model global test for heterogeneity(baseline fixed)")
colnames(test.result.second.mixed.ER) <- c("mixed model global test for association(baseline and ER fixed)","mixed model global test for heterogeneity(baseline and ER fixed)")

loglikelihood <- loglikelihood.icog+loglikelihood.onco
AIC <- 2*length(Heter.result.Onco[[1]])-2*loglikelihood

heter.result <- list(data.frame(test.result.second.wald,test.result.second.score, test.result.second.mixed,test.result.second.mixed.ER,loglikelihood = loglikelihood,AIC=AIC),
                     data.frame(test.result.first.wald))
#save(heter.result,file=paste0("./known_SNPs/known_SNPs_analysis_G_revised/additive_model/result/heter_result_",i1,".Rdata"))



save(heter.result,file="/data/zhangh24/breast_cancer_data_analysis/discovery_SNP/tp53/result/heter.result")


