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

data1 <- read.csv("./data/iCOGS_euro_v10_05242017.csv",header=T)
data2 <- read.csv("./data/Onco_euro_v10_05242017.csv",header=T)
x.test.all.mis.icog = read.csv("./data/LD.position.knownSNP.pruning.SNP.icog.value.csv",header=T)
x.test.all.mis.onco = read.csv("./data/LD.position.knownSNP.pruning.SNP.onco.value.csv",header=T)
y.pheno.mis.icog <- cbind(data1$Behaviour1,data1$PR_status1,data1$ER_status1,data1$HER2_status1)
y.pheno.mis.onco <- cbind(data2$Behaviour1,data2$PR_status1,data2$ER_status1,data2$HER2_status1)
snpname <- "chr17_7571752_G_T"
SNP.id.icog <- colnames(x.test.all.mis.icog)
idx.icog <- grep(snpname,SNP.id.icog)
i1 <- idx.icog
x.covar.mis1 <- data1[,5:14]
snpvalue <- x.test.all.mis.icog[,i1]



idx.control <- which(y.pheno.mis.icog[,1]==0)
snpvalue.control <- snpvalue[idx.control]
freq.control <- freq(snpvalue.control)
idx.luA <- which(y.pheno.mis.icog[,1]==1&(y.pheno.mis.icog[,4]==0)&
                   (y.pheno.mis.icog[,2]==1|y.pheno.mis.icog[,3]==1))
snpvalue.luA <- snpvalue[idx.luA]
freq.LuA <- freq(snpvalue.luA)
idx.luB <- which(y.pheno.mis.icog[,1]==1&(y.pheno.mis.icog[,4]==1)&
                   (y.pheno.mis.icog[,2]==1|y.pheno.mis.icog[,3]==1))
snpvalue.luB <- snpvalue[idx.luB]
freq.luB <- freq(snpvalue.luB)
idx.her2E <- which(y.pheno.mis.icog[,1]==1&(y.pheno.mis.icog[,4]==1)&
                   (y.pheno.mis.icog[,2]==0&y.pheno.mis.icog[,3]==0))
snpvalue.her2E <- snpvalue[idx.her2E]
freq.her2E <- freq(snpvalue.her2E)
idx.tri <- which(y.pheno.mis.icog[,1]==1&(y.pheno.mis.icog[,4]==0)&
                   (y.pheno.mis.icog[,2]==0&y.pheno.mis.icog[,3]==0))
snpvalue.tri <- snpvalue[idx.tri]
freq.tri <- freq(snpvalue.tri)

print(round(c(freq.control,freq.LuA,freq.luB,freq.her2E,freq.tri),3))

x.all.mis1 <- as.matrix(cbind(x.test.all.mis.icog[,i1],x.covar.mis1))
y.pheno.mis1 <- cbind(data1$Behaviour1,data1$PR_status1,data1$ER_status1,data1$HER2_status1)
colnames(y.pheno.mis1) = c("Behavior","PR","ER","HER2")

Heter.result.Icog = EMmvpoly(y.pheno.mis1,baselineonly = NULL,additive = x.all.mis1,pairwise.interaction = NULL,saturated = NULL,missingTumorIndicator = 888)
M <- 8
z.standard <- Heter.result.Icog[[12]]
z.additive.design <- as.matrix(cbind(1,z.standard))
M <- nrow(z.standard)
number.of.tumor <- ncol(z.standard)
log.odds.icog <- Heter.result.Icog[[1]][(M+1):(M+1+number.of.tumor)]

sigma.log.odds.icog <- Heter.result.Icog[[2]][(M+1):(M+1+number.of.tumor),(M+1):(M+1+number.of.tumor)]
beta.icog <- z.additive.design%*%log.odds.icog
beta.sigma.icog <- z.additive.design%*%sigma.log.odds.icog%*%t(z.additive.design)
loglikelihood.icog <- Heter.result.Icog[[8]]
score.test.support.icog <- ScoreTestSupport(
  y.pheno.mis1,
  baselineonly = NULL,
  additive = x.all.mis1[,2:11],
  pairwise.interaction = NULL,
  saturated = NULL,
  missingTumorIndicator = 888
)
score.test.icog<- ScoreTest(y=y.pheno.mis1,
                            x=x.all.mis1[,1,drop=F],
                            second.stage.structure="additive",
                            score.test.support=score.test.support.icog,
                            missingTumorIndicator=888)
z.design.score.baseline <- matrix(rep(1,8),ncol=1)
z.design.score.casecase <-z.standard

score.icog <- score.test.icog[[1]]
infor.icog <- score.test.icog[[2]]

score.test.support.icog.baseline <- ScoreTestSupport(
  y.pheno.mis1,
  baselineonly = NULL,
  additive = x.all.mis1[,2:11],
  pairwise.interaction = NULL,
  saturated = NULL,
  missingTumorIndicator = 888
)
score.test.icog.baseline<- ScoreTestSelfDesign(y=y.pheno.mis1,
                                               x=x.all.mis1[,1,drop=F],
                                               z.design=z.design.score.baseline,
                                               score.test.support=score.test.support.icog.baseline,
                                               missingTumorIndicator=888)

score.icog.baseline <- score.test.icog.baseline[[1]]
infor.icog.baseline <- score.test.icog.baseline[[2]]

score.test.support.icog.casecase <- ScoreTestSupport(
  y.pheno.mis1,
  baselineonly = x.all.mis1[,1,drop=F],
  additive = x.all.mis1[,2:11],
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







idxi1 = i1

y.pheno.mis2 <- cbind(data2$Behaviour1,data2$PR_status1,data2$ER_status1,data2$HER2_status1)
#y.pheno.mis2 <- cbind(data2$Behaviour1,data2$PR_status1,data2$ER_status1,data2$HER2_status1)
colnames(y.pheno.mis2) = c("Behaviour","PR",
                           "ER","HER2")


x.covar.mis2 <- data2[,5:14]
x.all.mis2 <- as.matrix(cbind(x.test.all.mis.onco[,idxi1],x.covar.mis2))
colnames(x.all.mis2)[1] = "gene"

snpvalue <- x.test.all.mis.onco[,idxi1]

idx.control <- which(y.pheno.mis.onco[,1]==0)
snpvalue.control <- snpvalue[idx.control]
freq.control <- freq(snpvalue.control)
idx.luA <- which(y.pheno.mis.onco[,1]==1&(y.pheno.mis.onco[,4]==0)&
                   (y.pheno.mis.onco[,2]==1|y.pheno.mis.onco[,3]==1))
snpvalue.luA <- snpvalue[idx.luA]
freq.LuA <- freq(snpvalue.luA)
idx.luB <- which(y.pheno.mis.onco[,1]==1&(y.pheno.mis.onco[,4]==1)&
                   (y.pheno.mis.onco[,2]==1|y.pheno.mis.onco[,3]==1))
snpvalue.luB <- snpvalue[idx.luB]
freq.luB <- freq(snpvalue.luB)
idx.her2E <- which(y.pheno.mis.onco[,1]==1&(y.pheno.mis.onco[,4]==1)&
                     (y.pheno.mis.onco[,2]==0&y.pheno.mis.onco[,3]==0))
snpvalue.her2E <- snpvalue[idx.her2E]
freq.her2E <- freq(snpvalue.her2E)
idx.tri <- which(y.pheno.mis.onco[,1]==1&(y.pheno.mis.onco[,4]==0)&
                   (y.pheno.mis.onco[,2]==0&y.pheno.mis.onco[,3]==0))
snpvalue.tri <- snpvalue[idx.tri]
freq.tri <- freq(snpvalue.tri)

print(round(c(freq.control,freq.LuA,freq.luB,freq.her2E,freq.tri),3))



Heter.result.Onco = EMmvpoly(y.pheno.mis2,baselineonly = NULL,additive = x.all.mis2,pairwise.interaction = NULL,saturated = NULL,missingTumorIndicator = 888)
z.standard <- Heter.result.Onco[[12]]
M <- nrow(z.standard)
number.of.tumor <- ncol(z.standard)
log.odds.onco <- Heter.result.Onco[[1]][(M+1):(M+1+number.of.tumor)]
sigma.log.odds.onco <- Heter.result.Onco[[2]][(M+1):(M+1+number.of.tumor),(M+1):(M+1+number.of.tumor)]
beta.onco <- z.additive.design%*%log.odds.onco
beta.sigma.onco <- z.additive.design%*%sigma.log.odds.onco%*%t(z.additive.design)
loglikelihood.onco <- Heter.result.Onco[[8]]


score.test.support.onco <- ScoreTestSupport(
  y.pheno.mis2,
  baselineonly = NULL,
  additive = x.all.mis2[,2:11],
  pairwise.interaction = NULL,
  saturated = NULL,
  missingTumorIndicator = 888
)
score.test.onco<- ScoreTest(y=y.pheno.mis2,
                            x=x.all.mis2[,1,drop=F],
                            second.stage.structure="additive",
                            score.test.support=score.test.support.onco,
                            missingTumorIndicator=888)
z.design.score.baseline <- matrix(rep(1,8),ncol=1)
z.design.score.casecase <-z.standard

score.onco <- score.test.onco[[1]]
infor.onco <- score.test.onco[[2]]

score.test.support.onco.baseline <- ScoreTestSupport(
  y.pheno.mis2,
  baselineonly = NULL,
  additive = x.all.mis2[,2:11],
  pairwise.interaction = NULL,
  saturated = NULL,
  missingTumorIndicator = 888
)
score.test.onco.baseline<- ScoreTestSelfDesign(y=y.pheno.mis2,
                                               x=x.all.mis2[,1,drop=F],
                                               z.design=z.design.score.baseline,
                                               score.test.support=score.test.support.onco.baseline,
                                               missingTumorIndicator=888)

score.onco.baseline <- score.test.onco.baseline[[1]]
infor.onco.baseline <- score.test.onco.baseline[[2]]

score.test.support.onco.casecase <- ScoreTestSupport(
  y.pheno.mis2,
  baselineonly = x.all.mis2[,1,drop=F],
  additive = x.all.mis2[,2:11],
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
test.result.second.mixed <- DisplayMixedScoreTestResult(score.meta.baseline,
                                                        infor.meta.baseline,
                                                        score.meta.casecase,
                                                        infor.meta.casecase)  
test.result.second.mixed <- data.frame(t(test.result.second.mixed))

colnames(test.result.second.mixed) <- c("mixed model global test for association","mixed model global test for heterogeneity")

loglikelihood <- loglikelihood.icog+loglikelihood.onco
AIC <- 2*length(Heter.result.Onco[[1]])-2*loglikelihood

heter.result <- list(data.frame(test.result.second.wald,test.result.second.score, test.result.second.mixed,loglikelihood = loglikelihood,AIC=AIC),
                     data.frame(test.result.first.wald))



result <-  NULL
first.stage <- NULL


  
  result <- rbind(result,heter.result[[1]])
  first.stage <- rbind(first.stage,heter.result[[2]])


tumor.characteristics <- c("PR","ER","HER2")
generate_second_stage_parameter_names(tumor.characteristics)

colnames(result) <- generate_second_stage_parameter_names(tumor.characteristics)

result <- as.data.frame(result)

colnames(first.stage) = generate_first_stage_parameter_names(tumor.characteristics,z.standard)

write.csv(result,file="./discovery_SNP/tp53_rs78378222/result/erprher2_additive_model.xlsx")





z.design <- matrix(c(
  c(0,1,1,1,0,0,0,0),
  c(0,0,0,0,0,1,1,1),
  c(0,0,0,0,1,0,0,0),
  c(1,0,0,0,0,0,0,0)
),ncol=4)

colnames(z.design) <- c("Luminial A",
                        "Luminal B",
                        "HER2 Enriched",
                        "Triple Negative")

Heter.result.Icog = EMmvpolySelfDesign(y.pheno.mis1,x.self.design = x.all.mis1[,1,drop=F],z.design=z.design,baselineonly = NULL,additive = x.all.mis1[,2:11],pairwise.interaction = NULL,saturated = NULL,missingTumorIndicator = 888)
z.standard <- Heter.result.Icog[[12]]
z.additive.design <- as.matrix(cbind(1,z.standard))
M <- nrow(z.standard)
number.of.tumor <- ncol(z.standard)
log.odds.icog <- Heter.result.Icog[[1]][(M+1):(M+1+number.of.tumor)]

sigma.log.odds.icog <- Heter.result.Icog[[2]][(M+1):(M+1+number.of.tumor),(M+1):(M+1+number.of.tumor)]
beta.icog <- z.design%*%log.odds.icog
beta.sigma.icog <- z.design%*%sigma.log.odds.icog%*%t(z.design)
loglikelihood.icog <- Heter.result.Icog[[8]]

score.test.support.icog <- ScoreTestSupport(
  y.pheno.mis1,
  baselineonly = NULL,
  additive = x.all.mis1[,2:11],
  pairwise.interaction = NULL,
  saturated = NULL,
  missingTumorIndicator = 888
)



score.test.icog<- ScoreTestSelfDesign(y=y.pheno.mis1,
                                      x=x.all.mis1[,1,drop=F],
                                      z.design=z.design,
                                      score.test.support=score.test.support.icog,
                                      missingTumorIndicator=888)
z.design.score.baseline <- matrix(rep(1,8),ncol=1)
z.design.score.casecase <-z.standard

score.icog <- score.test.icog[[1]]
infor.icog <- score.test.icog[[2]]




#analysis for Onco Array
#data2 <- read.csv("./V10/Onco_euro_v10_05242017.csv",header=T)
y.pheno.mis2 <- cbind(data2$Behaviour1,data2$PR_status1,data2$ER_status1,data2$HER2_status1)
#y.pheno.mis2 <- cbind(data2$Behaviour1,data2$PR_status1,data2$ER_status1,data2$HER2_status1)
colnames(y.pheno.mis2) = c("Behaviour","PR",
                           "ER","HER2")

x.test.all.mis2 <- data2[,c(27:212)]
x.covar.mis2 <- data2[,5:14]
x.all.mis2 <- as.matrix(cbind(x.test.all.mis.onco[,idxi1],x.covar.mis2))
colnames(x.all.mis2)[1] = "gene"


Heter.result.Onco = EMmvpolySelfDesign(y.pheno.mis2,x.self.design = x.all.mis2[,1,drop=F],z.design = z.design,baselineonly = NULL,additive = x.all.mis2[,2:11],pairwise.interaction = NULL,saturated = NULL,missingTumorIndicator = 888)
z.standard <- Heter.result.Onco[[12]]
M <- nrow(z.standard)
number.of.tumor <- ncol(z.standard)
log.odds.onco <- Heter.result.Onco[[1]][(M+1):(M+1+number.of.tumor)]
sigma.log.odds.onco <- Heter.result.Onco[[2]][(M+1):(M+1+number.of.tumor),(M+1):(M+1+number.of.tumor)]
beta.onco <- z.design%*%log.odds.onco
beta.sigma.onco <- z.design%*%sigma.log.odds.onco%*%t(z.design)
loglikelihood.onco <- Heter.result.Onco[[8]]


score.test.support.onco <- ScoreTestSupport(
  y.pheno.mis2,
  baselineonly = NULL,
  additive = x.all.mis2[,2:11],
  pairwise.interaction = NULL,
  saturated = NULL,
  missingTumorIndicator = 888
)

score.test.onco<- ScoreTestSelfDesign(y=y.pheno.mis2,
                                      x=x.all.mis2[,1,drop=F],
                                      z.design= z.design,
                                      score.test.support=score.test.support.onco,
                                      missingTumorIndicator=888)

score.onco <- score.test.onco[[1]]
infor.onco <- score.test.onco[[2]]



meta.result <- LogoddsMetaAnalysis(log.odds.icog,
                                   sigma.log.odds.icog,
                                   log.odds.onco,
                                   sigma.log.odds.onco)

second.stage.logodds.meta <- meta.result[[1]]
second.stage.sigma.meta <- meta.result[[2]]



test.result.second.wald <- DisplaySecondStageTestResult(second.stage.logodds.meta,second.stage.sigma.meta,self.design=T)


beta.meta <- z.design%*%second.stage.logodds.meta
beta.sigma.meta <- z.design%*%second.stage.sigma.meta%*%t(z.design)

test.result.first.wald <- DisplayFirstStageTestResult(beta.meta,beta.sigma.meta)

meta.result.score <- ScoreMetaAnalysis(score.icog,infor.icog,
                                       score.onco,infor.onco)
score.meta <- meta.result.score[[1]]
infor.meta <- meta.result.score[[2]]

test.result.second.score <- ScoreGlobalMixedTestForAssoc(score.meta,infor.meta)
test.result.second.score <- t(test.result.second.score)


colnames(test.result.second.score) <- c("fixed effect score test global test for association","random effect model global test for association")

loglikelihood <- loglikelihood.icog+loglikelihood.onco
AIC <- 2*length(Heter.result.Onco[[1]])-2*loglikelihood

heter.result <- list(data.frame(test.result.second.wald,test.result.second.score,loglikelihood = loglikelihood,AIC=AIC),
                     data.frame(test.result.first.wald))


generate_self_design_second_stage_parameter_names = function(tumor_characteristics){
  result = NULL
  for(i in 1:length(tumor_characteristics)){
    result = c(result,paste0(tumor_characteristics[i]," odds ratio(95%CI)"),
               paste0(tumor_characteristics[i]," P_Value"))
  }
  result = c(result,"Wald global test p value",
             "Wald global heterogneity test p value",
             "Score global test p value",
             "Mixed Model global test p value",
             "loglikelihood",
             "AIC")
  return(result)
}

result <-  NULL
first.stage <- NULL





  result <- rbind(result,heter.result[[1]])
  first.stage <- rbind(first.stage,heter.result[[2]])


tumor.characteristics <- c("Luminal A","Luminal B","HER2 Enriched","Triple Neg")
generate_self_design_second_stage_parameter_names(tumor.characteristics)

colnames(result) <- generate_self_design_second_stage_parameter_names(tumor.characteristics)

result <- as.data.frame(result)
tumor.characteristics <- c("PR","ER","HER2")

colnames(first.stage) = generate_first_stage_parameter_names(tumor.characteristics,z.standard)

write.csv(result,file="./discovery_SNP/tp53_rs78378222/result/intrinsic_subtypes_erprher2.csv")





