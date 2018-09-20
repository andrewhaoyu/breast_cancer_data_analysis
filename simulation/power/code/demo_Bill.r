###########simulate data with four tumor characteristics containing three binary tumor characteristics and one oridinal tumor characteristics.
###########one genotype with MAF 0.25 is simulated
###########one covariate with rnorm(n) simulated


SimulateDataPower <- function(theta_intercept,theta_test,theta_covar,n){
  a <- c(0,1)
  b <- c(0,1)
  c <- c(0,1)
  d <- c(1,2,3)
  ##number of the other covarites
  p_col <- 1
  z.standard <- as.matrix(expand.grid(a,b,c,d)) # orig
  M <- nrow(z.standard)
  #z <- as.matrix(expand.grid(a,b))
  
  
  #this z_design matrix is the second stage matrix 
  z.design.additive <- cbind(1,z.standard)
  M <- nrow(z.design.additive)
  additive.number <- 2
  additive.second.cat <- ncol(z.design.additive)
  total.covar.number <- 1+ additive.number
  z.all <- matrix(0,nrow=(M*total.covar.number),ncol = (M+                                                        additive.second.cat*additive.number))
  for(i in c("intercept",
             "additive"
  )){
    ##we always keep intercept as saturated model and to simply, we always use diagnonal matrix for intercept
    if(i=="intercept"){
      ###row start and column start point for this category
      row.start <- 0
      column.start <- 0
      for(j in 1:M){
        z.all[row.start+1+(j-1)*total.covar.number,(column.start+j)] = 1
      }
    }else if(i=="additive"){
      column.start <- M
      row.start <- 1
      if(additive.number!=0){
        for(j in 1:M){
          for(k in 1:additive.number){
            z.all[row.start+k+(j-1)*total.covar.number,
                  (column.start+(k-1)*additive.second.cat+1):
                    (column.start+k*additive.second.cat)] <- as.matrix(z.design.additive[j,])
          }
        }
      }
    }
  }
  
  
  K <- ncol(z.design.additive)
  
  # z <- kronecker(diag(2),z)
  #this theta is the true value
  theta <- c(theta_intercept,theta_test,theta_covar)
  
  #this is the true beta
  beta <- z.all%*%theta
  beta <- matrix(beta,nrow=p_col+2)
  #alpha <- c(0,rep(1,length(beta)-1))
  
  G <-  rbinom(n,2,0.25)
  x_covar <- rnorm(n)
  x_all <- cbind(1,G,x_covar) ##adding the intercept into the model
  
  predictor <- x_all%*%(beta)
  predictor <- cbind(0,predictor)
  #predictor <- sweep(predictor,2,alpha,"+")
  p <- exp(predictor)
  sp <- rowSums(p)
  #this standarize the probability sum into 1,since logit model:pr(D=1|predictor)=exp(predictor)/(1+exp(predictor))
  p <- sweep(p,1,sp,"/")
  y <- t(apply(p,1,function(x){rmultinom(1,1,x)}))
  
  y <- y[,-1] # this is the y matrix for the model
  y.case.control <- rowSums(y)
  y.tumor <- y%*%z.standard
  
  idx.control <- which(y.case.control==0)
  y.tumor[idx.control,] <- rep(NA,4)
  idx.case <- which(y.case.control==1)
  rate <- c(0.17,0.25,0.42,0.27)
  for(i in 1:4){
    idx.mis <-  sample(idx.case,size=length(idx.case)*rate[i])
    y.tumor[idx.mis,i] <- 888
  }
  y.pheno.mis <- cbind(y.case.control,y.tumor)
  return(list(y.pheno.mis=y.pheno.mis,
              G=G,
              x_covar=x_covar))  
}







install_github("andrewhaoyu/bc2")
library(bc2)

theta_intercept <- c(-6.51, -3.64, -3.71, -3.93, -4.74, -3.43, -4.45, -2.40, -3.60, -5.85,-1.20,-3.50, -4.51, -2.39, -4.46, -3.53, -5.95,-4.00, -3.62,-2.14,-5.14, -2.65, -3.88,-2.91)

theta_covar <- c(0.05,0,0,0,0)




theta_test <- c(0.05,0,0,0,0)

n <- 50000
######SimulateDataPower will simulate data with 4 tumor characteristics
######theta_test is the effect of Gene
######theta_covar is the effect of one other covariate
######n is the sample size
temp.simu <- SimulateDataPower(theta_intercept,theta_test,theta_covar,n)
#######y.pheno.mis is the phenotype file, first column is case control
#######2:5 columns are er, pr, her2,grade
#######G is Gene
#######x_covar is one other covariate  
y.pheno.mis <- temp.simu[[1]]
G <- temp.simu[[2]]
x_covar <- temp.simu[[3]]
z.standard <- GenerateZstandard(y.pheno.mis,missingTumorIndicator = 888)
M <- nrow(z.standard)
#####################standard fixed effect two-stage model
model.fixed <- TwoStageModel(y.pheno.mis,
                             additive = cbind(G,x_covar),
                             missingTumorIndicator = 888)
#####################stage1to2.model
#####################we want a self design second stage function for G,
#####################for x_covar1, we want additive structure
#####################for x_covar2, we want saturated structure

################### a scractch function for design stage1to2
Stage1To2 <- function(y,
                      x.self.design,
                      z.design,
                      baselineonly = NULL,
                      additive = NULL,
                      pairwise.interaction = NULL,
                      saturated = NULL,
                      missingTumorIndicator = 888
){
  missing.data.vec <- GenerateMissingPosition(y,missingTumorIndicator)
  y.pheno.complete <- y[-missing.data.vec,]
  y <- y.pheno.complete
  tumor.number <- ncol(y)-1
  y.case.control <- y[,1]
  y.tumor <- y[,2:(tumor.number+1)]
  freq.subtypes <- GenerateFreqTable(y.pheno.complete)
  if(CheckControlTumor(y.case.control,y.tumor)==1){
    return(print("ERROR:The tumor characteristics for control subtypes should put as NA"))
  }
  tumor.names <- colnames(y.tumor)
  if(is.null(tumor.names)){
    tumor.names <- paste0(c(1:tumor.number))
  }
  tumor.character.cat = GenerateTumorCharacterCat(y.pheno.complete)
  z.design.baselineonly <- GenerateZDesignBaselineonly(tumor.character.cat,
                                                       tumor.number,
                                                       tumor.names,
                                                       freq.subtypes)
  z.design.additive <- GenerateZDesignAdditive(tumor.character.cat,
                                               tumor.number,
                                               tumor.names,
                                               freq.subtypes)
  z.design.pairwise.interaction <- GenerateZDesignPairwiseInteraction(tumor.character.cat,
                                                                      tumor.number,
                                                                      tumor.names,
                                                                      freq.subtypes)
  z.design.saturated <- GenerateZDesignSaturated(tumor.character.cat,
                                                 tumor.number,
                                                 tumor.names,
                                                 freq.subtypes)
  full.second.stage.names <- colnames(z.design.saturated)
  covar.names <- GenerateSelfCovarName(x.self.design,
                                       baselineonly,
                                       additive,
                                       pairwise.interaction,
                                       saturated)
  z.all <- ZSelfDesigntoZall(x.self.design,
                             baselineonly,
                             additive,
                             pairwise.interaction,
                             saturated,
                             z.design,
                             z.design.baselineonly,
                             z.design.additive,
                             z.design.pairwise.interaction,
                             z.design.saturated)
  return(z.all)
}

z.standard <- GenerateZstandard(y.pheno.mis)
M <- nrow(z.standard)
x_covar1 <- x_covar
x_covar2 <- x_covar+rnorm(length(x_covar1))
####################all the ER positive as one group
####################all the ER negative as one group 
z.design.G <- matrix(0,M,2)
idx.1 <- which(z.standard[,1]==0)
z.design.G[idx.1,1] <- 1
idx.2 <- which(z.standard[,1]==1)
z.design.G[idx.1,2] <- 1
stage1to2.model <- Stage1To2(y.pheno.mis,
                          x.self.design = G,
                          z.design = z.design.G,
                          additive = x_covar1,
                          saturated = x_covar2)













##########Fixed effect global test for association################ 
##########Get the support function for fixed effect
#########Fixed effect supports, model fits under global null, G has no effect
score.test.support <- ScoreTestSupport(
  y.pheno.mis,
  baselineonly = NULL,
  additive = as.matrix(x_covar),
  pairwise.interaction = NULL,
  saturated = NULL,
  missingTumorIndicator = 888
)
#######Calculate the score for fixed effect
score.test = ScoreTest(y=y.pheno.mis,
                       x= G,
                       second.stage.structure = "additive",
                       score.test.support=score.test.support,
                       missingTumorIndicator=888)
######Get the score and information for fixed effect
score = score.test[[1]]
infor = score.test[[2]]
######Get p value for the mixed effect global test for association
p_global = DisplayFixedScoreTestResult(score,infor)

##########Mixed effect global test for association################ 
##########z.design fixed put the baseline effect and 
z.design.fixed <- cbind(rep(1,M),z.standard[,1])
z.design.random <-z.standard[,2:ncol(z.standard)]
#########Fixed effect supports, model fits under global null, G has no effect
score.test.support.fixed <- ScoreTestSupportMixedModel(
  y.pheno.mis,
  baselineonly = NULL,
  additive = as.matrix(x_covar),
  pairwise.interaction = NULL,
  saturated = NULL,
  missingTumorIndicator = 888,
  c(theta_intercept,
    theta_covar)
)


#######Calculate the score for fixed effect
score.test.fixed<- ScoreTestMixedModel(y=y.pheno.mis,
                                       x=as.matrix(G),
                                       z.design=z.design.fixed,
                                       score.test.support=score.test.support.fixed,
                                       missingTumorIndicator=888
)
######Get the score and information for fixed effect
score.fixed <- score.test.fixed[[1]]
infor.fixed <- score.test.fixed[[2]]

######Random effect supports, model fits under the null that only random effect is 0. G has the fixed effect
score.test.support.random <- ScoreTestSupportMixedModelSelfDesign(
  y.pheno.mis,
  x.self.design  = G,
  z.design = z.design.fixed,
  additive =  as.matrix(x_covar),
  pairwise.interaction = NULL,
  saturated = NULL,
  missingTumorIndicator = 888,
  delta0 = c(theta_intercept,theta_test[1:ncol(z.design.fixed)],
             theta_covar)
)
#####Calculate the score for random effect
score.test.random<- ScoreTestMixedModel(y=y.pheno.mis,
                                        x=as.matrix(G),
                                        z.design=z.design.random,
                                        score.test.support=score.test.support.random,
                                        missingTumorIndicator=888)
######Get the score and information for random effect
score.random <- score.test.random[[1]]
infor.random <- score.test.random[[2]]
######Get p value for the mixed effect global test for association
p_mglobal <- DisplayMixedScoreTestResult(score.fixed,
                                         infor.fixed,
                                         score.random,
                                         infor.random)[1]
