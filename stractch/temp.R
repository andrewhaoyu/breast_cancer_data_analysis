# OneStepMLE <- function(y,
#                        baselineonly,
#                        additive,
#                        pairwise.interaction,
#                        saturated,
#                        missingTumorIndicator){
#                        
library(devtools)
#install_github("andrewhaoyu/bc2", ref='development',args = c('--library="/home/zhangh24/R/x86_64-pc-linux-gnu-library/3.4"'))
i1 = 1

install_github("andrewhaoyu/bc2",ref='development')
library(bc2)
library(data.table)
setwd("/data/zhangh24/breast_cancer_data_analysis/")
data1 <- fread("./data/iCOGS_euro_v10_10232017.csv",header=T)
data1 <- as.data.frame(data1)
y.pheno.mis1 <- cbind(data1$Behaviour1,data1$ER_status1,data1$PR_status1,data1$HER2_status1,data1$Grade1)
colnames(y.pheno.mis1) = c("Behavior","ER","PR","HER2","Grade")

# Grade1.fake <- data1$Grade1
# Grade1.fake[data1$Grade1==2|data1$Grade1==3] <- 1
# Grade1.fake[data1$Grade1==1] <- 0
#y.pheno.mis1 <- cbind(data1$Behaviour1,data1$PR_status1,data1$ER_status1,data1$HER2_status1,Grade1.fake)
# y.pheno.mis1 <- cbind(data1$Behaviour1,data1$PR_status1,data1$ER_status1,data1$HER2_status1)

x.test.all.mis1 <- data1[,c(27:203)]
###pc1-10 and age
x.covar.mis1 <- data1[,c(5:14)]



x.all.mis1 <- as.matrix(cbind(x.test.all.mis1[,i1],x.covar.mis1))


y <- y.pheno.mis1
baselineonly=NULL
additive=x.all.mis1
pairwise.interaction=NULL
saturated=NULL
missingTumorIndicator = 888

missing.data.vec <- GenerateMissingPosition(y,missingTumorIndicator)
y.pheno.complete <- y[-missing.data.vec,]
x.all.complete <- x.all[-missing.data.vec,]
initial.set <- InitialSetup(y.pheno.complete,
                            baselineonly,
                            additive,
                            pairwise.interaction,
                            saturated
)
###z standard matrix means the additive model z design matrix without baseline effect
###z standard matrix is used to match the missing tumor characteristics to the complete subtypes

delta0 = initial.set$delta0
z.all = initial.set$z.all
z.standard = initial.set$z.standard
z.deisign.baselineonly = initial.set$z.design.baseline.only
z.design.additive = initial.set$z.design.additive
z.design.pairwise.interaction = initial.set$z.design.pairwise.interaction
z.design.saturated = initial.set$z.design.saturated
x.all <- as.matrix(GenerateXAll(y,baselineonly,additive,pairwise.interaction,saturated))
covar.names <- initial.set$covar.names
x.all.complete <- x.all[-missing.data.vec,]


prob.fit.result <- ProbFitting(delta0,y.pheno.complete,x.all.complete,z.standard,z.all,missingTumorIndicator=NULL)
y.fit.complete <- prob.fit.result[[1]]
M <- as.integer(nrow(z.standard))
p.main <- ncol(z.standard)+1

tol <- as.numeric(1e-04)


delta_old <- delta0

N <- as.integer(nrow(x.all.complete))

M <- as.integer(nrow(z.standard))

NCOV   <- as.integer(ncol(x.all.complete))
NM     <- N*M
nparm  <- as.integer(length(delta0))
deltai <- as.numeric(delta0)

NITER  <- as.integer(500)
Y <- as.numeric(as.vector(y.fit.complete))
X <- as.numeric(as.vector(x.all.complete))
ZallVec = as.numeric(as.vector(z.all))
Znr = as.integer(nrow(z.all))
Znc = as.integer(ncol(z.all))
debug     <- as.integer(1)
ret_rc    <- as.integer(1)
ret_delta <- as.numeric(rep(-9999, nparm))
ret_info <- as.numeric(rep(-9999,nparm^2))
ret_p <- as.numeric(rep(0,NM))
ret_lxx <- as.numeric(rep(0,NM))
loglikelihood <- as.numeric(-1);



temp <- .C("Mvpoly_complete",deltai, nparm, Y=Y, X, ZallVec,Znr,Znc, N, M, NCOV, NITER, tol,
           debug, ret_rc=ret_rc, ret_delta=ret_delta,ret_info=ret_info,ret_p=ret_p,loglikelihood = loglikelihood)


delta0 <- temp$ret_delta


prob.fit.result <- ProbFitting(delta0,y,x.all,z.standard,z.all,missingTumorIndicator)
y.fit <- prob.fit.result[[1]]

missing.vec <- as.numeric(as.vector(prob.fit.result[[2]]))
missing.mat <- prob.fit.result[[3]]
missing.mat.vec <- as.numeric(as.vector(missing.mat))
missing.number <- as.integer(length(missing.vec))
complete.vec <- prob.fit.result[[4]]

N <- as.integer(nrow(x.all))
NM     <- N*M

deltai <- as.numeric(delta0)



Y <- as.numeric(as.vector(y.fit))
X <- as.numeric(as.vector(x.all))
ret_p <- as.numeric(rep(0,NM))
ret_lxx <- as.numeric(rep(0,NM))



temp <- .C("OneStepMLE",deltai, nparm, Y=Y, X, ZallVec,Znr,Znc, N, M, NCOV, NITER, tol,
           debug, ret_rc=ret_rc, ret_delta=ret_delta,ret_info=ret_info,ret_p=ret_p,missing.vec,
           missing.mat.vec,missing.number,loglikelihood = loglikelihood)


info <- matrix(unlist(temp$ret_info),nparm,nparm)
result <- list(temp$ret_delta,info,
               temp$ret_p)
y_em <- matrix(unlist(temp$Y),N,M)

# infor_mis_c <- infor_mis(y_em,x.all,z.all)
#infor_obs <- result[[2]]-infor_mis_c
delta=result[[1]]
infor_obs=result[[2]]
p=result[[3]]
loglikelihood = temp$loglikelihood
AIC = 2*nparm - 2*loglikelihood

OneStep.result <- list(delta=delta,
                       infor_obs=infor_obs,
                       p=p,y_em=y_em,
                       M=M,
                       NumberofTumor=ncol(z.standard),
                       loglikelihood = loglikelihood,
                       AIC = AIC
)

covariance.delta <- solve(OneStep.result$infor_obs)
loglikelihood <- OneStep.result$loglikelihood
AIC <- OneStep.result$AIC
second.stage.mat <-
  GenerateSecondStageMat(baselineonly,
                         additive,
                         pairwise.interaction,
                         saturated,
                         M,
                         full.second.stage.names,
                         covar.names,
                         delta,
                         z.design.additive,
                         z.design.pairwise.interaction,
                         z.design.saturated)
##take out the intercept from second stage parameters

takeout.intercept.result <- TakeoutIntercept(delta,covariance.delta,
                                             M,
                                             tumor.names,
                                             z.all,covar.names)
beta <- takeout.intercept.result$beta
covariance.beta <- takeout.intercept.result$covariance.beta
delta.no.inter <- takeout.intercept.result$delta.no.inter
covariance.delta.no.inter <-
  takeout.intercept.result$covariance.delta.no.inter
beta.no.inter <- takeout.intercept.result$beta.no.inter
covariance.beta.no.inter <- takeout.intercept.result$covariance.beta.no.inter



second.stage.test <- SecondStageTest(delta.no.inter,covariance.delta.no.inter,M,second.stage.mat)
global.test <- GenerateGlobalTest(delta.no.inter,
                                  covariance.delta.no.inter,
                                  M,
                                  second.stage.mat)
##beta represent first stage parameters

subtypes.names <- GenerateSubtypesName(z.design.additive,M,
                                       tumor.names)
first.stage.mat <- GenerateFirstStageMat(beta,
                                         covar.names,
                                         subtypes.names)

first.stage.test <- FirstStageTest(beta.no.inter,
                                   covariance.beta.no.inter,
                                   M,
                                   first.stage.mat)

















