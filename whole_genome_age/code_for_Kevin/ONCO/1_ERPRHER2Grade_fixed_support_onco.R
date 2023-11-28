
library(devtools)
library(bc2, 
        lib.loc ="/home/zhangh24/R/x86_64-pc-linux-gnu-library/4.2/")

library(data.table)
setwd("/data/zhangh24/breast_cancer_data_analysis/")
data2 <- fread("./data/Onco_euro_v10_10232017.csv",header=T)
data2 <- as.data.frame(data2)
y.pheno.mis2 <- cbind(data2$Behaviour1,data2$ER_status1,data2$PR_status1,data2$HER2_status1,data2$Grade1)
#y.pheno.mis2 <- cbind(data2$Behaviour1,data2$PR_status1,data2$ER_status1,data2$HER2_status1)
colnames(y.pheno.mis2) = c("Behaviour","ER",
                           "PR","HER2","Grade")

x.covar.mis2 <- data2[,c(5:14,204)]
ages <- data2[,204]
idx.complete <- which(ages!=888)

y.pheno.mis2 <- y.pheno.mis2[idx.complete,]
x.covar.mis2 <- x.covar.mis2[idx.complete,]




score.test.support.onco.ERPRHER2Grade <- ScoreTestSupport(
  y.pheno.mis2,
  baselineonly = NULL,
  additive = x.covar.mis2,
  pairwise.interaction = NULL,
  saturated = NULL,
  missingTumorIndicator = 888
)

save(score.test.support.onco.ERPRHER2Grade,file="./whole_genome_age/ONCO/ERPRHER2GRADE_fixed_baseline/result/score.test.support.onco.ERPRHER2Grade.Rdata")



y = y.pheno.mis2;
baselineonly = NULL;
additive = x.covar.mis2;
pairwise.interaction = NULL;
saturated = NULL;
missingTumorIndicator = 888
y <- as.matrix(y)
tumor.number <- ncol(y)-1
y.case.control <- y[,1]
y.tumor <- y[,2:(tumor.number+1)]
y.pheno.complete <- GenerateCompleteYPheno(y,missingTumorIndicator)
freq.subtypes <- GenerateFreqTable(y.pheno.complete)
tumor.names <- colnames(y.tumor)
if(is.null(tumor.names)){
  tumor.names <- paste0(c(1:tumor.number))
}
tumor.character.cat = GenerateTumorCharacterCat(y.pheno.complete)
cutoff = 10
z.design.baselineonly <- GenerateZDesignBaselineonly(tumor.character.cat,
                                                     tumor.number,
                                                     tumor.names,
                                                     freq.subtypes,
                                                     cutoff)
z.design.additive <- GenerateZDesignAdditive(tumor.character.cat,
                                             tumor.number,
                                             tumor.names,
                                             freq.subtypes,
                                             cutoff)
z.design.pairwise.interaction <- GenerateZDesignPairwiseInteraction(tumor.character.cat,
                                                                    tumor.number,
                                                                    tumor.names,
                                                                    freq.subtypes,
                                                                    cutoff)
z.design.saturated <- GenerateZDesignSaturated(tumor.character.cat,
                                               tumor.number,
                                               tumor.names,
                                               freq.subtypes,
                                               cutoff)
z.all <- ZDesigntoZall(baselineonly,
                       additive,
                       pairwise.interaction,
                       saturated,
                       z.design.baselineonly,
                       z.design.additive,
                       z.design.pairwise.interaction,
                       z.design.saturated)
delta0 <-StartValueFunction(freq.subtypes,y.case.control,z.all)
x.all <- GenerateXAll(y,baselineonly,additive,pairwise.interaction,saturated)
z.standard <- z.design.additive[,-1]
tol <- as.numeric(1e-04)
tolMaxstep <- as.numeric(1e-03)
#delta_old <- rep(0,length(delta0))
delta_old <- delta0
##EM algorithm
##first E step
#print(paste0("Begin EM algorithm"))
#print(paste0("EM round: 1"))
prob.fit.result <- ProbFitting(delta_old,as.matrix(y),x.all,
                               z.standard,z.all,missingTumorIndicator)
y_em <- prob.fit.result[[1]]
missing.vec <- as.numeric(as.vector(prob.fit.result[[2]]))
missing.mat <- prob.fit.result[[3]]
missing.mat.vec <- as.numeric(as.vector(missing.mat))
missing.number <- as.integer(length(missing.vec))
idx.drop = prob.fit.result[[4]]
if(length(idx.drop)!=0){
  x.all <- x.all[-idx.drop,,drop=F]
  y_em <- y_em[-idx.drop,,drop=F]
}
for(k in 1:length(missing.vec)){
  missing.vec[k] <- missing.vec[k]-sum(missing.vec[k]>=idx.drop)
}
N <- as.integer(nrow(x.all))
#x <- cbind(1,x)
#p <- ncol(x)
M <- as.integer(nrow(z.standard))

NCOV   <- as.integer(ncol(x.all))
NM     <- N*M
nparm  <- as.integer(length(delta0))
deltai <- as.numeric(delta0)

NITER  <- as.integer(500)
Y <- as.numeric(as.vector(y_em))
X <- as.numeric(as.vector(x.all))
ZallVec = as.numeric(as.vector(z.all))
Znr = as.integer(nrow(z.all))
Znc = as.integer(ncol(z.all))
debug     <- as.integer(1)
ret_rc    <- as.integer(1)
ret_delta <- as.numeric(rep(-9999, nparm))
ret_info <- as.numeric(rep(-9999,nparm^2))
ret_p <- as.numeric(rep(0,NM))
ret_Inv_info_vec <- as.numeric(as.vector(matrix(0,Znc,Znc)))
YminusP <- Y
W_obs <- as.numeric(rep(0,N*M*M))
WXZ_vec <- as.numeric(rep(0,N*M*Znc))
WX_vec <- as.numeric(rep(0,N*M*Znr))
temp <- .C("EMStepScoreSupport",
           deltai,
           nparm,
           Y=Y,
           X,
           ZallVec,
           Znr,Znc, N, M, NCOV, NITER,
           tol,
           tolMaxstep,
           debug,
           ret_rc=ret_rc,
           ret_delta=ret_delta,
           ret_info=ret_info,
           ret_p=ret_p,
           missing.vec,
           missing.mat.vec,missing.number,
           ret_Inv_info_vec=ret_Inv_info_vec,
           YminusP=YminusP,
           W_obs = W_obs,
           WXZ_vec = WXZ_vec,
           WX_vec = WX_vec)
print(paste0("EM Algorithm Converged"))
inv_info_vec <- temp$ret_Inv_info_vec
YminusP <- temp$YminusP
W_obs <- temp$W_obs
WXZ_vec <- temp$WXZ_vec
WX_vec <- temp$WX_vec
