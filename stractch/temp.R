# OneStepMLE <- function(y,
#                        baselineonly,
#                        additive,
#                        pairwise.interaction,
#                        saturated,
#                        missingTumorIndicator){
#                        
i1 = 1
library(devtools)
install_github("andrewhaoyu/bc2",ref='development')
library(bc2)
library(data.table)
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
x.covar.mis1 <- data1[,c(5:14,204)]

age <- data1[,204]
idx.complete <- which(age!=888)
x.all.mis1 <- as.matrix(cbind(x.test.all.mis1[,i1],x.covar.mis1))
colnames(x.all.mis1)[1] <- "gene"
y.pheno.mis1 <- y.pheno.mis1[idx.complete,]
x.all.mis1 <- x.all.mis1[idx.complete,]















missingTumorIndicator <- 888
y.pheno.complete <- GenerateCompleteYPheno(y.pheno.mis1,missingTumorIndicator)
x.all.complete1 <- GenerateCompleteXCovariates(y.pheno.mis1,x.all.mis1,missingTumorIndicator)
  
  y <- y.pheno.complete
  baselineonly <- NULL
  additive <- x.all.complete1
  pairwise.interaction <- NULL
  saturated <- NULL
  
  
  
  tumor.number <- ncol(y)-1
  y.case.control <- y[,1]
  y.tumor <- y[,2:(tumor.number+1)]
  y.pheno.complete <- GenerateCompleteYPheno(y,missingTumorIndicator)
  
  y <- as.matrix(y)
  tumor.number <- ncol(y)-1
  y.case.control <- y[,1]
  y.tumor <- y[,2:(tumor.number+1)]
  y.pheno.complete <- y
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
  z.all <- ZDesigntoZall(baselineonly,
                         additive,
                         pairwise.interaction,
                         saturated,
                         z.design.baselineonly,
                         z.design.additive,
                         z.design.pairwise.interaction,
                         z.design.saturated)
  delta0 <-StartValueFunction(freq.subtypes,y.case.control,z.all)
  #x.all has no intercept yet
  #we will add the intercept in C code
  x.all <- GenerateXAll(y,baselineonly,additive,pairwise.interaction,saturated)
  ###z standard matrix means the additive model z design matrix without baseline effect
  ###z standard matrix is used to match the missing tumor characteristics to the complete subtypes
  
  z.standard <- z.design.additive[,-1]
  
  tol <- as.numeric(1e-04)
  
  #delta_old <- rep(0,length(delta0))
  delta_old <- delta0
  ##EM algorithm
  ##first E step
  #print(paste0("Begin EM algorithm"))
  #print(paste0("EM round: 1"))
  y.fit <- ProbFitting(delta0,y,x.all,z.standard,z.all,missingTumorIndicator=NULL)[[1]]
  
  # sof <- "try5.so"
  # dyn.load(sof)
  
  N <- as.integer(nrow(x.all))
  #x <- cbind(1,x)
  #p <- ncol(x)
  M <- as.integer(nrow(z.standard))
  
  NCOV   <- as.integer(ncol(x.all))
  NM     <- N*M
  nparm  <- as.integer(length(delta0))
  deltai <- as.numeric(delta0)
  
  NITER  <- as.integer(500)
  Y <- as.numeric(as.vector(y.fit))
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
  
  
  temp <- .C("CompleteCasesScoreTestSupport",
             deltai,
             nparm,
             Y=Y,
             X,
             ZallVec,
             Znr,Znc, N, M, NCOV, NITER,
             tol,
             debug,
             ret_rc=ret_rc,
             ret_delta=ret_delta,
             ret_info=ret_info,
             ret_p=ret_p,
             ret_Inv_info_vec=ret_Inv_info_vec,
             YminusP=YminusP,
             W_obs = W_obs,
             WXZ_vec = WXZ_vec,
             WX_vec = WX_vec)
  
  inv_info_vec <- temp$ret_Inv_info_vec
  YminusP <- temp$YminusP
  W_obs <- temp$W_obs
  WXZ_vec <- temp$WXZ_vec
  WX_vec <- temp$WX_vec
  
  
  result <- list(inv_info_vec=inv_info_vec,YminusP=YminusP,W_obs=W_obs,WXZ_vec = WXZ_vec,zc=z.all)
  
  
  # score_support_result <- score_support(pxx,x.all,baselineonly,z.all,z.standard,y_em)
  #score_test_mis <- score_test_mis(y_em,baselineonly,score_support_result)
  #return(list(score_c=score_test_mis$score_c,infor_c = score_test_mis$infor_c))
  
  result[[6]] <- z.design.baselineonly
  result[[7]] <- z.design.additive
  result[[8]] <- z.design.pairwise.interaction
  result[[9]] <- z.design.saturated
  result[[10]] <- z.standard
























