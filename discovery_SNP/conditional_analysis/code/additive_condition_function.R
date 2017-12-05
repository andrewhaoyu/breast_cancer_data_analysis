condition_additive_model <- function(y.pheno.mis1,
                                     x.covar.mis1,
                                     snp.name.icog,
                                     snp.icog,
                                     y.pheno.mis2,
                                     x.covar.mis2,
                                     snp.name.onco,
                                     snp.onco,
                                     known.flag,
                                     known.all.mis1,
                                     known.all.mis2,
                                     z.standard,
                                     z.additive.design,
                                     M,
                                     number.of.tumor,
                                     z.design.score.baseline,
                                     z.design.score.casecase,
                                     z.design.score.baseline.ER,
                                     z.design.score.casecase.ER){
  if(is.na(snp.name.icog)&is.na(snp.name.onco)){
    p.value <- 1  
    return(p.value)
  }else if((is.na(snp.name.icog)==F)&(is.na(snp.name.onco)==T)){
    idx.control <- which(y.pheno.mis1[,1]==0)
    known.snp.value.icog <- known.all.mis1[,known.flag]
    known.snp.value.icog.control <- known.snp.value.icog[idx.control]
    
    snp.icog.control <- snp.icog[idx.control]
    
    if(cor(snp.icog.control,known.snp.value.icog.control)^2>=0.8){
      p.value <- 1
      return(p.value)
    }else{
      known.snp.value.icog <- known.all.mis1[,known.flag]
      known.snp.value.onco <- known.all.mis2[,known.flag]
      x.all.mis1 <- cbind(snp.icog,known.snp.value.icog,
                          x.covar.mis1)
      x.all.mis2 <- cbind(snp.onco,known.snp.value.onco,
                          x.covar.mis2)
      
      score.test.support.icog <- ScoreTestSupport(
        y.pheno.mis1,
        baselineonly = NULL,
        additive = x.all.mis1[,2:ncol(x.all.mis1)],
        pairwise.interaction = NULL,
        saturated = NULL,
        missingTumorIndicator = 888
      )
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
      
      mixed.baseline.er <- DisplayMixedScoreTestResult(score.icog.baseline.ER,
                                                       infor.icog.baseline.ER,
                                                       score.icog.casecase.ER,
                                                       infor.icog.casecase.ER)  
      p.value <- mixed.baseline.er[1]
      return(p.value)
    }
    
    
    
    
    
    
  }else if(((is.na(snp.name.icog)==T)&(is.na(snp.name.onco)==F))|known.flag==178|known.flag==207){
    
    idx.control <- which(y.pheno.mis2[,1]==0)
    known.snp.value.onco <- known.all.mis2[,known.flag]
    known.snp.value.onco.control <- known.snp.value.onco[idx.control]
    
    snp.onco.control <- snp.onco[idx.control]
    if(cor(snp.onco.control,known.snp.value.onco.control)^2>=0.8){
      p.value <- 1
      return(p.value)
    }else{
      
      known.snp.value.onco <- known.all.mis2[,known.flag]
      x.all.mis2 <- cbind(snp.onco,known.snp.value.onco,
                          x.covar.mis2)
      
      
      
      
      
      
      score.test.support.onco <- ScoreTestSupport(
        y.pheno.mis2,
        baselineonly = NULL,
        additive = x.all.mis2[,2:ncol(x.all.mis2)],
        pairwise.interaction = NULL,
        saturated = NULL,
        missingTumorIndicator = 888
      )
      
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
      
      test.result.second.mixed.ER <- DisplayMixedScoreTestResult(score.onco.baseline.ER,
                                                                 infor.onco.baseline.ER,
                                                                 score.onco.casecase.ER,
                                                                 infor.onco.casecase.ER)  
      p.value <- as.numeric(test.result.second.mixed.ER[1])
      return(p.value)
      
    }
    
    
    
    
    
  }else{
    idx.control <- which(y.pheno.mis2[,1]==0)
    known.snp.value.onco <- known.all.mis2[,known.flag]
    known.snp.value.onco.control <- known.snp.value.onco[idx.control]
    snp.onco.control <- snp.onco[idx.control]
    
    if(cor(snp.onco.control,known.snp.value.onco.control)^2>=0.8){
      p.value <- 1
      return(p.value)
    }else{
      known.snp.value.icog <- known.all.mis1[,known.flag]
      known.snp.value.onco <- known.all.mis2[,known.flag]
      x.all.mis1 <- cbind(snp.icog,known.snp.value.icog,
                          x.covar.mis1)
      x.all.mis2 <- cbind(snp.onco,known.snp.value.onco,
                          x.covar.mis2)
      score.test.support.icog <- ScoreTestSupport(
        y.pheno.mis1,
        baselineonly = NULL,
        additive = x.all.mis1[,2:ncol(x.all.mis1)],
        pairwise.interaction = NULL,
        saturated = NULL,
        missingTumorIndicator = 888
      )
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
      
      
      
      
      score.test.support.onco <- ScoreTestSupport(
        y.pheno.mis2,
        baselineonly = NULL,
        additive = x.all.mis2[,2:ncol(x.all.mis2)],
        pairwise.interaction = NULL,
        saturated = NULL,
        missingTumorIndicator = 888
      )
      z.design.score.baseline <- matrix(rep(1,M),ncol=1)
      z.design.score.casecase <-z.standard
      z.design.score.baseline.ER <- cbind(z.design.score.baseline,z.standard[,1])
      z.design.score.casecase.ER <- z.standard[,2:ncol(z.standard)]
      
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
      
      meta.result <- LogoddsMetaAnalysis(log.odds.icog,
                                         sigma.log.odds.icog,
                                         log.odds.onco,
                                         sigma.log.odds.onco)
      
      second.stage.logodds.meta <- meta.result[[1]]
      second.stage.sigma.meta <- meta.result[[2]]
      test.result.second.wald <- DisplaySecondStageTestResult(second.stage.logodds.meta,second.stage.sigma.meta)
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
      test.result.second.mixed.ER <- DisplayMixedScoreTestResult(score.meta.baseline.ER,
                                                                 infor.meta.baseline.ER,
                                                                 score.meta.casecase.ER,
                                                                 infor.meta.casecase.ER)  
      test.result.second.mixed.ER <- data.frame(t(test.result.second.mixed.ER))
      
      p.value <- as.numeric(test.result.second.mixed.ER[1])
      return(p.value)
      
      
    }
    
    
    
    
    
    
    
    
    
    
    
    
  }
  
}








