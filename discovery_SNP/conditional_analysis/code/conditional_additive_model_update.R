condition_additive_model_update <- function(y.pheno.mis1,
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
                                     z.design.score.casecase.ER,
                                     score.test.support.icog = NULL,
                                     score.test.support.onco = NULL,
                                     conditional.snps.icog = NULL,
                                     conditional.snps.onco = NULL,
                                     region.all
){
  if(is.na(snp.name.icog)&is.na(snp.name.onco)){
    p.value <- 1  
    return(p.value)
  }else if((is.na(snp.name.icog)==F)&(is.na(snp.name.onco)==T)){
    if(known.flag==178|known.flag==207|known.flag==122|known.flag==172){
      p.value <- 1
      return(p.value)
    }else{
      idx.control <- which(y.pheno.mis1[,1]==0)
      idx.known <- which(region.all==region.all[known.flag])  
      
      known.snp.value.icog <- as.matrix(known.all.mis1[,idx.known])
      
      
      known.snp.value.icog.control <- known.snp.value.icog[idx.control]
      
      snp.icog.control <- snp.icog[idx.control]
      conditional.snps.icog <- as.matrix(conditional.snps.icog)
      conditional.snps.icog.control <- conditional.snps.icog[idx.control,]
      cor.icog.control <- cor(snp.icog.control,conditional.snps.icog.control)
      cor.max <- max(cor.icog.control)
      
      if(max(cor(snp.icog.control,known.snp.value.icog.control)^2)>=0.8|cor.max^2>=0.8){
        p.value <- 1
        return(p.value)
      }else{
        
        idx.known <- which(region.all==region.all[known.flag])  
        
        known.snp.value.icog <- as.matrix(known.all.mis1[,idx.known])
        
        x.all.mis1 <- cbind(snp.icog,known.snp.value.icog,
                            conditional.snps.icog,
                            x.covar.mis1)
        
        score.test.icog.baseline.ER<- ScoreTestSelfDesign(y=y.pheno.mis1,
                                                          x=x.all.mis1[,1,drop=F],
                                                          z.design=z.design.score.baseline.ER,
                                                          score.test.support=score.test.support.icog,
                                                          missingTumorIndicator=888)
        score.icog.baseline.ER <- score.test.icog.baseline.ER[[1]]
        infor.icog.baseline.ER <- score.test.icog.baseline.ER[[2]]
        
        
        
        
        
        score.test.support.icog.casecase <- ScoreTestSupportMixedModelSelfDesign(y=y.pheno.mis1,
                                                                                 x.self.design = x.all.mis1[,1,drop=F],
                                                                                 z.design = z.design.score.baseline.ER,
                                                                                 additive=x.all.mis1[,2:ncol(x.all.mis1)],
                                                                                 missingTumorIndicator = 888)
        
        score.test.icog.casecase.ER<- ScoreTestMixedModel(y=y.pheno.mis1,
                                                          x=x.all.mis1[,1,drop=F],
                                                          z.design = z.design.score.casecase.ER,
                                                          
                                                          score.test.support= score.test.support.icog.casecase,
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
    }
    
    
    
    
    
    
  }else if(((is.na(snp.name.icog)==T)&(is.na(snp.name.onco)==F))|(known.flag==178&is.na(snp.name.onco)==F)|(known.flag==207&is.na(snp.name.onco)==F)|(known.flag==172&is.na(snp.name.onco)==F)|(known.flag==122&is.na(snp.name.onco)==F)){
    
    idx.control <- which(y.pheno.mis2[,1]==0)
    idx.known <- which(region.all==region.all[known.flag])  
    known.snp.value.onco <- as.matrix(known.all.mis2[,idx.known])
    
    known.snp.value.onco.control <- known.snp.value.onco[idx.control]
    
   
    snp.onco.control <- snp.onco[idx.control]
    
    conditional.snps.onco <- as.matrix( conditional.snps.onco)
    conditional.snps.onco.control <- conditional.snps.onco[idx.control,]
    
    
    cor.onco.control <- cor(snp.onco.control,conditional.snps.onco.control)
    cor.max <- max(cor.onco.control)
    
    if(max(cor(snp.onco.control,known.snp.value.onco.control)^2)>=0.8|cor.max^2>=0.8){
      p.value <- 1
      return(p.value)
    }else{
      idx.known <- which(region.all==region.all[known.flag])  
      known.snp.value.onco <- as.matrix(known.all.mis2[,idx.known])
      x.all.mis1 <- cbind(snp.icog,known.snp.value.icog,
                          conditional.snps.icog,
                          x.covar.mis1)
      x.all.mis2 <- cbind(snp.onco,known.snp.value.onco,
                          conditional.snps.onco,
                          x.covar.mis2) 
      score.test.onco.baseline.ER<- ScoreTestSelfDesign(y=y.pheno.mis2,
                                                        x=x.all.mis2[,1,drop=F],
                                                        z.design=z.design.score.baseline.ER,
                                                        score.test.support=score.test.support.onco,
                                                        missingTumorIndicator=888)
      
      score.onco.baseline.ER <- score.test.onco.baseline.ER[[1]]
      infor.onco.baseline.ER <- score.test.onco.baseline.ER[[2]]
      
      
      score.test.support.onco.casecase <- ScoreTestSupportMixedModelSelfDesign(y=y.pheno.mis2,
                                                                               x.self.design = x.all.mis2[,1,drop=F],
                                                                               z.design = z.design.score.baseline.ER,
                                                                               additive=x.all.mis2[,2:ncol(x.all.mis2)],
                                                                               missingTumorIndicator = 888)
      
      score.test.onco.casecase.ER<- ScoreTestMixedModel(y=y.pheno.mis2,
                                                        x=x.all.mis2[,1,drop=F],
                                                        z.design = z.design.score.casecase.ER,
                                                        
                                                        score.test.support= score.test.support.onco.casecase,
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
    
    idx.known <- which(region.all==region.all[known.flag])  
    
    
    
    known.snp.value.onco <- as.matrix(known.all.mis2[,idx.known])
    
    known.snp.value.onco.control <- known.snp.value.onco[idx.control]
    snp.onco.control <- snp.onco[idx.control]
    

    
    conditional.snps.onco <- as.matrix( conditional.snps.onco)
    conditional.snps.onco.control <- conditional.snps.onco[idx.control,]
    
    
    cor.onco.control <- cor(snp.onco.control,conditional.snps.onco.control)
    cor.max <- max(cor.onco.control)
    
    
    
    if(max(cor(snp.onco.control,known.snp.value.onco.control)^2)>=0.8|cor.max^2>=0.8){
      p.value <- 1
      return(p.value)
    }else{
      idx.known <- which(region.all==region.all[known.flag])  
      
      known.snp.value.icog <- as.matrix(known.all.mis1[,idx.known])
      
      known.snp.value.onco <- as.matrix(known.all.mis2[,idx.known])
      x.all.mis1 <- cbind(snp.icog,known.snp.value.icog,
                          conditional.snps.icog,
                          x.covar.mis1)
      x.all.mis2 <- cbind(snp.onco,known.snp.value.onco,
                          conditional.snps.onco,
                          x.covar.mis2)
      score.test.icog.baseline.ER<- ScoreTestSelfDesign(y=y.pheno.mis1,
                                                        x=x.all.mis1[,1,drop=F],
                                                        z.design=z.design.score.baseline.ER,
                                                        score.test.support=score.test.support.icog,
                                                        missingTumorIndicator=888)
      score.icog.baseline.ER <- score.test.icog.baseline.ER[[1]]
      infor.icog.baseline.ER <- score.test.icog.baseline.ER[[2]]
      
      score.test.support.icog.casecase <- ScoreTestSupportMixedModelSelfDesign(y=y.pheno.mis1,
                                                                               x.self.design = x.all.mis1[,1,drop=F],
                                                                               z.design = z.design.score.baseline.ER,
                                                                               additive=x.all.mis1[,2:ncol(x.all.mis1)],
                                                                               missingTumorIndicator = 888)
      
      score.test.icog.casecase.ER<- ScoreTestMixedModel(y=y.pheno.mis1,
                                                        x=x.all.mis1[,1,drop=F],
                                                        z.design = z.design.score.casecase.ER,
                                                        
                                                        score.test.support= score.test.support.icog.casecase,
                                                        missingTumorIndicator=888)
      
      score.icog.casecase.ER <- score.test.icog.casecase.ER[[1]]
      infor.icog.casecase.ER <- score.test.icog.casecase.ER[[2]]
      
      
      score.test.onco.baseline.ER<- ScoreTestSelfDesign(y=y.pheno.mis2,
                                                        x=x.all.mis2[,1,drop=F],
                                                        z.design=z.design.score.baseline.ER,
                                                        score.test.support=score.test.support.onco,
                                                        missingTumorIndicator=888)
      
      score.onco.baseline.ER <- score.test.onco.baseline.ER[[1]]
      infor.onco.baseline.ER <- score.test.onco.baseline.ER[[2]]
      
      score.test.support.onco.casecase <- ScoreTestSupportMixedModelSelfDesign(y=y.pheno.mis2,
                                                                               x.self.design = x.all.mis2[,1,drop=F],
                                                                               z.design = z.design.score.baseline.ER,
                                                                               additive=x.all.mis2[,2:ncol(x.all.mis2)],
                                                                               missingTumorIndicator = 888)
      
      score.test.onco.casecase.ER<- ScoreTestMixedModel(y=y.pheno.mis2,
                                                        x=x.all.mis2[,1,drop=F],
                                                        z.design = z.design.score.casecase.ER,
                                                        
                                                        score.test.support= score.test.support.onco.casecase,
                                                        missingTumorIndicator=888)
      
      
      score.onco.casecase.ER <- score.test.onco.casecase.ER[[1]]
      infor.onco.casecase.ER <- score.test.onco.casecase.ER[[2]]
      
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








