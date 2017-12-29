all_additive_support <- function(y.pheno.mis1,
                                 y.pheno.mis2,
                                 x.all.mis1,
                                 x.all.mis2,
                                 idx.known,
                                 conditional.round,
                                 all.idx){

  Heter.result.Icog = EMmvpoly(y.pheno.mis1,baselineonly = NULL,additive = x.all.mis1,pairwise.interaction = NULL,saturated = NULL,missingTumorIndicator = 888)
  z.standard <- Heter.result.Icog[[12]]
  M <- nrow(z.standard)
  number.of.tumor <- ncol(z.standard)
  Heter.result.Onco = EMmvpoly(y.pheno.mis2,baselineonly = NULL,additive = x.all.mis2,pairwise.interaction = NULL,saturated = NULL,missingTumorIndicator = 888)
  result.all <- NULL
  for(i in 1:length(idx.known)){
    log.odds.icog <- Heter.result.Icog[[1]][M+(1:(1+number.of.tumor))+(i-1)*(1+number.of.tumor)]
    sigma.log.odds.icog <- Heter.result.Icog[[2]][M+(1:(1+number.of.tumor))+(i-1)*(1+number.of.tumor),M+(1:(1+number.of.tumor))+(i-1)*(1+number.of.tumor)]
    log.odds.onco <- Heter.result.Onco[[1]][M+(1:(1+number.of.tumor))+(i-1)*(1+number.of.tumor)]
    sigma.log.odds.onco <- Heter.result.Onco[[2]][M+(1:(1+number.of.tumor))+(i-1)*(1+number.of.tumor),M+(1:(1+number.of.tumor))+(i-1)*(1+number.of.tumor)]
    meta.result <- LogoddsMetaAnalysis(log.odds.icog,
                                       sigma.log.odds.icog,
                                       log.odds.onco,
                                       sigma.log.odds.onco)
    
    second.stage.logodds.meta <- meta.result[[1]]
    second.stage.sigma.meta <- meta.result[[2]]
    
    
    test.result.second.wald <- DisplaySecondStageTestResult(second.stage.logodds.meta,second.stage.sigma.meta)
    snp.infor <- fine_mapping[idx.known[i],c(1,3,4,5,6)]
    known.flag <- idx.known[i]
    mark <- "known snp"
    conditional.p.value <- NA
    result <- data.frame(snp.infor,conditional.p.value,test.result.second.wald,known.flag,mark)
    colnames(result) <- c("rs id","CHR","position","Alleles","MAF","conditional p value","Baseline OR(95% CI)","Baseline P",
                          "ER OR(95% CI)","ER P",
                          "PR OR(95% CI)","PR P",
                          "HER2 OR(95% CI)","HER2 P",
                          "Grade OR(95% CI)","Grade P",
                          "Global Association P",
                          "Global Heterogeneity P",
                          "known flag",
                          "mark")
    result.all <- rbind(result.all,result)
  }
  
  
  for(i in 1:conditional.round){
    log.odds.icog <- Heter.result.Icog[[1]][M+(1:(1+number.of.tumor))+(i-1)*(1+number.of.tumor)+length(idx.known)*(1+number.of.tumor)]
    sigma.log.odds.icog <- Heter.result.Icog[[2]][M+(1:(1+number.of.tumor))+(i-1)*(1+number.of.tumor)+length(idx.known)*(1+number.of.tumor),M+(1:(1+number.of.tumor))+(i-1)*(1+number.of.tumor)+length(idx.known)*(1+number.of.tumor)]
    log.odds.onco <- Heter.result.Onco[[1]][M+(1:(1+number.of.tumor))+(i-1)*(1+number.of.tumor)+length(idx.known)*(1+number.of.tumor)]
    sigma.log.odds.onco <- Heter.result.Onco[[2]][M+(1:(1+number.of.tumor))+(i-1)*(1+number.of.tumor)+length(idx.known)*(1+number.of.tumor),M+(1:(1+number.of.tumor))+(i-1)*(1+number.of.tumor)+length(idx.known)*(1+number.of.tumor)]
    meta.result <- LogoddsMetaAnalysis(log.odds.icog,
                                       sigma.log.odds.icog,
                                       log.odds.onco,
                                       sigma.log.odds.onco)
    
    second.stage.logodds.meta <- meta.result[[1]]
    second.stage.sigma.meta <- meta.result[[2]]
    
    
    test.result.second.wald <- DisplaySecondStageTestResult(second.stage.logodds.meta,second.stage.sigma.meta)
    
    idx.temp <- all.idx[[i]]
    conditional.results.temp <- all.condition.results[[i]]
    
    SNP.ICOGS.temp <- conditional.results.temp[idx.temp,1] 
    SNP.ONCO.temp <-  conditional.results.temp[idx.temp,2]  
    
    
    snp.infor <- meta_result_shared_1p%>%filter(SNP.ICOGS==SNP.ICOGS.temp&
                                                  SNP.ONCO==SNP.ONCO.temp)
    snp.infor <- snp.infor[,c(1,11,3,2,4)]
    conditional.p.value <- conditional.results.temp[idx.temp,4]
    known.flag <- conditional.results.temp[idx.temp,3]
    mark <- paste0("condition analysis ",i," round")
    
    result <- data.frame(snp.infor,conditional.p.value,test.result.second.wald,known.flag,mark)
    colnames(result) <- c("rs id","CHR","position","Alleles","MAF","conditional p value","Baseline OR(95% CI)","Baseline P",
                          "ER OR(95% CI)","ER P",
                          "PR OR(95% CI)","PR P",
                          "HER2 OR(95% CI)","HER2 P",
                          "Grade OR(95% CI)","Grade P",
                          "Global Association P",
                          "Global Heterogeneity P",
                          "known flag",
                          "mark")
    result.all <- rbind(result.all,result)
  }
  return(result.all)
}




all_additive_support_onco <- function(y.pheno.mis2,
                                 x.all.mis2,
                                 idx.known,
                                 conditional.round,
                                 all.idx){

  Heter.result.Onco = EMmvpoly(y.pheno.mis2,baselineonly = NULL,additive = x.all.mis2,pairwise.interaction = NULL,saturated = NULL,missingTumorIndicator = 888)
  z.standard <- Heter.result.Icog[[12]]
  M <- nrow(z.standard)
  number.of.tumor <- ncol(z.standard)
  
  result.all <- NULL
  for(i in 1:length(idx.known)){
    log.odds.onco <- Heter.result.Onco[[1]][M+(1:(1+number.of.tumor))+(i-1)*(1+number.of.tumor)]
    sigma.log.odds.onco <- Heter.result.Onco[[2]][M+(1:(1+number.of.tumor))+(i-1)*(1+number.of.tumor),M+(1:(1+number.of.tumor))+(i-1)*(1+number.of.tumor)]
     
    
    
    test.result.second.wald <- DisplaySecondStageTestResult(og.odds.onco, sigma.log.odds.onco)
    snp.infor <- fine_mapping[idx.known[i],c(1,3,4,5,6)]
    known.flag <- idx.known[i]
    mark <- "known snp"
    conditional.p.value <- NA
    result <- data.frame(snp.infor,conditional.p.value,test.result.second.wald,known.flag,mark)
    colnames(result) <- c("rs id","CHR","position","Alleles","MAF","conditional p value","Baseline OR(95% CI)","Baseline P",
                          "ER OR(95% CI)","ER P",
                          "PR OR(95% CI)","PR P",
                          "HER2 OR(95% CI)","HER2 P",
                          "Grade OR(95% CI)","Grade P",
                          "Global Association P",
                          "Global Heterogeneity P",
                          "known flag",
                          "mark")
    result.all <- rbind(result.all,result)
  }
  
  
  for(i in 1:conditional.round){
    log.odds.onco <- Heter.result.Onco[[1]][M+(1:(1+number.of.tumor))+(i-1)*(1+number.of.tumor)+length(idx.known)*(1+number.of.tumor)]
    sigma.log.odds.onco <- Heter.result.Onco[[2]][M+(1:(1+number.of.tumor))+(i-1)*(1+number.of.tumor)+length(idx.known)*(1+number.of.tumor),M+(1:(1+number.of.tumor))+(i-1)*(1+number.of.tumor)+length(idx.known)*(1+number.of.tumor)]
    
    test.result.second.wald <- DisplaySecondStageTestResult(log.odds.onco,sigma.log.odds.onco)
    
    idx.temp <- all.idx[[i]]
    conditional.results.temp <- all.condition.results[[i]]
    
    
    SNP.ONCO.temp <-  conditional.results.temp[idx.temp,2]  
    
    
    snp.infor <- meta_result_shared_1p%>%filter(SNP.ONCO==SNP.ONCO.temp)
    snp.infor <- snp.infor[,c(1,11,3,2,4)]
    conditional.p.value <- conditional.results.temp[idx.temp,4]
    known.flag <- conditional.results.temp[idx.temp,3]
    mark <- paste0("condition analysis ",i," round")
    
    result <- data.frame(snp.infor,conditional.p.value,test.result.second.wald,known.flag,mark)
    colnames(result) <- c("rs id","CHR","position","Alleles","MAF","conditional p value","Baseline OR(95% CI)","Baseline P",
                          "ER OR(95% CI)","ER P",
                          "PR OR(95% CI)","PR P",
                          "HER2 OR(95% CI)","HER2 P",
                          "Grade OR(95% CI)","Grade P",
                          "Global Association P",
                          "Global Heterogeneity P",
                          "known flag",
                          "mark")
    result.all <- rbind(result.all,result)
  }
  return(result.all)
}












