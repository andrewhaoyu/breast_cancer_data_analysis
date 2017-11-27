#' Title
#'
#' @param baselineonly
#' @param additive
#' @param pairwise.interaction
#' @param saturated
#' @param M
#' @param full.second.stage.names
#' @param covar.names
#'
#' @return
#' @export
#'
#' @examples
GenerateSecondStageMat <- function(z.all,
                                   delta,
                                   baselineonly,
                                   additive,
                                   pairwise.interaction,
                                   saturated,
                                   M,
                                   full.second.stage.names,
                                   covar.names,
                                   z.design.additive,
                                   z.design.pairwise.interaction,
                                   z.design.saturated){
  baselineonly.number <- CountCovarNumber(baselineonly)
  additive.number <- CountCovarNumber(additive)
  pairwise.interaction.number <- CountCovarNumber(pairwise.interaction)
  saturated.number <- CountCovarNumber(saturated)
  ###second.stage.category for different model structures
  baselineonly.second.cat <- 1
  additive.second.cat <- ncol(z.design.additive)
  pairwise.interaction.second.cat <- ncol(z.design.pairwise.interaction)
  saturated.second.cat <- ncol(z.design.saturated)
  ###1 for intercept
  total.covar.number <- 1+ baselineonly.number+additive.number+
    pairwise.interaction.number+saturated.number
  
  max.number.second.stage.parameter <- 0
  if(baselineonly.number!=0){
    max.number.second.stage.parameter <- baselineonly.second.cat
  }
  if(additive.number!=0){
    max.number.second.stage.parameter <- additive.second.cat
  }
  if(pairwise.interaction.number!=0){
    max.number.second.stage.parameter <- pairwise.interaction.second.cat
  }
  if(saturated.number!=0){
    max.number.second.stage.parameter <- saturated.second.cat
  }
  
  
  second.stage.mat <- matrix(NA,max.number.second.stage.parameter,(total.covar.number-1))
  
  delta.no.inter <- delta[(M+1):length(delta)]
  ind.delta <- 0
  ind.covar <- 0
  if(baselineonly.number!=0){
    baselineonly.second.stage.delta <- matrix(NA,max.number.second.stage.parameter,baselineonly.number)
    baselineonly.second.stage.delta[1,(ind.covar+1):(ind.covar+baselineonly.number)] <- delta.no.inter[(ind.covar+1):(ind.covar+baselineonly.number)]
    second.stage.mat[,(ind.covar+1):(ind.covar+baselineonly.number)]<- baselineonly.second.stage.delta
    ind.covar <- ind.covar+baselineonly.number
    ind.delta <- ind.delta + ind.covar*baselineonly.number
  }
  
  if(additive.number!=0){
    
    for(i in 1:additive.number){
      ind.covar <- ind.covar+1
      second.stage.mat[1:additive.second.cat,ind.covar]<-
        delta.no.inter[(ind.delta+1):(ind.delta+additive.second.cat)]
      ind.delta <- ind.delta + additive.second.cat
    }
  }
  
  
  if(pairwise.interaction.number!=0){
    
    for(i in 1:pairwise.interaction.number){
      ind.covar <- ind.covar+1
      second.stage.mat[1:pairwise.interaction.second.cat,ind.covar]<-
        delta.no.inter[(ind.delta+1):(ind.delta+pairwise.interaction.second.cat)]
      ind.delta <- ind.delta + pairwise.interaction.second.cat
    }
  }
  
  if(saturated.number!=0){
    
    for(i in 1:saturated.number){
      ind.covar <- ind.covar+1
      second.stage.mat[1:saturated.second.cat,ind.covar]<-
        delta.no.inter[(ind.delta+1):(ind.delta+saturated.second.cat)]
      ind.delta <- ind.delta + saturated.second.cat
    }
  }
  
  colnames(second.stage.mat) <- covar.names
  rownames(second.stage.mat) <- full.second.stage.names[1:max.number.second.stage.parameter]
  places = 3
  second.stage.mat <- round(second.stage.mat,digits = 3)
  return(second.stage.mat)
  
}
