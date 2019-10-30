Heter.result.Icog = EMmvpolySelfDesign(y.pheno.mis1.train,x.self.design = x.snp.all.train1,z.design=z.design,baselineonly = NULL,additive = x.covar.train1,pairwise.interaction = NULL,saturated = NULL,missingTumorIndicator = 888)

y = y.pheno.mis1.train;
x.self.design = x.snp.all.train1[,i1:2];
z.design=z.design;
baselineonly = NULL;
additive = x.covar.train1;
pairwise.interaction = NULL;
saturated = NULL;
missingTumorIndicator = 888
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


M <- nrow(z.design.additive)
##the number of covariates in different potential model structures
self.design.number <- CountCovarNumber(x.self.design)
baselineonly.number <- CountCovarNumber(baselineonly)
additive.number <- CountCovarNumber(additive)
pairwise.interaction.number <- CountCovarNumber(pairwise.interaction)
saturated.number <- CountCovarNumber(saturated)
###second.stage.category for different model structures
self.design.second.cat <- ncol(z.design)
baselineonly.second.cat <- 1
additive.second.cat <- ncol(z.design.additive)
pairwise.interaction.second.cat <- ncol(z.design.pairwise.interaction)
saturated.second.cat <- ncol(z.design.saturated)
###1 for intercept
total.covar.number <- 1+ self.design.number+baselineonly.number+additive.number+
  pairwise.interaction.number+saturated.number

z.all <- matrix(0,nrow=(M*total.covar.number),ncol = (M+self.design.second.cat*self.design.number+baselineonly.number*baselineonly.second.cat+
                                                        additive.second.cat*additive.number+
                                                        pairwise.interaction.second.cat*pairwise.interaction.number)+saturated.second.cat*saturated.number)

for(i in c("intercept","self design","baselineonly",
           "additive","pairwise.interaction",
           "satuared")){
  ##we always keep intercept as saturated model and to simply, we always use diagnonal matrix for intercept
  if(i=="intercept"){
    ###row start and column start point for this category
    row.start <- 0
    column.start <- 0
    for(j in 1:M){
      z.all[row.start+1+(j-1)*total.covar.number,(column.start+j)] = 1
    }
  }else if(i=="self design"){
    column.start <- M
    row.start <- 1
    if(self.design.number!=0){
      for(j in 1:M){
        for(k in 1:self.design.number){
          z.all[row.start+k+(j-1)*total.covar.number,
                (column.start+(k-1)*self.design.second.cat+1):
                  (column.start+k*self.design.second.cat)] <- as.matrix(z.design[j,])
        }
      }
    }
  }else if(i=="baselineonly"){
    column.start = M+self.design.number*self.design.second.cat
    row.start <- 1+self.design.number
    ###test whether there is any baselineonly variable
    if(baselineonly.number!=0){
      for(j in 1:M){
        for(k in 1:baselineonly.number){
          z.all[row.start+k+(j-1)*total.covar.number,
                (column.start+(k-1)*baselineonly.second.cat+1):
                  (column.start+k*baselineonly.second.cat)] <- as.matrix(z.design.baselineonly[j,])
        }
      }
    }
  }else if(i=="additive"){
    column.start <- M+self.design.number*self.design.second.cat+baselineonly.number
    row.start <- 1+self.design.number+baselineonly.number
    if(additive.number!=0){
      for(j in 1:M){
        for(k in 1:additive.number){
          z.all[row.start+k+(j-1)*total.covar.number,
                (column.start+(k-1)*additive.second.cat+1):
                  (column.start+k*additive.second.cat)] <- as.matrix(z.design.additive[j,])
        }
      }
    }
  }else if(i == "pairwise.interaction"){
    column.start <- M+self.design.number*self.design.second.cat+baselineonly.number+additive.number*additive.second.cat
    row.start <- 1+self.design.number+baselineonly.number+additive.number
    if(pairwise.interaction.number!=0){
      for(j in 1:M){
        for(k in 1:pairwise.interaction.number){
          z.all[row.start+k+(j-1)*total.covar.number,
                (column.start+(k-1)*pairwise.interaction.second.cat+1):
                  (column.start+k*pairwise.interaction.second.cat)] <- as.matrix(z.design.pairwise.interaction[j,])
        }
      }
    }
  }else {
    column.start <- M+self.design.number*self.design.second.cat+baselineonly.number+additive.number*additive.second.cat+
      pairwise.interaction.number*pairwise.interaction.second.cat
    row.start <- 1+self.design.number+baselineonly.number+additive.number+pairwise.interaction.number
    if(saturated.number!=0){
      for(j in 1:M){
        for(k in 1:saturated.number){
          z.all[row.start+k+(j-1)*total.covar.number,
                (column.start+(k-1)*saturated.second.cat+1):
                  (column.start+k*saturated.second.cat)] <- as.matrix(z.design.saturated[j,])
        }
      }
    }
  }
  
  
  
}
return(z.all)












delta0 <-StartValueFunction(freq.subtypes,y.case.control,z.all)
z.standard <- z.design.additive[,-1]
delta0 <-StartValueFunction(freq.subtypes,y.case.control,z.all)
