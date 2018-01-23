y=y.pheno.mis1.sub
baselineonly=NULL
additive=x.all.mis1.sub
pairwise.interaction=NULL
saturated=NULL
missingTumorIndicator = 888
delta0= NULL

y=y.pheno.mis2.sub
baselineonly=NULL
additive=x.all.mis2.sub
pairwise.interaction=NULL
saturated=NULL
missingTumorIndicator = 888
delta0= NULL


missing.data.vec <- GenerateMissingPosition(y,missingTumorIndicator)
y.pheno.complete <- y[-missing.data.vec,]
y.pheno.mis <- y[missing.data.vec,]
head(y.pheno.mis)
x.all.mis1.sub.complete <- x.all.mis1.sub[-missing.data.vec,]

initial.set <- InitialSetup(y.pheno.complete,
                            baselineonly,
                            additive,
                            pairwise.interaction,
                            saturated
)
###z standard matrix means the additive model z design matrix without baseline effect
###z standard matrix is used to match the missing tumor characteristics to the complete subtypes
if(is.null(delta0)==T){
  delta0 = initial.set$delta0
}else{
  delta0 = delta0
}

z.all = initial.set$z.all
z.standard = initial.set$z.standard
z.deisign.baselineonly = initial.set$z.design.baseline.only
z.design.additive = initial.set$z.design.additive
z.design.pairwise.interaction = initial.set$z.design.pairwise.interaction
z.design.saturated = initial.set$z.design.saturated
x.all <- as.matrix(GenerateXAll(y,baselineonly,additive,pairwise.interaction,saturated))
covar.names <- initial.set$covar.names
tumor.names <- initial.set$tumor.names


model.result = EMStep(delta0,as.matrix(y),x.all,z.standard,z.all,missingTumorIndicator)

summary.result <- SummaryResult(model.result,
                                baselineonly,
                                additive,
                                pairwise.interaction,
                                saturated,
                                z.standard,
                                covar.names,
                                delta,
                                z.design.additive,
                                z.design.pairwise.interaction,
                                z.design.saturated,
                                tumor.names,
                                z.all
)

