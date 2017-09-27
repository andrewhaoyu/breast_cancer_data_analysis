library(bc2)
##tp53

data1 <- read.csv("./data/iCOGS_euro_v10_05242017.csv",header=T)
data2 <- read.csv("./data/Onco_euro_v10_05242017.csv",header=T)
x.test.all.mis.icog = read.csv("./data/LD.position.knownSNP.pruning.SNP.icog.value.csv",header=T)
x.test.all.mis.onco = read.csv("./data/LD.position.knownSNP.pruning.SNP.onco.value.csv",header=T)
y.pheno.mis.icog <- cbind(data1$Behaviour1,data1$PR_status1,data1$ER_status1,data1$HER2_status1,data1$Grade1)
y.pheno.mis.onco <- cbind(data2$Behaviour1,data2$PR_status1,data2$ER_status1,data2$HER2_status1,data2$Grade1)
snpname <- "chr17_7571752_G_T"
SNP.id.icog <- colnames(x.test.all.mis.icog)
idx.icog <- grep(snpname,SNP.id.icog)
i1 <- idx.icog
x.covar.mis1 <- data1[,5:14]
x.all.mis1 <- as.matrix(cbind(x.test.all.mis.icog[,i1],x.covar.mis1))
y.pheno.mis1 <- cbind(data1$Behaviour1,data1$PR_status1,data1$ER_status1,data1$HER2_status1,data1$Grade1)
colnames(y.pheno.mis1) = c("Behavior","PR","ER","HER2","Grade")

Heter.result.Icog = EMmvpoly(y.pheno.mis1,baselineonly = NULL,additive = x.all.mis1,pairwise.interaction = NULL,saturated = NULL,missingTumorIndicator = 888)
M <- 23
z.standard <- Heter.result.Icog[[12]]
z.additive.design <- as.matrix(cbind(1,z.standard))
M <- nrow(z.standard)
number.of.tumor <- ncol(z.standard)
log.odds.icog <- Heter.result.Icog[[1]][(M+1):(M+1+number.of.tumor)]

sigma.log.odds.icog <- Heter.result.Icog[[2]][(M+1):(M+1+number.of.tumor),(M+1):(M+1+number.of.tumor)]
beta.icog <- z.additive.design%*%log.odds.icog
beta.sigma.icog <- z.additive.design%*%sigma.log.odds.icog%*%t(z.additive.design)
loglikelihood.icog <- Heter.result.Icog[[8]]







