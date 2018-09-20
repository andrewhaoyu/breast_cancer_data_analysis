z.design <- matrix(c(
  c(0,1,1,1,0,0,0,0,1,1,1,0,0,0,0,0,1,1,0,0,0,0,0),
  c(0,0,0,0,0,1,1,0,0,0,0,0,1,1,1,0,0,0,0,0,1,1,1),
  c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0),
  c(0,0,0,0,1,0,0,0,0,0,0,1,0,0,0,0,0,0,0,1,0,0,0),
  c(1,0,0,0,0,0,0,1,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0)
),ncol=5)


i1 = 1
library(bc2)
install_github("andrewhaoyu/bc2", ref = "development",args = c('--library="/home/zhangh24/R/x86_64-pc-linux-gnu-library/3.4"'))
library(bc2)
data1 <- read.csv("./data/iCOGS_euro_v10_05242017.csv",header=T)
y.pheno.mis1 <- cbind(data1$Behaviour1,data1$PR_status1,data1$ER_status1,data1$HER2_status1,data1$Grade1)
colnames(y.pheno.mis1) = c("Behavior","PR","ER","HER2","Grade")
# Grade1.fake <- data1$Grade1
# Grade1.fake[data1$Grade1==2|data1$Grade1==3] <- 1
# Grade1.fake[data1$Grade1==1] <- 0
#y.pheno.mis1 <- cbind(data1$Behaviour1,data1$PR_status1,data1$ER_status1,data1$HER2_status1,Grade1.fake)
# y.pheno.mis1 <- cbind(data1$Behaviour1,data1$PR_status1,data1$ER_status1,data1$HER2_status1)

x.test.all.mis1 <- data1[,c(27:206)]

x.covar.mis1 <- data1[,5:14]
x.all.mis1 <- as.matrix(cbind(x.test.all.mis1[,i1],x.covar.mis1))
colnames(x.all.mis1)[1] <- "gene"

Heter.result.Icog = EMmvpoly(y.pheno.mis1,baselineonly = NULL,additive = x.all.mis1[,2:11],pairwise.interaction = NULL,saturated = NULL,missingTumorIndicator = 888)

M <- 23
delta0.icog.temp <- Heter.result.Icog[[1]]
delta0.icog <- c(delta0.icog.temp[1:M],0,delta0.icog.temp[(M+1):length(delta0.icog.temp)])
z.standard <- Heter.result.Icog[[12]]

save(delta0.icog,file = "./whole_genome/ICOG/ERPRHER2GRADE_fixed_baseline/result/delta0.icog.Rdata")
save(z.standard,file = "./whole_genome/ICOG/ERPRHER2GRADE_fixed_baseline/result/z.standard.Rdata")


z.standard <- Heter.result.Icog[[12]]
z.additive.design <- as.matrix(cbind(1,z.standard))
M <- nrow(z.standard)
number.of.tumor <- ncol(z.standard)


score.test.support.icog.mix <- ScoreTestSupportMixedModel(
  y.pheno.mis1,
  z.design = 
  baselineonly = NULL,
  additive = x.all.mis1[,2:3],
  pairwise.interaction = NULL,
  saturated = NULL,
  missingTumorIndicator = 888
)
ScoreTestMixedModel(y.pheno.mis1,x.all.mis1[,1,drop=F],z.standard,score.test.support = score.test.support.icog.mix,missingTumorIndicator = 888)



score.test.support.icog <- ScoreTestSupportMixedModel(
  y.pheno.mis1,
  z.design = z.design,
  x.self.design = x.all.mis1[,1,drop=F],
  baselineonly = NULL,
  additive = x.all.mis1[,2:3],
  pairwise.interaction = NULL,
  saturated = NULL,
  missingTumorIndicator = 888
)
ScoreTestMixedModel(y.pheno.mis1,x.all.mis1[,1,drop=F],second.stage.structure = "additive",score.test.support = score.test.support.icog,missingTumorIndicator = 888)




score.test.icog<- ScoreTestSelfDesign(y=y.pheno.mis1,
                                      x=x.all.mis1[,1,drop=F],
                                      z.design=z.design,
                                      score.test.support=score.test.support.icog,
                                      missingTumorIndicator=888)


data2 <- read.csv("./data/Onco_euro_v10_rs554219.csv",header = T)
#data2 <- read.csv("/data/zhangh20/breast_cancer/known_SNPs_anlysis/Onco_euro_v10_rs554219.csv",1)

# names1 = colnames(data1)[27:206]
#rm(data1)
names2 = colnames(data2)

idxi1 = which(names2=="rs554219")
y.pheno.mis2 <- cbind(data2$Behaviour1,data2$PR_status1,data2$ER_status1,data2$HER2_status1,data2$Grade1)
colnames(y.pheno.mis2) <- c("Bahavior1","PR","ER","HER2","Grade")
x.test.all.mis2 <- data2

x.covar.mis2 <- data2[,5:14]
x.all.mis2 <- as.matrix(cbind(x.test.all.mis2[,idxi1],x.covar.mis2))



Heter.result.Onco = EMmvpoly(y.pheno.mis2,baselineonly = NULL,additive = x.all.mis2[,2:11],pairwise.interaction = NULL,saturated = NULL,missingTumorIndicator = 888)
z.standard <- Heter.result.Onco[[12]]
M <- 23
delta0.onco.temp <- Heter.result.Onco[[1]]
delta0.onco <- c(delta0.onco.temp[1:M],0,delta0.onco.temp[(M+1):length(delta0.onco.temp)])
z.standard <- Heter.result.Onco[[12]]

save(delta0.onco,file = "./whole_genome/ONCO/ERPRHER2GRADE_fixed_baseline/result/delta0.onco.Rdata")
save(z.standard,file = "./whole_genome/ONCO/ERPRHER2GRADE_fixed_baseline/result/z.standard.Rdata")
