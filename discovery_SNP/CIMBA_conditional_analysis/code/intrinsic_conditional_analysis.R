arg <- commandArgs(trailingOnly=T)
i1 <- as.numeric(arg[[1]])
print(i1)
setwd("/data/zhangh24/breast_cancer_data_analysis/")
library(bc2)

z.design <- matrix(c(
  c(0,1,1,1,0,0,0,0,1,1,1,0,0,0,0,0,1,1,0,0,0,0,0),
  c(0,0,0,0,0,1,1,0,0,0,0,0,1,1,1,0,0,0,0,0,1,1,1),
  c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0),
  c(0,0,0,0,1,0,0,0,0,0,0,1,0,0,0,0,0,0,0,1,0,0,0),
  c(1,0,0,0,0,0,0,1,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0)
),ncol=5)
rowSums(z.design)



colnames(z.design) <- c("Luminial A","Luminal B",
                        "Luminal B HER2Neg",
                        "HER2 Enriched",
                        "Triple Negative")

#Icog.order <- read.table(gzfile(subject.file))
load("./discovery_SNP/CIMBA_conditional_analysis/result/ICOG/CIMBA_icog_snpvalue.Rdata")
load("./discovery_SNP/CIMBA_conditional_analysis/result/ONCO/CIMBA_onco_snpvalue.Rdata")

library(data.table)
setwd("/data/zhangh24/breast_cancer_data_analysis/")
data1 <- fread("./data/iCOGS_euro_v10_10232017.csv",header=T)
data1 <- as.data.frame(data1)
y.pheno.mis1 <- cbind(data1$Behaviour1,data1$ER_status1,data1$PR_status1,data1$HER2_status1,data1$Grade1)
colnames(y.pheno.mis1) = c("Behavior","ER","PR","HER2","Grade")
#x.test.all.mis1 <- data1[,c(27:206)]
SG_ID <- data1$SG_ID
x.covar.mis1 <- data1[,c(5:14,153,204)]
age <- data1[,204]
idx.complete1 <- which(age!=888)
y.pheno.mis1 <- y.pheno.mis1[idx.complete1,]
x.covar.mis1 <- x.covar.mis1[idx.complete1,]
SG_ID <- SG_ID[idx.complete1]


data2 <- fread("./data/Onco_euro_v10_10232017.csv",header=T)
data2 <- as.data.frame(data2)
y.pheno.mis2 <- cbind(data2$Behaviour1,data2$ER_status1,data2$PR_status1,data2$HER2_status1,data2$Grade1)
#y.pheno.mis2 <- cbind(data2$Behaviour1,data2$PR_status1,data2$ER_status1,data2$HER2_status1)
colnames(y.pheno.mis2) = c("Behaviour","ER",
                           "PR","HER2","Grade")

x.covar.mis2 <- data2[,c(5:14,153,204)]
ages <- data2[,204]
idx.complete2 <- which(ages!=888)

y.pheno.mis2 <- y.pheno.mis2[idx.complete2,]
x.covar.mis2 <- x.covar.mis2[idx.complete2,]
Onc_ID <- data2$Onc_ID
Onc_ID <- Onc_ID[idx.complete2]

conditional_result <- NULL

for(i2 in (5*i1-4):(5*i1)){
  if(i2<=4120){
    snpvalue = snpvalue.result.icog[idx.complete1,i2]
    Heter.result.Icog = EMmvpolySelfDesign(y.pheno.mis1,x.self.design = snpvalue,z.design=z.design,baselineonly = NULL,additive = x.covar.mis1,pairwise.interaction = NULL,saturated = NULL,missingTumorIndicator = 888)
    z.standard <- Heter.result.Icog[[12]]
    M <- nrow(z.standard)
    number.of.tumor <- ncol(z.standard)
    log.odds.icog <- Heter.result.Icog[[1]][(M+1):(M+1+number.of.tumor)]
    nparm <- length(Heter.result.Icog[[1]])  
    sigma.log.odds.icog <- Heter.result.Icog[[2]][(M+1):(M+1+number.of.tumor),(M+1):(M+1+number.of.tumor)]
    snpvalue = snpvalue.result.onco[idx.complete2,i2]
    Heter.result.Onco = EMmvpolySelfDesign(y.pheno.mis2,x.self.design = snpvalue,z.design = z.design,baselineonly = NULL,additive = x.covar.mis2,pairwise.interaction = NULL,saturated = NULL,missingTumorIndicator = 888)
    z.standard <- Heter.result.Onco[[12]]
    M <- nrow(z.standard)
    number.of.tumor <- ncol(z.standard)
    log.odds.onco <- Heter.result.Onco[[1]][(M+1):(M+1+number.of.tumor)]
    nparm <- length(Heter.result.Onco[[1]])
    sigma.log.odds.onco <- Heter.result.Onco[[2]][(M+1):(M+1+number.of.tumor),(M+1):(M+1+number.of.tumor)]
    
    
    
    
    
    meta.result <- LogoddsMetaAnalysis(log.odds.icog,
                        sigma.log.odds.icog,
                        log.odds.onco,
                        sigma.log.odds.onco)
    temp <- c(as.numeric(meta.result[[1]]),as.numeric(meta.result[[2]]))
    
    conditional_result <- rbind(conditional_result,temp)
    
    
    
    }else if(i2==4121){
      snpvalue = snpvalue.result.onco[idx.complete2,i2]
      Heter.result.Onco = EMmvpolySelfDesign(y.pheno.mis2,x.self.design = snpvalue,z.design = z.design,baselineonly = NULL,additive = x.covar.mis2,pairwise.interaction = NULL,saturated = NULL,missingTumorIndicator = 888)
      z.standard <- Heter.result.Onco[[12]]
      M <- nrow(z.standard)
      number.of.tumor <- ncol(z.standard)
      log.odds.onco <- Heter.result.Onco[[1]][(M+1):(M+1+number.of.tumor)]
      nparm <- length(Heter.result.Onco[[1]])
      sigma.log.odds.onco <- Heter.result.Onco[[2]][(M+1):(M+1+number.of.tumor),(M+1):(M+1+number.of.tumor)]
      temp <- c(as.numeric(log.odds.onco),as.numeric(sigma.log.odds.onco))
      conditional_result <- rbind(conditional_result,temp)
      
  }else{}
   
}
save(conditional_result,file = paste0("/data/zhangh24/breast_cancer_data_analysis/discovery_SNP/CIMBA_conditional_analysis/result/conditional_result",i1,".Rdata"))
