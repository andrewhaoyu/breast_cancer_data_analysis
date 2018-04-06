library(data.table)
setwd('/spin1/users/zhangh24/breast_cancer_data_analysis/')
discovery_snp <- as.data.frame(fread("./data/discovery_snps_annotated_clean.csv",header=T))
fine_mapping <- as.data.frame(fread("./data/fine_mapping_annotated_clean.csv"))

dis <- 2*10^6
check.data <- NULL
check.id <- NULL
for(i in 11:28){
  pos <- discovery_snp$position[i]
  CHR <- discovery_snp$CHR[i]
  for(j in 1:178){
    pos.known <- fine_mapping$position[j]
    CHR.known <- fine_mapping$CHR[j]
    if(CHR==CHR.known&pos>=(pos.known-dis)&pos<=(pos.known+dis)){
      print(c(i-10,j))
      check.id <- rbind(check.id,c(i,j))
      temp1 <- discovery_snp[i,c(1,3,2)]
      colnames(temp1) <- c("SNP","CHR","Position")
      temp2 <- fine_mapping[j,c(1,3,4)]
      colnames(temp2) <- c("SNP","CHR","Position")
      result <- rbind(temp1,temp2)
      check.data <- rbind(check.data,result)      
    }
      
  }
}

write.csv(check.data,file="./data/check_SNPs.csv",row.names = F,col.names = T)

data1 <- fread("./data/iCOGS_euro_v10_10232017.csv",header=T)
data1 <- as.data.frame(data1)
discovery.snp.icog <- as.data.frame(fread("/spin1/users/zhangh24/breast_cancer_data_analysis/data/discovery_icog_data.csv",header=T))
y.pheno.mis1 <- cbind(data1$Behaviour1,data1$ER_status1,data1$PR_status1,data1$HER2_status1,data1$Grade1)
colnames(y.pheno.mis1) = c("Behavior","ER","PR","HER2","Grade")
# Grade1.fake <- data1$Grade1
# Grade1.fake[data1$Grade1==2|data1$Grade1==3] <- 1
# Grade1.fake[data1$Grade1==1] <- 0
#y.pheno.mis1 <- cbind(data1$Behaviour1,data1$PR_status1,data1$ER_status1,data1$HER2_status1,Grade1.fake)
# y.pheno.mis1 <- cbind(data1$Behaviour1,data1$PR_status1,data1$ER_status1,data1$HER2_status1)

#x.test.all.mis1 <- data1[,c(27:203)]
###pc1-10 and age
i1 <- 1
idx.known <- which(colnames(data1)==check.data[2*i1,1])
x.covar1 <- cbind(data1[,c(5:14)],data1[,idx.known])

gene1 <- discovery.snp.icog[,check.id[i1,1]]
age <- data1[,204]
x.covar1 <- cbind(x.covar1,age)
idx.complete <- which(age!=888)
y.pheno.mis1 <- y.pheno.mis1[idx.complete,]
x.covar1 <- x.covar1[idx.complete,]
gene1 <- gene1[idx.complete]


data2 <- fread("./data/Onco_euro_v10_10232017.csv",header=T)
discovery.snp.onco <- as.data.frame(fread("/spin1/users/zhangh24/breast_cancer_data_analysis/data/discovery_onco_data.csv"))

data2 <- as.data.frame(data2)
y.pheno.mis2 <- cbind(data2$Behaviour1,data2$ER_status1,data2$PR_status1,data2$HER2_status1,data2$Grade1)
#y.pheno.mis2 <- cbind(data2$Behaviour1,data2$PR_status1,data2$ER_status1,data2$HER2_status1)
colnames(y.pheno.mis2) = c("Behaviour","ER",
                           "PR","HER2","Grade")

# x.test.all.mis2 <- data2[,c(27:203)]
idx.known <- which(colnames(data2)==check.data[2*i1,1])
x.covar2 <- cbind(data2[,c(5:14)],data2[,idx.known])
age <- data2[,204]
x.covar2 <- cbind(x.covar2,age)
idx.complete <- which(age!=888)
gene2 <- discovery.snp.onco[,check.id[i1,1]]
y.pheno.mis2 <- y.pheno.mis2[idx.complete,]
x.covar2 <- x.covar2[idx.complete,]
gene2 <- gene2[idx.complete]




two_data_two_stage_fixed <- function(
  y.pheno.mis1,
  gene1,
  x.covar1,
  y.pheno.mis2,
  gene2,
  x.covar2
){
  
  Heter.result1 = EMmvpoly(y.pheno.mis1,baselineonly = NULL,additive = cbind(gene1,x.covar1),pairwise.interaction = NULL,saturated = NULL,missingTumorIndicator = 888)
  z.standard <- Heter.result1[[12]]
  M <- nrow(z.standard)
  number.of.tumor <- ncol(z.standard)
  log.odds1 <- Heter.result1[[1]][(M+1):(M+1+number.of.tumor)]
  sigma.log.odds1 <- Heter.result1[[2]][(M+1):(M+1+number.of.tumor),(M+1):(M+1+number.of.tumor)]
  Heter.result2 = EMmvpoly(y.pheno.mis2,baselineonly = NULL,additive = cbind(gene2,x.covar2),pairwise.interaction = NULL,saturated = NULL,missingTumorIndicator = 888)
  
  M <- nrow(z.standard)
  number.of.tumor <- ncol(z.standard)
  log.odds2 <- Heter.result2[[1]][(M+1):(M+1+number.of.tumor)]
  sigma.log.odds2 <- Heter.result2[[2]][(M+1):(M+1+number.of.tumor),(M+1):(M+1+number.of.tumor)]
  meta.result <- LogoddsMetaAnalysis(log.odds1,
                                     sigma.log.odds1,
                                     log.odds2,
                                     sigma.log.odds2)
  second.stage.logodds.meta <- meta.result[[1]]
  second.stage.sigma.meta <- meta.result[[2]]
  test.result.second.wald <- DisplaySecondStageTestResult(second.stage.logodds.meta,second.stage.sigma.meta)
  p.value <- test.result.second.wald[length(test.result.second.wald)-1]
  return(list(meta.result,globalp = p.value))
}
two_data_two_stage_random <- function(y.pheno.mis1,
                                       gene1,
                                       x.covar1,
                                       y.pheno.mis2,
                                       gene2,
                                       x.covar2){
  
  z.standard <- GenerateZstandard(y.pheno.mis1)
  number.of.tumor <- ncol(z.standard)
  z.random.support <- cbind(1,z.standard[,1])
  z.random.test <- cbind(1,z.standard[,2:number.of.tumor])
  score.test.support.fixed1 <- ScoreTestSupport(
    y.pheno.mis1,
    baselineonly = NULL,
    additive = x.covar1,
    pairwise.interaction = NULL,
    saturated = NULL,
    missingTumorIndicator = 888
  )
  score.test.fixed1<- ScoreTest(y=y.pheno.mis1,
                              x=as.matrix(gene1),
                              second.stage.structure="additive",
                              score.test.support=  score.test.support.fixed1,
                              missingTumorIndicator=888)
  score.fixed1 <- score.test.fixed1[[1]]
  infor.fixed1 <- score.test.fixed1[[2]]
 
  score.test.support.random1 <- ScoreTestSupportSelfDesign(
    y.pheno.mis1,
    x.self.design  = as.matrix(gene1),
    z.design = z.random.support,
    additive = x.covar1,
    pairwise.interaction = NULL,
    saturated = NULL,
    missingTumorIndicator = 888
  )
  score.test.random1<- ScoreTestSelfDesign(y=y.pheno.mis1,
                                                    x= as.matrix(gene1),
                                                    z.design=z.random.test,
                                                    score.test.support=score.test.support.random1,
                                                    missingTumorIndicator=888)
  
}
