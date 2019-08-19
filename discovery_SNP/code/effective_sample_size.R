#Goal: Generate effective sample size


setwd('/dcl01/chatterj/data/hzhang1/breast_cancer_data_analysis')
#load summary level statistics for subtypes
load("/dcl01/chatterj/data/hzhang1/breast_intrinsic/whole_genome_breast_cancer_results/BCAC_subtypes_result.Rdata")

EffectiveCases <- function(N_control,effective_sample_size){
  return(N_control*effective_sample_size/(N_control-effective_sample_size))
}

idx <- grep("effect_sample_size",colnames(BCAC_subtypes_result))
colnames(BCAC_subtypes_result)[idx]
N_control <- 91477 


M <- 5
result <- rep(0,M)
for(i in 1:M){
temp <- EffectiveCases(N_control,BCAC_subtypes_result[,idx[i]])  
result[i] <- mean(temp)
}
names(result) <- colnames(BCAC_subtypes_result)[idx]



IntrinsicSampleSize <- function(y.pheno.tumor1){
  idx.1 <- which((y.pheno.tumor1[,1]==1|y.pheno.tumor1[,2]==1)
                 &y.pheno.tumor1[,3]==0
                 &(y.pheno.tumor1[,4]==1|y.pheno.tumor1[,4]==2))
  
  #for second subtype HR+_HER2-_highgrade 
  idx.2 <- which((y.pheno.tumor1[,1]==1|y.pheno.tumor1[,2]==1)
                 &y.pheno.tumor1[,3]==0
                 &y.pheno.tumor1[,4]==3)
  
  #for third subtype HR+_HER2+
  idx.3 <- which((y.pheno.tumor1[,1]==1|y.pheno.tumor1[,2]==1)
                 &y.pheno.tumor1[,3]==1)
  
  
  #for third subtype HR-_HER2+
  idx.4 <- which(y.pheno.tumor1[,1]==0&y.pheno.tumor1[,2]==0
                 &y.pheno.tumor1[,3]==1)
  
  #for third subtype HR-_HER2-
  idx.5 <- which(y.pheno.tumor1[,1]==0&y.pheno.tumor1[,2]==0
                 &y.pheno.tumor1[,3]==0)
  result <- c(length(idx.1),length(idx.2),length(idx.3),
              length(idx.4),length(idx.5))
  return(result)
}







data1 <- fread("./data/iCOGS_euro_v10_10232017.csv",header=T)
data1 <- as.data.frame(data1)
y.pheno.mis1 <- cbind(data1$Behaviour1,data1$ER_status1,data1$PR_status1,data1$HER2_status1,data1$Grade1)
colnames(y.pheno.mis1) = c("Behavior","ER","PR","HER2","Grade")
age <- data1[,204]
age <- data1[,204]
idx.complete <- which(age!=888)
y.pheno.mis1 <- y.pheno.mis1[idx.complete,]
y.pheno.tumor1 <- y.pheno.mis1[,2:5]

IntrinsicSampleSize(y.pheno.tumor1)


data2 <- fread("./data/Onco_euro_v10_10232017.csv",header=T)
data2 <- as.data.frame(data2)
y.pheno.mis2 <- cbind(data2$Behaviour1,data2$ER_status1,data2$PR_status1,data2$HER2_status1,data2$Grade1)
colnames(y.pheno.mis2) = c("Behaviour","ER",
                           "PR","HER2","Grade")
ages <- data2[,204]
idx.complete <- which(ages!=888)

y.pheno.mis2 <- y.pheno.mis2[idx.complete,]
y.pheno.tumor2 <- y.pheno.mis2[,2:5]
IntrinsicSampleSize(y.pheno.tumor2)
IntrinsicSampleSize(y.pheno.tumor2)+IntrinsicSampleSize(y.pheno.tumor1)


IntrinsicSampleSize(y.pheno.tumor2)+IntrinsicSampleSize(y.pheno.tumor1)