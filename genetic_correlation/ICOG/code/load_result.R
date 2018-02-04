setwd('/Users/zhangh24/GoogleDrive/breast_cancer_data_analysis/genetic_correlation/ICOG/result/')
load("/Users/zhangh24/GoogleDrive/breast_cancer_data_analysis/genetic_correlation/ICOG/result/ICOG_ldsc.rda")
covariance.matrix <- ICOG_ldsc[[1]]
covariance.matrix.se <- ICOG_ldsc[[3]]
correlation.matrix <- ICOG_ldsc[[2]]
correlation.matrix.se <- ICOG_ldsc[[4]]








covariance.matrix.new <- rep("c",25)

temp <- 1
for(i in 1:5){
  for(j in 1:5){
    covariance.matrix.new[temp] <- paste0(covariance.matrix[i,j]," (",covariance.matrix[i,j]-1.96*covariance.matrix.se[i,j],"-",covariance.matrix[i,j]+1.96*covariance.matrix.se[i,j],")") 
    temp <- temp+1
  }
}
covariance.matrix.new <- matrix(covariance.matrix.new,5,5)



correlation.matrix.new <- rep("c",25)

temp <- 1
for(i in 1:5){
  for(j in 1:5){
    correlation.matrix.new[temp] <- paste0(correlation.matrix[i,j]," (",correlation.matrix[i,j]-1.96*correlation.matrix.se[i,j],"-",correlation.matrix[i,j]+1.96*correlation.matrix.se[i,j],")") 
    temp <- temp+1
  }
}
correlation.matrix.new <- matrix(correlation.matrix.new,5,5)





write.csv(covariance.matrix.new,file="./covariance_matrix_icog.csv",quote=F)
write.csv(correlation.matrix.new,file="./correlation_matrix_icog.csv",quote=F)

library(tidyverse)
library(reshape2)
correlation.matrix.com <- melt(correlation.matrix)
ggplot(data = correlation.matrix.com, aes(Var1,Var2,fill=value))+
  geom_tile()+
  theme_minimal()+
 scale_fill_gradient(low="white",high="dodgerblue4")
library(ggcorrplot)
ggcorrplot(correlation.matrix)

  




