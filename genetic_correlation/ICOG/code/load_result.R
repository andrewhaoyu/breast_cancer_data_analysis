get_cor_95 <- function(correlation.matrix,correlation.matrix.se){
  correlation.matrix.transform <- 0.5*log((1+correlation.matrix)/(1-correlation.matrix))
  correlation.matrix.se.transform <- abs(correlation.matrix.se/((1+correlation.matrix)*(1-correlation.matrix)))
  
  correlation.matrix.transform.low <- correlation.matrix.transform-1.96*correlation.matrix.se.transform
  correlation.matrix.transform.high <- correlation.matrix.transform+1.96*correlation.matrix.se.transform
  
  correlation.matrix.low <- (exp(2*correlation.matrix.transform.low)-1)/(exp(2*correlation.matrix.transform.low)+1)
  correlation.matrix.high <- (exp(2*correlation.matrix.transform.high)-1)/(exp(2*correlation.matrix.transform.high)+1)
  
  (exp(2*correlation.matrix.transform)-1)/(exp(2*correlation.matrix.transform)+1)
  
  correlation.matrix <- round(correlation.matrix,places)
  correlation.matrix.low <- round(correlation.matrix.low,places)
  correlation.matrix.high <- round(correlation.matrix.high,places)
  
  M <- nrow(correlation.matrix)
  
  correlation.matrix.new <- rep("c",M^2)
  
  temp <- 1
  for(i in 1:M){
    for(j in 1:M){
      correlation.matrix.new[temp] <- paste0(correlation.matrix[i,j]," (",correlation.matrix.low[i,j],"-",correlation.matrix.high[i,j],")") 
      temp <- temp+1
    }
  }
  correlation.matrix.new <- matrix(correlation.matrix.new,M,M)
  return(correlation.matrix.new)  
}

setwd('/Users/zhangh24/GoogleDrive/breast_cancer_data_analysis/genetic_correlation/ICOG/result/breast_cancer_ldsc/two-stage-model-intrinsic_subtype/')
# load("./breast_cancer_ldsc/ICOG_ldsc_result.rda")
# load("./breast_cancer_ldsc/ONCO_ldsc_result.rda")
#load("./breast_cancer_ldsc/two-stage-model-intrinsic_subtype/meta_ldsc_result.rda")
intrinsic_result <- read.csv("two_stage_model_genetic_correlation_result.csv",header=F)
correlation.matrix <- intrinsic_result[1:5,2:6]
correlation.matrix.se <- intrinsic_result[6:10,2:6]

correlation.95 <- get_cor_95(correlation.matrix,correlation.matrix.se)
colnames(correlation.matrix) <- intrinsic_result[1:5,1]
rownames(correlation.matrix) <- intrinsic_result[1:5,1]

colnames(correlation.95) <- intrinsic_result[1:5,1]
rownames(correlation.95) <- intrinsic_result[1:5,1]
write.csv(correlation.95,file="covariance_95_two_stage.csv",quote=F)
places <- 3



# covariance.matrix.new <- rep("c",25)
# 
# 
# covariance.matrix.transform <- log(covariance.matrix)
# covariance.matrix.transform.low <- covariance.matrix.transform - 1.96*covariance.matrix.se/covariance.matrix
# covariance.matrix.transform.high <- covariance.matrix.transform + 1.96*covariance.matrix.se/covariance.matrix
# covariance.matrix.high <- exp(covariance.matrix.transform.high)
# covariance.matrix.low  <- exp(covariance.matrix.transform.low)
# covariance.matrix.high <- round(covariance.matrix.high,places)
# covariance.matrix.low <- round(covariance.matrix.low,places)
# 
# temp <- 1
# for(i in 1:5){
#   for(j in 1:5){
#     covariance.matrix.new[temp] <- paste0(covariance.matrix[i,j]," (",covariance.matrix.low[i,j],"-",covariance.matrix.high[i,j],")") 
#     temp <- temp+1
#   }
# }
# covariance.matrix.new <- matrix(covariance.matrix.new,5,5)







#covariance.matrix.icog <- covariance.matrix.new



# write.csv(covariance.matrix.icog,file="./covariance_matrix_icog.csv",quote=F)
# write.csv(correlation.matrix.new,file="./correlation_matrix_icog.csv",quote=F)

library(tidyverse)
library(reshape2)
library(ggplot2)
library(gplots)



# correlation.matrix <- as.matrix(as.data.frame(fread("genetic_correlation_meta.csv")[,1:5]))
#rownames(correlation.matrix) <- colnames(correlation.matrix)
correlation.matrix <- correlation.matrix[c(1,2,4,5),c(1,2,4,5)]
rownames(correlation.matrix) <- colnames(correlation.matrix)
correlation.matrix <- as.matrix(correlation.matrix)
#correlation.matrix <- correlation.matrix[c(2,4,3,1),c(2,4,3,1)]

library(corrplot)


pal.breaks <- seq(0.5,1,0.01)
col <- colorRampPalette(c("white","red"))(length(pal.breaks)-1)

# heatmap(correlation.matrix,col = col,symm=T,margins=c(10,4),key.title="",key.ylab="",cexRow=1,cexCol=1)
png(filename="./meta_heatmap_two_stage.png",width=10,heigh=10,units="in",res=600)
corrplot(correlation.matrix, method = "circle",order="hclust",
         addrect=2, tl.col = "black", tl.srt = 30)
# heatmap.2(correlation.matrix,tracecol=NA,cexRow=1,cexCol=1,margins=c(10,12),col = col,breaks=pal.breaks,key.ylab="",key.title = "",
#           main=" Genetic Correlation Heatmap",dendrogram="row",density.info="none",lwid = c(1.5,4))
dev.off()


ER_result <- as.data.frame(fread("ER_genetic_correlation_result.csv",header=F))
correlation.matrix <- ER_result[1:2,1:2]
correlation.matrix.se <- ER_result[3:4,1:2]

correlation.95 <- get_cor_95(correlation.matrix,correlation.matrix.se)
colnames(correlation.matrix) <- c("ER pos","ER neg")
rownames(correlation.matrix) <- c("ER pos","ER neg")

colnames(correlation.95) <- c("ER pos","ER neg")
rownames(correlation.95) <- c("ER pos","ER neg")
write.csv(correlation.95,file="covariance_95_ER.csv",quote=F)
places <- 3





# pal.breaks <- seq(-1,1,0.01)
# col <- colorRampPalette(c("dodgerblue4","white","red"))(length(pal.breaks)-1)
# correlation.matrix.icog <- correlation.matrix
# # heatmap(correlation.matrix,col = col,symm=T,margins=c(10,4),key.title="",key.ylab="",cexRow=1,cexCol=1)
# png(filename="./meta_heatmap2.png",width=10,heigh=10,units="in",res=300)
# heatmap.2(correlation.matrix,tracecol=NA,cexRow=1,cexCol=1,margins=c(10,12),col = col,breaks=pal.breaks,key.ylab="",key.title = "",
#           main=" Genetic Correlation Heatmap",dendrogram="row",density.info="none",lwid = c(1.5,4))
# dev.off()
# 
# 
# # correlation.matrix.com <- melt(correlation.matrix))
# # ggplot(data = correlation.matrix.com, aes(Var1,Var2,fill=value))+
# #   geom_tile()+
# #   theme_minimal()+
# #  scale_fill_gradient(low="white",high="dodgerblue4")
# # 
# # cluster.p <- as.dendrogram(hclust(d = dist(x = correlation.matrix)))
# # ggdendrogram(data=cluster.p,rotate=F)
# # library(ggcorrplot)
# # ggcorrplot(correlation.matrix)
# # 
# #   
# # 
# # 
# # 
# #rm(list=ls())
# setwd('/Users/zhangh24/GoogleDrive/breast_cancer_data_analysis/genetic_correlation/ICOG/result/')
# load("/Users/zhangh24/GoogleDrive/breast_cancer_data_analysis/genetic_correlation/ICOG/result/breast_cancer_ldsc/ONCO_ldsc_result.rda")
# covariance.matrix <- ldsc_result[[1]]
# covariance.matrix.se <- ldsc_result[[3]]
# correlation.matrix <- ldsc_result[[2]]
# correlation.matrix.se <- ldsc_result[[4]]
# 
# 
# 
# 
# places <- 3
# 
# 
# 
# covariance.matrix.new <- rep("c",25)
# 
# covariance.matrix.transform <- log(covariance.matrix)
# covariance.matrix.transform.low <- covariance.matrix.transform - 1.96*covariance.matrix.se/covariance.matrix
# covariance.matrix.transform.high <- covariance.matrix.transform + 1.96*covariance.matrix.se/covariance.matrix
# covariance.matrix.high <- exp(covariance.matrix.transform.high)
# covariance.matrix.transform.low  <- exp(covariance.matrix.transform.low)
# covariance.matrix.high <- round(covariance.matrix.high,places)
# covariance.matrix.low <- round(covariance.matrix.low,places)
# 
# temp <- 1
# for(i in 1:5){
#   for(j in 1:5){
#     covariance.matrix.new[temp] <- paste0(covariance.matrix[i,j]," (",covariance.matrix.low[i,j],"-",covariance.matrix.high[i,j],")") 
#     temp <- temp+1
#   }
# }
# covariance.matrix.new <- matrix(covariance.matrix.new,5,5)
# 
# 
# 
# 
# 
# 
# correlation.matrix.transform <- 0.5*log((1+correlation.matrix)/(1-correlation.matrix))
# correlation.matrix.se.transform <- abs(correlation.matrix.se/((1+correlation.matrix)*(1-correlation.matrix)))
# 
# correlation.matrix.transform.low <- correlation.matrix.transform-1.96*correlation.matrix.se.transform
# correlation.matrix.transform.high <- correlation.matrix.transform+1.96*correlation.matrix.se.transform
# 
# correlation.matrix.low <- (exp(2*correlation.matrix.transform.low)-1)/(exp(2*correlation.matrix.transform.low)+1)
# correlation.matrix.high <- (exp(2*correlation.matrix.transform.high)-1)/(exp(2*correlation.matrix.transform.high)+1)
# 
# (exp(2*correlation.matrix.transform)-1)/(exp(2*correlation.matrix.transform)+1)
# 
# correlation.matrix <- round(correlation.matrix,places)
# correlation.matrix.low <- round(correlation.matrix.low,places)
# correlation.matrix.high <- round(correlation.matrix.high,places)
# 
# 
# correlation.matrix.new <- rep("c",25)
# 
# temp <- 1
# for(i in 1:5){
#   for(j in 1:5){
#     correlation.matrix.new[temp] <- paste0(correlation.matrix[i,j]," (",correlation.matrix.low[i,j],"-",correlation.matrix.high[i,j],")") 
#     temp <- temp+1
#   }
# }
# correlation.matrix.new <- matrix(correlation.matrix.new,5,5)
# 
# covariance.matrix.onco <- covariance.matrix.new
# 
# 
# 
# write.csv(covariance.matrix.onco,file="./covariance_matrix_onco.csv",quote=F)
# write.csv(correlation.matrix.new,file="./correlation_matrix_onco.csv",quote=F)
# 
# library(tidyverse)
# library(reshape2)
# library(ggplot2)
# library(gplots)
# pal.breaks <- seq(-1,1,0.01)
# col <- colorRampPalette(c("dodgerblue4","white","red"))(length(pal.breaks)-1)
# # heatmap(correlation.matrix,col = col,symm=T,margins=c(10,4),key.title="",key.ylab="",cexRow=1,cexCol=1)
# correlation.matrix.onco<- correlation.matrix
# png(filename="./ONCO_heatmap.png",width=800,heigh=600,units="px")
# heatmap.2(correlation.matrix.onco,tracecol=NA,cexRow=1,cexCol=1,margins=c(10,12),col = col,breaks=pal.breaks,key.ylab="",key.title = "",
#           main=" Genetic Correlation Heatmap Onco",dendrogram="row",density.info="none",lwid = c(1.5,4))
# dev.off()
# 
# 
# 
# 
# #rm(list=ls())
# setwd('/Users/zhangh24/GoogleDrive/breast_cancer_data_analysis/genetic_correlation/ICOG/result/')
# load("/Users/zhangh24/GoogleDrive/breast_cancer_data_analysis/genetic_correlation/ICOG/result/breast_cancer_ldsc/meta_ldsc_result.rda")
# covariance.matrix <- ldsc_result[[1]]
# covariance.matrix.se <- ldsc_result[[3]]
# correlation.matrix <- ldsc_result[[2]]
# correlation.matrix.se <- ldsc_result[[4]]
# 
# 
# 
# 
# places <- 3
# 
# 
# 
# covariance.matrix.new <- rep("c",25)
# 
# covariance.matrix.transform <- log(covariance.matrix)
# covariance.matrix.transform.low <- covariance.matrix.transform - 1.96*covariance.matrix.se/covariance.matrix
# covariance.matrix.transform.high <- covariance.matrix.transform + 1.96*covariance.matrix.se/covariance.matrix
# covariance.matrix.high <- exp(covariance.matrix.transform.high)
# covariance.matrix.transform.low  <- exp(covariance.matrix.transform.low)
# covariance.matrix.high <- round(covariance.matrix.high,places)
# covariance.matrix.low <- round(covariance.matrix.low,places)
# 
# temp <- 1
# for(i in 1:5){
#   for(j in 1:5){
#     covariance.matrix.new[temp] <- paste0(covariance.matrix[i,j]," (",covariance.matrix.low[i,j],"-",covariance.matrix.high[i,j],")") 
#     temp <- temp+1
#   }
# }
# covariance.matrix.new <- matrix(covariance.matrix.new,5,5)
# 
# 
# 
# 
# 
# 
# correlation.matrix.transform <- 0.5*log((1+correlation.matrix)/(1-correlation.matrix))
# correlation.matrix.se.transform <- abs(correlation.matrix.se/((1+correlation.matrix)*(1-correlation.matrix)))
# 
# correlation.matrix.transform.low <- correlation.matrix.transform-1.96*correlation.matrix.se.transform
# correlation.matrix.transform.high <- correlation.matrix.transform+1.96*correlation.matrix.se.transform
# 
# correlation.matrix.low <- (exp(2*correlation.matrix.transform.low)-1)/(exp(2*correlation.matrix.transform.low)+1)
# correlation.matrix.high <- (exp(2*correlation.matrix.transform.high)-1)/(exp(2*correlation.matrix.transform.high)+1)
# 
# (exp(2*correlation.matrix.transform)-1)/(exp(2*correlation.matrix.transform)+1)
# 
# correlation.matrix <- round(correlation.matrix,places)
# correlation.matrix.low <- round(correlation.matrix.low,places)
# correlation.matrix.high <- round(correlation.matrix.high,places)
# 
# 
# correlation.matrix.new <- rep("c",25)
# 
# temp <- 1
# for(i in 1:5){
#   for(j in 1:5){
#     correlation.matrix.new[temp] <- paste0(correlation.matrix[i,j]," (",correlation.matrix.low[i,j],"-",correlation.matrix.high[i,j],")") 
#     temp <- temp+1
#   }
# }
# correlation.matrix.new <- matrix(correlation.matrix.new,5,5)
# 
# covariance.matrix.meta <- covariance.matrix.new
# 
# 
# 
# write.csv(covariance.matrix.meta,file="./covariance_matrix_meta.csv",quote=F)
# write.csv(correlation.matrix.new,file="./correlation_matrix_meta.csv",quote=F)
# 
# library(tidyverse)
# library(reshape2)
# library(ggplot2)
# library(gplots)
# pal.breaks <- seq(-1,1,0.01)
# col <- colorRampPalette(c("dodgerblue4","white","red"))(length(pal.breaks)-1)
# correlation.matrix.meta <- correlation.matrix
# # heatmap(correlation.matrix,col = col,symm=T,margins=c(10,4),key.title="",key.ylab="",cexRow=1,cexCol=1)
# png(filename="./meta_heatmap.png",width=800,heigh=600,units="px")
# heatmap.2(correlation.matrix.meta,tracecol=NA,cexRow=1,cexCol=1,margins=c(10,12),col = col,breaks=pal.breaks,key.ylab="",key.title = "",
#           main=" Genetic Correlation Heatmap Meta",dendrogram="row",density.info="none",lwid = c(1.5,4))
# dev.off()
# 
# p <- c(5,2,3,4,1)
# par(mfrow=c(1,3))
# png(filename="./meta_heatmap_original.png",width=800,heigh=600,units="px")
# heatmap.2(correlation.matrix.meta[p,p],tracecol=NA,cexRow=1,cexCol=1,margins=c(10,12),col = col,breaks=pal.breaks,key.ylab="",key.title = "",
#           main=" Genetic Correlation Heatmap Meta",dendrogram="none",density.info="none",lwid = c(1.5,4),
#           Rowv=FALSE,Colv=FALSE)
# dev.off()
# png(filename="./onco_heatmap_original.png",width=800,heigh=600,units="px")
# 
# heatmap.2(correlation.matrix.onco[p,p],tracecol=NA,cexRow=1,cexCol=1,margins=c(10,12),col = col,breaks=pal.breaks,key.ylab="",key.title = "",
#           main=" Genetic Correlation Heatmap Onco",dendrogram="none",density.info="none",lwid = c(1.5,4),
#           Rowv=FALSE,Colv=FALSE)
# dev.off()
# png(filename="./icog_heatmap_original.png",width=800,heigh=600,units="px")
# 
# heatmap.2(correlation.matrix.icog[p,p],tracecol=NA,cexRow=1,cexCol=1,margins=c(10,12),col = col,breaks=pal.breaks,key.ylab="",key.title = "",
#           main=" Genetic Correlation Heatmap ICOGs",dendrogram="none",density.info="none",lwid = c(1.5,4),
#           Rowv=FALSE,Colv=FALSE)
# 
# dev.off()
# 
# com <- combn(colnames(correlation.matrix.icog),2)
# co.heri.name <- paste0(com[1,]," ",com[2,])
# 
# covariance.matrix.icog[lower.tri(covariance.matrix.icog)]
# covariance.matrix.onco[lower.tri(covariance.matrix.icog)]
# covariance.matrix.meta[lower.tri(covariance.matrix.icog)]
# 
# 
# co.heri <- data.frame(co.heri.name,
#                       covariance.matrix.icog[lower.tri(covariance.matrix.icog)],
#                       covariance.matrix.onco[lower.tri(covariance.matrix.icog)],
#                       covariance.matrix.meta[lower.tri(covariance.matrix.icog)])
# 
# heri <- data.frame(colnames(correlation.matrix.icog),
#                    diag(covariance.matrix.icog),
#                    diag(covariance.matrix.onco),
#                    diag(covariance.matrix.meta))
# 
# write.csv(heri,file=paste0("./heritability.csv"))
# write.csv(co.heri,file=paste0("./coheritability.csv"))
