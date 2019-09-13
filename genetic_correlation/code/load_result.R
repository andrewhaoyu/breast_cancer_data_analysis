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


get_cor_se <- function(correlation.matrix,correlation.matrix.se){

  correlation.matrix <- round(correlation.matrix,places)
  correlation.matrix.se <- round(correlation.matrix,places)
  M <- nrow(correlation.matrix)
  
  correlation.matrix.new <- rep("c",M^2)
  
  temp <- 1
  for(i in 1:M){
    for(j in 1:M){
      correlation.matrix.new[temp] <- paste0(correlation.matrix[i,j]," (",correlation.matrix.se[i,j],")") 
      temp <- temp+1
    }
  }
  correlation.matrix.new <- matrix(correlation.matrix.new,M,M)
  return(correlation.matrix.new)  
}









get_cor_95_list <- function(correlation.matrix,correlation.matrix.se){
  correlation.matrix.transform <- 0.5*log((1+correlation.matrix)/(1-correlation.matrix))
  correlation.matrix.se.transform <- abs(correlation.matrix.se/((1+correlation.matrix)*(1-correlation.matrix)))
  
  correlation.matrix.transform.low <- correlation.matrix.transform-1.96*correlation.matrix.se.transform
  correlation.matrix.transform.high <- correlation.matrix.transform+1.96*correlation.matrix.se.transform
  
  correlation.matrix.low <- (exp(2*correlation.matrix.transform.low)-1)/(exp(2*correlation.matrix.transform.low)+1)
  correlation.matrix.high <- (exp(2*correlation.matrix.transform.high)-1)/(exp(2*correlation.matrix.transform.high)+1)
  diag(correlation.matrix.se.transform) = 0.1
  diag(correlation.matrix.transform) = 1000
  correlation.matrix.se.transform = as.matrix(correlation.matrix.se.transform)
  correlation.matrix.transform = as.matrix(correlation.matrix.transform)
  p.value = 2*pnorm(abs(correlation.matrix.transform/correlation.matrix.se.transform),lower.tail = F)
  
 
  return(list(correlation.matrix.low, correlation.matrix.high,p.value))
}



setwd('/Users/zhangh24/GoogleDrive/breast_cancer_data_analysis')
# load("./breast_cancer_ldsc/ICOG_ldsc_result.rda")
# load("./breast_cancer_ldsc/ONCO_ldsc_result.rda")
#load("./breast_cancer_ldsc/two-stage-model-intrinsic_subtype/meta_ldsc_result.rda")
#
#
#
##################only based on BCAC
places <- 2
load("./genetic_correlation/result/ldsc_result_meta_091219.rda")

diag(ldsc_result[[1]][c(2,5,4,3,1,6),c(2,5,4,3,1,6)])
diag(ldsc_result[[3]][c(2,5,4,3,1,6),c(2,5,4,3,1,6)])
write.csv(ldsc_result[[1]],file = "./data/BCAC_heritability.csv")

correlation.matrix <- ldsc_result[[2]]
correlation.matrix.se <- ldsc_result[[4]]
correlation.matrix <- correlation.matrix[c(2,5,4,3,1,6),c(2,5,4,3,1,6)]
correlation.matrix.se <- correlation.matrix.se[c(2,5,4,3,1,6),c(2,5,4,3,1,6)]

names <- colnames(ldsc_result[[2]])[c(2,5,4,3,1,6)]
correlation.95 <- get_cor_95(correlation.matrix,correlation.matrix.se)
correlation.se <- get_cor_se(correlation.matrix,correlation.matrix.se)
names <- colnames(ldsc_result[[2]])[c(2,5,4,3,1,6)]
colnames(correlation.matrix) <- names
rownames(correlation.matrix) <- names

colnames(correlation.95) <- names
rownames(correlation.95) <- names
colnames(correlation.se) <- names
rownames(correlation.se) <- names

correlation.matrix.low <- get_cor_95_list(correlation.matrix,correlation.matrix.se)[[1]]
correlation.matrix.high<- get_cor_95_list(correlation.matrix,correlation.matrix.se)[[2]]
write.csv(correlation.95,file="./genetic_correlation/result/correlation_95_BCAC_CIMBA_091219.csv",quote=F)
write.csv(correlation.se,file="./genetic_correlation/result/correlation_se_BCAC_CIMBA_091219.csv",quote=F)


# #get results for CIMBA+ BCAC
#  setwd("/Users/zhangh24/GoogleDrive/breast_cancer_data_analysis/genetic_correlation/result/")
#  intrinsic_result <- read.csv("/Users/zhangh24/GoogleDrive/breast_cancer_data_analysis/genetic_correlation/result/CIMBA_BCAC.csv",header=F)
#  correlation.matrix <- intrinsic_result[1:5,2:6]
#  correlation.matrix.se <- intrinsic_result[6:10,2:6]
#  correlation.matrix <- correlation.matrix[c(2,5,4,3,1),c(2,5,4,3,1)]
#  correlation.matrix.se <- correlation.matrix.se[c(2,5,4,3,1),c(2,5,4,3,1)]
#  correlation.95 <- get_cor_95(correlation.matrix,correlation.matrix.se)
#  names <- as.character(intrinsic_result[c(2,5,4,3,1),1])
#  colnames(correlation.matrix) <- names
#  rownames(correlation.matrix) <- names
# 
#  colnames(correlation.95) <- names
#  rownames(correlation.95) <- names
# 
#  correlation.matrix.low <- get_cor_95_list(correlation.matrix,correlation.matrix.se)[[1]]
#  correlation.matrix.high<- get_cor_95_list(correlation.matrix,correlation.matrix.se)[[2]]


#
# write.csv(correlation.95,file="correlation_95_CIMBA_BCAC.csv",quote=F)
# places <- 3




p.value<- get_cor_95_list(correlation.matrix,correlation.matrix.se)[[3]]

M <- nrow(correlation.matrix)
combn.list <- combn(M,2)
n <- ncol(combn.list)
cor.vec.high <- cor.vec.low <- cor.vec <- rep(0,n)
subtypes <- rep("c",n)
names <- c("Luminal A-like",
           "Luminal B, HER2-negative-like",
           "Luminal B-like ",
           "HER2-enriched-like ",
           "TN",
           "CIMBA BRCA1")
for(i in n:1){
  subtypes[i] <- paste0(names[combn.list[1,i]],
                     " vs ",
                     names[combn.list[2,i]])
  cor.vec[i] <- correlation.matrix[combn.list[2,i],combn.list[1,i]]
  cor.vec.low[i] <- correlation.matrix.low[combn.list[2,i],combn.list[1,i]]
  cor.vec.high[i] <- correlation.matrix.high[combn.list[2,i],combn.list[1,i]]
}
subtypes.f <- factor(subtypes,levels=subtypes)

cor.data <-data.frame(subtypes.f,cor.vec,
                      cor.vec.low,
                      cor.vec.high)

cor.data <- cor.data[order(cor.data[,2]),]
cor.data[,1] <- factor(cor.data[,1],
                       levels=as.character(cor.data[,1]))
library(ggplot2)

png(filename="./genetic_correlation/result/genetic_correlation_plot_ci.png",width=13,heigh=9.5,units="in",res=600)
ggplot(cor.data,aes(x=subtypes.f,y=cor.vec))+
  geom_point(size=4)+
  geom_errorbar(aes(ymax = cor.vec.high, ymin = cor.vec.low))+
  theme_minimal()+  theme(axis.text.x = element_text(angle = 90, hjust = 1)) + ylab("Genetic correlation") + 
  ggtitle("Genetic correlation between five intrinsic subtypes")  + 
  theme(text = element_text(size=18),plot.title = element_text(hjust = 0.5,face = "bold"),axis.text.x = element_text(face = "bold"),
        axis.text.y = element_text(face = "bold.italic"),
        axis.title.x = element_text(face = "bold"),
        axis.title.y = element_text(face = "bold"))+
  xlab("Intrinsic subtypes")+
  scale_y_continuous(limits=c(0,1))+
  coord_flip()
dev.off()
  #geom_hline (yintercept = -log10(0.05/220), color = "red")+
  
  # coord_flip()+
  # facet_grid(.~subtypes)+
  # theme(strip.text = element_text(face = "bold"))






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

library(dplyr)
library(reshape2)
library(ggplot2)
#library(gplots)



# correlation.matrix <- as.matrix(as.data.frame(fread("genetic_correlation_meta.csv")[,1:5]))
#rownames(correlation.matrix) <- colnames(correlation.matrix)
# correlation.matrix <- correlation.matrix[c(2,5,4,3,1),c(2,5,4,3,1)]
# rownames(correlation.matrix) <- colnames(correlation.matrix)
# correlation.matrix <- as.matrix(correlation.matrix)
#correlation.matrix <- correlation.matrix[c(2,4,3,1),c(2,4,3,1)]

library(corrplot)




pal.breaks <- seq(0.5,1,0.01)
col <- colorRampPalette(c("white","red"))(length(pal.breaks)-1)

# heatmap(correlation.matrix,col = col,symm=T,margins=c(10,4),key.title="",key.ylab="",cexRow=1,cexCol=1)
paletteLength <- 50
load("./genetic_correlation/result/ldsc_result_meta_091219.rda")

correlation.matrix <- ldsc_result[[2]]
correlation.matrix.se <- ldsc_result[[4]]
correlation.matrix <- correlation.matrix[c(2,5,4,3,1,6),c(2,5,4,3,1,6)]
correlation.matrix.se <- correlation.matrix.se[c(2,5,4,3,1,6),c(2,5,4,3,1,6)]

my.color <- colorRampPalette(c("white", "dodgerblue4"))(paletteLength)
myBreaks <- c(seq(0, max(correlation.matrix), length.out=floor(paletteLength)))

png(filename="./genetic_correlation/result/meta_heatmap_two_stage.png",width=10,heigh=10,units="in",res=600)
colnames(correlation.matrix) <- c("Luminal A-like",
                                  "Luminal B, HER2-negative-like",
                                  "Luminal B-like ",
                                  "HER2-enriched-like ",
                                  "TN",
                                  "CIMBA BRCA1")
rownames(correlation.matrix) <- c("Luminal A-like",
                                  "Luminal B, HER2-negative-like",
                                  "Luminal B-like ",
                                  "HER2-enriched-like ",
                                  "TN",
                                  "CIMBA BRCA1")
#corrplot.mixed(as.matrix(correlation.matrix),
             #  lower.col = "black", number.cex = .7)
# correlation.matrix <- correlation.matrix[c(3,2,1,4,5,6),c(3,2,1,4,5,6)]
# correlation.matrix.se <- correlation.matrix.se[c(3,2,1,4,5,6),c(3,2,1,4,5,6)]


M <- 6
correlation.result <- matrix(rep("c",M^2),M,M)
for(i in 1:M){
  for(j in 1:M){
    correlation.result[i,j] <- paste0(
      round(correlation.matrix[i,j],2)," (",
      round(correlation.matrix.se[i,j],2),")")
  }
}
colnames(correlation.result) <- colnames(correlation.matrix) 
rownames(correlation.result) <- rownames(correlation.matrix)
  
  
write.csv(correlation.result,file="./genetic_correlation/result/plot_order_BCAC_CIMBA_gc.csv")
corrplot(as.matrix(correlation.matrix),
         method = "circle",
         #p.mat = p.value,
         #insig = "label_sig",
         #low=as.matrix(correlation.matrix.low),
         #upp= as.matrix(correlation.matrix.high),
         
         #rect.col = "navy", plotC = "rect", cl.pos = "n",
         #col = my.color,
         # method = "color",
         tl.cex=1.2,
         addrect=2, 
         order = "hclust",
         tl.col = "black", 
         tl.srt = 30,
         cl.lim = c(0, 1),
         #sig.level = c(0.005,0.01,0.05),
         pch.cex = 0.9, pch.col = "white"
)
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
