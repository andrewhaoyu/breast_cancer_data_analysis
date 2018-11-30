#-------------------------------------------------------------------
# Update Date: 11/24/2018
# Create Date: 11/24/2018
# Goal: analyze auc
# Author: Haoyu Zhang
#-------------------------------------------------------------------
#---------------------------------------#---------------------------------------
load("./risk_prediction/LD_clumping/result/auc.result.rdata")
head(auc.result)

auc.result %>% filter(subtypes=="Luminal_B")
names.subtypes <- c("Luminal_A",
                    "Luminal_B",
                    "Luminal_B_HER2Neg",
                    "HER2Enriched",
                    "TripleNeg")
method.names <- c("standard","two-stage","eb")

total <- length(names.subtypes)*length(method.names)
method.out <- rep("c",total)
subtypes.out <- rep("c",total)
auc.out <- rep(0,total)
ind <- 1
for(i in 1:length(names.subtypes)){
  for(j in 1:length(method.names)){
    auc.out[ind] <- auc.result %>% filter(subtypes==names.subtypes[i]&
                            method==method.names[j]) %>%
      select(auc) %>% 
      max
    method.out[ind] <- method.names[j]
    subtypes.out[ind] <- names.subtypes[i]
    ind <- ind + 1
  }
}
result.temp <- data.frame(method.out,subtypes.out,auc.out)
colnames(result.temp) <- c("method","subtypes","auc")
n.method <- 3
n.subtypes <- 5

auc.table <- matrix(result.temp[,3],n.subtypes,n.method,byrow=T)
row.names(auc.table) <- names.subtypes
colnames(auc.table) <- method.names
write.csv(auc.table,"./risk_prediction/LD_clumping/result/auc.table.csv")
source("./plot_support/publish_theme.R")
plot.list <- list()
temp <- 1
library(dplyr)
library(ggplot2)
library(gridExtra)
temp <- 1
for(i in 1:length(names.subtypes)){
  for(j in 1:length(method.names)){
    auc.temp <- auc.result %>% filter(method==method.names[j]&
                            subtypes==names.subtypes[i])
    p.temp <- ggplot(auc.temp,aes(as.factor(p),auc))+geom_point(size=3)+scale_colour_Publication()+theme_Publication()+
      xlab("p-value")+
      ylab("auc")+
      ggtitle(paste0(names.subtypes[i],"_",method.names[j]))+
    theme(axis.text.x = element_text(angle = 75, hjust = 1))
    plot.list[[temp]] <- p.temp
    temp <- temp+1
  }
}

temp <- 1
for(i in 1:length(names.subtypes)){
  png(filename = paste0("./risk_prediction/LD_clumping/result/",names.subtypes[i],".png"),
      width=40,height=20,res=300,units="cm")
  grid.arrange(plot.list[[temp]],plot.list[[temp+1]],plot.list[[temp+2]],
               nrow=1)
  dev.off()
  temp <- temp+3
  
}
