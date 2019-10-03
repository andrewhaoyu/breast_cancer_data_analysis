#load four tumor characteristics results with low effect size = 0.05, alpha = 1E-03
setwd('/Users/zhangh24/GoogleDrive/breast_cancer_data_analysis/simulation/power/result')
#data <- read.csv("power.simulation.result_0.25.csv")
data <- read.csv("power.simulation.result.csv")
data <- data[,-1,drop=F]
colnames(data) <- c("FTOP",
                    "MTOP",
                    "Stanadrd logistic regression",
                    "FTOP with only complete data",
                    "Polytomous model")

library(ggplot2)
Sample_size = rep(c(5000,25000,50000,100000),3)
Heter_pattern = c(rep("No Heterogeneity",4),
                  rep("One Marker Drove the Heterogeneity",4),
                  rep("Multiple Markers Drove the Heterogneity",4))
data = cbind(Sample_size,Heter_pattern,data)
data[,1] = as.character(data[,1])
Heter_pattern = c("No Heterogeneity",
                  "One Marker Drove the Heterogeneity",
                  "Multiple Markers Drove the Heterogneity")
p <- list()
for(i in 1:length(Heter_pattern)){
  Heter_pattern_temp = Heter_pattern[i]
  idx <- which(data$Heter_pattern==Heter_pattern_temp)  
  data.temp = data[idx,-2,drop=F]
  data.m = melt(data.temp,id="Sample_size")
  p[[i]] = ggplot(data.m,aes(x=Sample_size,y=value,fill=variable))+
    geom_bar(stat="identity",position=position_dodge(),
             alpha=0.8)+
    scale_fill_Publication()+
    theme_Publication()+
    scale_x_discrete(name = "Sample size", limits=c("5000","25000","50000",
                                                    "1e+05"))+
    scale_y_continuous(name="Power",limits = c(0,1))+
    ggtitle(Heter_pattern_temp)+
    theme(legend.position="none")
  
  }


library(reshape2)
library(gridExtra)
grid.arrange(p[[1]],p[[2]],p[[3]],nrow=1)
  
  
#load six tumor characteristics results with low effect size = 0.05, alpha = 1E-03
data <- read.csv("power_high.result.csv")

data <- data[,-1,drop=F]
colnames(data) <- c("MTOP",
                    "FTOP",
                    "Stanadrd logistic regression",
                    "FTOP with only complete data",
                    "Polytomous model")

library(ggplot2)
Sample_size = rep(c(5000,25000,50000,100000),3)
Heter_pattern = c(rep("No Heterogeneity",4),
                  rep("One Marker Drove the Heterogeneity",4),
                  rep("Multiple Markers Drove the Heterogneity",4))
data = cbind(Sample_size,Heter_pattern,data)
data[,1] = as.character(data[,1])
Heter_pattern = c("No Heterogeneity",
                  "One Marker Drove the Heterogeneity",
                  "Multiple Markers Drove the Heterogneity")
for(i in 1:length(Heter_pattern)){
  Heter_pattern_temp = Heter_pattern[i]
  idx <- which(data$Heter_pattern==Heter_pattern_temp)  
  data.temp = data[idx,-2,drop=F]
  data.m = melt(data.temp,id="Sample_size")
  #first three plots for four tumor markers
  #then 4-6 plots for six tumor markers
  p[[i+3]] = ggplot(data.m,aes(x=Sample_size,y=value,fill=variable))+
    geom_bar(stat="identity",position=position_dodge(),
             alpha=0.8)+
    scale_fill_Publication()+
    theme_Publication()+
    scale_x_discrete(name = "Sample size", limits=c("5000","25000","50000",
                                                    "1e+05"))+
    scale_y_continuous(name="Power",limits = c(0,1))+
    ggtitle(Heter_pattern_temp)+
    theme(legend.position="none")
  
}


library(reshape2)
library(gridExtra)
png(file = paste0("power_plot_0.05.png"),width=16,height=8,units = "in",res =300)
grid.arrange(p[[1]],p[[2]],p[[3]],
             p[[4]],p[[5]],p[[6]],
             nrow=2)
dev.off()
#generate one plot for legend
colnames(data.m)[2] <- "Method"
png(file = paste0("power_plot_lengend.png"),width=10,height = 8, units = "in", res = 300)
ggplot(data.m,aes(x=Sample_size,y=value,fill=Method))+
  geom_bar(stat="identity",position=position_dodge(),
           alpha=0.8)+
  scale_fill_Publication()+
  theme_Publication()+
  scale_x_discrete(name = "Sample size", limits=c("5000","25000","50000",
                                                  "1e+05"))+
  scale_y_continuous(name="Power",limits = c(0,1))+
  ggtitle(Heter_pattern_temp)
dev.off()
  
#plot(data.new[,3],data.new[,4])

#load four tumor characteristics results with low effect size = 0.05, alpha = 5E-08
setwd('/Users/zhangh24/GoogleDrive/breast_cancer_data_analysis/simulation/power/result')
data <- read.csv("power.simulation.result_5E8.csv")

data <- data[,-1,drop=F]
colnames(data) <- c("MTOP",
                    "FTOP",
                    "Stanadrd logistic regression",
                    "FTOP with only complete data",
                    "Polytomous model")

library(ggplot2)
Sample_size = rep(c(5000,25000,50000,100000),3)
Heter_pattern = c(rep("No Heterogeneity",4),
                  rep("One Marker Drove the Heterogeneity",4),
                  rep("Multiple Markers Drove the Heterogneity",4))
data = cbind(Sample_size,Heter_pattern,data)
data[,1] = as.character(data[,1])
Heter_pattern = c("No Heterogeneity",
                  "One Marker Drove the Heterogeneity",
                  "Multiple Markers Drove the Heterogneity")
p <- list()
for(i in 1:length(Heter_pattern)){
  Heter_pattern_temp = Heter_pattern[i]
  idx <- which(data$Heter_pattern==Heter_pattern_temp)  
  data.temp = data[idx,-2,drop=F]
  data.m = melt(data.temp,id="Sample_size")
  p[[i]] = ggplot(data.m,aes(x=Sample_size,y=value,fill=variable))+
    geom_bar(stat="identity",position=position_dodge(),
             alpha=0.8)+
    scale_fill_Publication()+
    theme_Publication()+
    scale_x_discrete(name = "Sample size", limits=c("5000","25000","50000",
                                                    "1e+05"))+
    scale_y_continuous(name="Power",limits = c(0,1))+
    ggtitle(Heter_pattern_temp)+
    theme(legend.position="none")
  
}


library(reshape2)
library(gridExtra)
grid.arrange(p[[1]],p[[2]],p[[3]],nrow=1)


#load six tumor characteristics results with low effect size = 0.05, alpha = 5E-08
data <- read.csv("power_high.result5E8.csv")

data <- data[,-1,drop=F]
colnames(data) <- c("MTOP",
                    "FTOP",
                    "Stanadrd logistic regression",
                    "FTOP with only complete data",
                    "Polytomous model")

library(ggplot2)
Sample_size = rep(c(5000,25000,50000,100000),3)
Heter_pattern = c(rep("No Heterogeneity",4),
                  rep("One Marker Drove the Heterogeneity",4),
                  rep("Multiple Markers Drove the Heterogneity",4))
data = cbind(Sample_size,Heter_pattern,data)
data[,1] = as.character(data[,1])
Heter_pattern = c("No Heterogeneity",
                  "One Marker Drove the Heterogeneity",
                  "Multiple Markers Drove the Heterogneity")
for(i in 1:length(Heter_pattern)){
  Heter_pattern_temp = Heter_pattern[i]
  idx <- which(data$Heter_pattern==Heter_pattern_temp)  
  data.temp = data[idx,-2,drop=F]
  data.m = melt(data.temp,id="Sample_size")
  #first three plots for four tumor markers
  #then 4-6 plots for six tumor markers
  p[[i+3]] = ggplot(data.m,aes(x=Sample_size,y=value,fill=variable))+
    geom_bar(stat="identity",position=position_dodge(),
             alpha=0.8)+
    scale_fill_Publication()+
    theme_Publication()+
    scale_x_discrete(name = "Sample size", limits=c("5000","25000","50000",
                                                    "1e+05"))+
    scale_y_continuous(name="Power",limits = c(0,1))+
    ggtitle(Heter_pattern_temp)+
    theme(legend.position="none")
  
}


library(reshape2)
library(gridExtra)
png(file = paste0("power_plot_0.05_5E8.png"),width=16,height=8,units = "in",res =300)
grid.arrange(p[[1]],p[[2]],p[[3]],
             p[[4]],p[[5]],p[[6]],
             nrow=2)
dev.off()
