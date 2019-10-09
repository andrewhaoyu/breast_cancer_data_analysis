#load four tumor characteristics results computing time information
setwd('/Users/zhangh24/GoogleDrive/breast_cancer_data_analysis/simulation/power/result')
data <- read.csv("time.simulation.result.csv")
#data <- read.csv("power.simulation.result.csv")
data <- data[,-1,drop=F]
colnames(data) <- c("FTOP",
                    "MTOP",
                    "Standard logistic regression",
                    "Two-stage model with only complete data",
                    "Polytomous model")

library(ggplot2)
sizes <- c(5000,25000,50000,100000)
Sample_size = rep(sizes,1)
n.sizes <- length(sizes)


data = cbind(Sample_size,data)
data[,1] = as.character(data[,1])

 
Method1 <- c("FTOP","Standard logistic regression",
             "Two-stage model with only complete data")
Method2 <- c("MTOP","Polytomous model")

Method <- list(Method1,Method2)
  
  data.m = melt(data,id="Sample_size")
   colnames(data.m)[2] <- "Method"
p <- list()
   k <- 1
  idx <- which(data.m$Method%in%Method[[k]])
  data.m.temp <- data.m[idx,]
  p[[k]] <- ggplot(data.m.temp,aes(x=Sample_size,y=value,fill=Method))+
    geom_bar(stat="identity",position=position_dodge(),
             alpha=0.8)+
    theme_Publication()+
    scale_x_discrete(name = "Sample size", limits=c("5000","25000","50000",
                                                    "1e+05"))+
    scale_y_continuous(name="Computing time (s)")+
    ggtitle(Heter_pattern_temp)+
    theme(legend.position="none")+
    discrete_scale("fill","Publication",manual_pal(values = c("#386cb0","#FFD042","#7fc97f","#004953","#662506","#a6cee3","#fb9a99","#984ea3","#ffff33")))+
    ggtitle(NULL)
  
    
  
k <- 2
idx <- which(data.m$Method%in%Method[[k]])
data.m.temp <- data.m[idx,]
p[[k]] <- ggplot(data.m.temp,aes(x=Sample_size,y=value,fill=Method))+
  geom_bar(stat="identity",position=position_dodge(),
           alpha=0.8)+
  theme_Publication()+
  scale_x_discrete(name = "Sample size", limits=c("5000","25000","50000",
                                                  "1e+05"))+
  scale_y_continuous(name="Computing time (s)")+
  theme(legend.position="none")+
  discrete_scale("fill","Publication",manual_pal(values = c("#EF7E3D","#004953","#662506","#a6cee3","#fb9a99","#984ea3","#ffff33")))+
  ggtitle(NULL)
png(file = "Computation_time.png",width = 16, height = 8, units = "in", res = 300)
grid.arrange(p[[1]],p[[2]],nrow=1
             )
 dev.off()



library(reshape2)
library(gridExtra)
grid.arrange(p[[1]],p[[2]],p[[3]],nrow=1)


#load six tumor characteristics results with low effect size = 0.05, alpha = 1E-03
data <- read.csv("power_high.result.csv")

data <- data[,-1,drop=F]
colnames(data) <- c("FTOP",
                    "MTOP",
                    "Stanadrd logistic regression",
                    "Two-stage model with only complete data",
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
idx <- which(data$Sample_size==5000)
data <- data[-idx,]

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
    scale_x_discrete(name = "Sample size", limits=c("25000","50000",
                                                    "1e+05"))+
    scale_y_continuous(name="Power",limits = c(0,1))+
    ggtitle(Heter_pattern_temp)+
    theme(legend.position="none")
  
}


library(reshape2)
library(gridExtra)
png(file = paste0("power_plot_0.08.png"),width=16,height=8,units = "in",res =300)
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
  ggtitle(Heter_pattern_temp)+
  theme(legend.position = "bottom")
dev.off()

#plot(data.new[,3],data.new[,4])





#plot for high effect size 0.25
setwd('/Users/zhangh24/GoogleDrive/breast_cancer_data_analysis/simulation/power/result')
data <- read.csv("power.simulation.result_0.25.csv")
#data <- read.csv("power.simulation.result.csv")
data <- data[,-1,drop=F]
colnames(data) <- c("FTOP",
                    "MTOP",
                    "Standard logistic regression",
                    "Two-stage model with only complete data",
                    "Polytomous model")

library(ggplot2)
sizes = c(5000)
Sample_size = rep(sizes,3)
n.sizes <- length(sizes)
Heter_pattern = c(rep("No Heterogeneity",n.sizes),
                  rep("One Marker Drove the Heterogeneity",n.sizes),
                  rep("Multiple Markers Drove the Heterogneity",n.sizes))
data = cbind(Sample_size,Heter_pattern,data)
data[,1] = as.character(data[,1])
Heter_pattern = c("No Heterogeneity",
                  "One Marker Drove the Heterogeneity",
                  "Multiple Markers Drove the Heterogneity")
library(reshape2)
p <- list()
data <- data[,-1]
data.m <- melt(data,id = "Heter_pattern")
colnames(data.m)[2] <- "Method"
p[[1]]=ggplot(data.m,aes(x=Heter_pattern,y=value,fill=Method))+
  geom_bar(stat="identity",position=position_dodge())+
  theme_Publication()+
  scale_fill_Publication()+
  scale_x_discrete(name="Heterogenous pattern",limits=c("No Heterogeneity",
                                                        "One Marker Drove the Heterogeneity",
                                                        "Multiple Markers Drove the Heterogneity"))+
  scale_y_continuous(name="Power",limits=c(0,1))+
  theme(axis.text.x = element_text(angle=0,vjust =0.5))+
  theme(legend.position = "none")



#load six tumor characteristics results with low effect size = 0.25, alpha = 5E-08
data <- read.csv("power_high_0.25_result.csv")

data <- data[,-1,drop=F]
colnames(data) <- c("MTOP",
                    "FTOP",
                    "Stanadrd logistic regression",
                    "Two-stage model with only complete data",
                    "Polytomous model")

library(ggplot2)
sizes = c(5000)
Sample_size = rep(sizes,3)
n.sizes <- length(sizes)
Heter_pattern = c(rep("No Heterogeneity",n.sizes),
                  rep("One Marker Drove the Heterogeneity",n.sizes),
                  rep("Multiple Markers Drove the Heterogneity",n.sizes))
data = cbind(Sample_size,Heter_pattern,data)
data <- data[,-1]
data.m <- melt(data,id = "Heter_pattern")
colnames(data.m)[2] = "Method"
p[[2]] = ggplot(data.m,aes(x= Heter_pattern,y=value,fill=Method))+
  geom_bar(stat="identity",position = position_dodge())+
  theme_Publication()+
  scale_fill_Publication()+
  scale_x_discrete(name = "Heterogeneity pattern",limits = c("No Heterogeneity",
                                                             "One Marker Drove the Heterogeneity",
                                                             "Multiple Markers Drove the Heterogneity"))+
  scale_y_continuous(name="Power",limits=c(0,1))+
  theme(axis.text.x = element_text(angle= 0,vjust = 0.5))+
  theme(legend.position = "none")
library(reshape2)
library(gridExtra)
png(file = paste0("power_plot_0.25.png"),width=16,height=8,units = "in",res =300)
grid.arrange(p[[1]],p[[2]],
             nrow=2)
dev.off()
