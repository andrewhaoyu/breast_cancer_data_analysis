#load individual tumor characteristics test result size = 0.08, alpha = 8E-05
setwd('/Users/zhangh24/GoogleDrive/breast_cancer_data_analysis/simulation/power/result')
data <- read.csv("indi.simulation.result.csv")
#data <- read.csv("power.simulation.result.csv")
data <- data[,-1,drop=F]
colnames(data) <- c("FTOP with all markers",
"Polytomous model only uses ER")

library(ggplot2)
sizes <- c(25000,50000,100000)
Sample_size = rep(sizes,3)
n.sizes <- length(sizes)


Missing_rate = c(rep(0.17,3),
                rep(0.3,3),
                  rep(0.5,3))
data = cbind(Sample_size,Missing_rate,data)
data[,1] = as.character(data[,1])
p <- list()
missing_rate_vec <- c(0.17,0.30,0.50)
for(i in 1:length(missing_rate_vec)){
  missing_rate_temp = missing_rate_vec[i]
  idx <- which(data$Missing_rate==missing_rate_temp)  
  data.temp = data[idx,-2,drop=F]
  data.m = melt(data.temp,id="Sample_size")
  p[[i]] = ggplot(data.m,aes(x=Sample_size,y=value,fill=variable))+
    geom_bar(stat="identity",position=position_dodge(),
             alpha=0.8)+
    scale_fill_Publication()+
    theme_Publication()+
    scale_x_discrete(name = "Sample size", limits=c("25000","50000",
                                                    "1e+05"))+
    scale_y_continuous(name="Power",limits = c(0,0.35))+
    ggtitle(paste0 ("Missing rate of ER = ",missing_rate_temp))+
    theme(legend.position="none")
  
}


library(reshape2)
library(gridExtra)
png(paste0("plot_individual_test.png"),
    height = 8,width = 16,
    units = "in", res = 300)
grid.arrange(p[[1]],p[[2]],p[[3]],nrow=1)
dev.off()


colnames(data.m)[2] <- "Method"
data.m[,2] <- as.character(data.m[,2])
png(paste0("plot_individual_test_legend.png"),
    height = 8,width = 16,
    units = "in", res = 300)
ggplot(data.m,aes(x=Sample_size,y=value,fill=Method))+
  geom_bar(stat="identity",position=position_dodge(),
           alpha=0.8)+
  scale_fill_Publication()+
  theme_Publication()+
  scale_x_discrete(name = "Sample size", limits=c("25000","50000",
    "1e+05"))+
  scale_y_continuous(name="Power",limits = c(0,0.35))+
  ggtitle(paste0 ("Missing rate of ER = ",missing_rate_temp))+
  theme(legend.position="bottom")
dev.off()







setwd('/Users/zhangh24/GoogleDrive/breast_cancer_data_analysis/simulation/power/result')
data <- read.csv("indi.simulation.result_others.csv")

data <- data[,-1]

colnames(data) <- c(paste0("FTOP with all the markers ", c("PR","HER2","grade")),
                    
                    c(paste0("Polytomous model with only one marker ", c("PR","HER2","grade"))))


library(ggplot2)
sizes <- c(25000,50000,100000)
Sample_size = rep(sizes,3)
n.sizes <- length(sizes)


Missing_rate = c(rep(0.17,3),
                 rep(0.3,3),
                 rep(0.5,3))
data = cbind(Sample_size,Missing_rate,data)
data[,1] = as.character(data[,1])
p <- list()
missing_rate_vec <- c(0.17,0.30,0.50)
temp = 1
for(i in 1:length(missing_rate_vec)){
  missing_rate_temp = missing_rate_vec[i]
  idx <- which(data$Missing_rate==missing_rate_temp)  
  data.temp = data[idx,-2,drop=F]
  for(k in 2:2){
    
    data.m = melt(data.temp[,c(1,k+1,k+4)],id="Sample_size")
    p[[temp]] = ggplot(data.m,aes(x=Sample_size,y=value*1000,fill=variable))+
      geom_bar(stat="identity",position=position_dodge(),
               alpha=0.8)+
      scale_fill_Publication()+
      theme_Publication()+
      scale_x_discrete(name = "Sample size", limits=c("25000","50000",
                                                      "1e+05"))+
      scale_y_continuous(name="Type I error *1000",limits=c(0,105))+
      theme(legend.position="none")+
      geom_hline(yintercept = 1,color = "red",
                 linetype = 2)+
      ggtitle(paste0 ("Missing rate of ER = ",missing_rate_temp))
    temp <- temp+1
  }
  
}

  
  
  
  
for(i in 1:length(missing_rate_vec)){
  missing_rate_temp = missing_rate_vec[i]
  idx <- which(data$Missing_rate==missing_rate_temp)  
  for(k in 1:1){
    
    data.m = melt(data.temp[,c(1,k+1,k+4)],id="Sample_size")
    p[[temp]] = ggplot(data.m,aes(x=Sample_size,y=value*1000,fill=variable))+
      geom_bar(stat="identity",position=position_dodge(),
               alpha=0.8)+
      scale_fill_Publication()+
      theme_Publication()+
      scale_x_discrete(name = "Sample size", limits=c("25000","50000",
                                                      "1e+05"))+
      scale_y_continuous(name="Type I error *1000",limits=c(0,40))+
      theme(legend.position="none")+
      geom_hline(yintercept = 1,color = "red",
                 linetype = 2)+
      ggtitle(paste0 ("Missing rate of ER = ",missing_rate_temp))
    temp <- temp+1
  }
  
}

for(i in 1:length(missing_rate_vec)){
  missing_rate_temp = missing_rate_vec[i]
  idx <- which(data$Missing_rate==missing_rate_temp)  
  data.temp = data[idx,-2,drop=F]
  for(k in 3:3){
    
    data.m = melt(data.temp[,c(1,k+1,k+4)],id="Sample_size")
    p[[temp]] = ggplot(data.m,aes(x=Sample_size,y=value*1000,fill=variable))+
      geom_bar(stat="identity",position=position_dodge(),
               alpha=0.8)+
      scale_fill_Publication()+
      theme_Publication()+
      scale_x_discrete(name = "Sample size", limits=c("25000","50000",
                                                      "1e+05"))+
      scale_y_continuous(name="Type I error *1000",limits=c(0,5))+
      theme(legend.position="none")+
      geom_hline(yintercept = 1,color = "red",
                 linetype = 2)+
      ggtitle(paste0 ("Missing rate of ER = ",missing_rate_temp))
    temp <- temp+1
  }
  
}

library(reshape2)
library(gridExtra)
png(paste0("plot_individual_test_all_others.png"),
    height = 8,width = 16,
    units = "in", res = 300)
grid.arrange(p[[1]],p[[2]],p[[3]],
             p[[4]],p[[5]],p[[6]],
             p[[7]],p[[8]],p[[9]],
             nrow=3)
dev.off()

png(paste0("plot_individual_all_other_markers_lengend.png"),width=16,height=8,units="in",res = 300)
data.temp2 <- data.temp[,c(1,k+1,k+4)]
colnames(data.temp2)[c(2,3)] <- c("FTOP with all the markers","Polytomous model with only one marker")

data.m = melt(data.temp2,id="Sample_size")
colnames(data.m)[2] <- "Method"
ggplot(data.m,aes(x=Sample_size,y=value*1000,fill=Method))+
  geom_bar(stat="identity",position=position_dodge(),
           alpha=0.8)+
  scale_fill_Publication()+
  theme_Publication()+
  scale_x_discrete(name = "Sample size", limits=c("25000","50000",
                                                  "1e+05"))+
  scale_y_continuous(name="Type I error *1000",limits=c(0,5))+
  theme(legend.position="bottom")+
  geom_hline(yintercept = 1,color = "red",
             linetype = 2)+
  ggtitle(paste0 ("Missing rate of ER = ",missing_rate_temp))
dev.off()
