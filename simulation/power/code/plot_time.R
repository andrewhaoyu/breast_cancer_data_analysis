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

 
Method1 <- c("FTOP","MTOP","Standard logistic regression")
Method2 <- c("Two-stage model with only complete data","Polytomous model")

Method <- list(Method1,Method2)
  
  data.m = melt(data,id="Sample_size")
   colnames(data.m)[2] <- "Method"
p <- list()
   
k <- 1
idx <- which(data.m$Method%in%Method[[k]])
data.m.temp <- data.m[idx,]
p[[k]] <- ggplot(data.m.temp,aes(x=Sample_size,y=log(1000*value),fill=Method))+
  geom_bar(stat="identity",position=position_dodge(),
           alpha=0.8)+
  theme_Publication()+
  scale_x_discrete(name = "Sample size", limits=c("5000","25000","50000",
                                                  "1e+05"))+
  scale_y_continuous(name="Log (computing time)",limits = c(0,14.5))+
  theme(legend.position="none")+
  discrete_scale("fill","Publication",manual_pal(values = c("#386cb0","#EF7E3D","#FFD042","#7fc97f","#004953","#662506","#a6cee3","#fb9a99","#984ea3","#ffff33")))+
  ggtitle("Methods working on all the data")
  
k <- 2
idx <- which(data.m$Method%in%Method[[k]])
data.m.temp <- data.m[idx,]
p[[k]] <- ggplot(data.m.temp,aes(x=Sample_size,y=log(1000*value),fill=Method))+
  geom_bar(stat="identity",position=position_dodge(),
           alpha=0.8)+
  theme_Publication()+
  scale_x_discrete(name = "Sample size", limits=c("5000","25000","50000",
                                                  "1e+05"))+
  scale_y_continuous(name="Log (computing time)",limits = c(0,14.5))+
  theme(legend.position="none")+
  discrete_scale("fill","Publication",manual_pal(values = c("#7fc97f","#004953","#662506","#a6cee3","#fb9a99","#984ea3","#ffff33")))+
  ggtitle("Methods only working on complete data")

data <- read.csv("time.simulation.result_high.csv")
#data <- read.csv("power.simulation.result.csv")
data <- data[,-1,drop=F]
colnames(data) <- c("FTOP",
                    "MTOP",
                    "Standard logistic regression",
                    "Two-stage model with only complete data",
                    "Polytomous model")


sizes <- c(5000,25000,50000,100000)
Sample_size = rep(sizes,1)
n.sizes <- length(sizes)


data = cbind(Sample_size,data)
data[,1] = as.character(data[,1])


Method1 <- c("FTOP","MTOP","Standard logistic regression")
Method2 <- c("Two-stage model with only complete data","Polytomous model")

Method <- list(Method1,Method2)

data.m = melt(data,id="Sample_size")
colnames(data.m)[2] <- "Method"


k <- 3
idx <- which(data.m$Method%in%Method[[k-2]])
data.m.temp <- data.m[idx,]
p[[k]] <- ggplot(data.m.temp,aes(x=Sample_size,y=log(1000*value),fill=Method))+
  geom_bar(stat="identity",position=position_dodge(),
           alpha=0.8)+
  theme_Publication()+
  scale_x_discrete(name = "Sample size", limits=c("5000","25000","50000",
                                                  "1e+05"))+
  scale_y_continuous(name="Log (computing time)",limits = c(0,14.5))+
  theme(legend.position="none")+
  discrete_scale("fill","Publication",manual_pal(values = c("#386cb0","#EF7E3D","#FFD042","#7fc97f","#004953","#662506","#a6cee3","#fb9a99","#984ea3","#ffff33")))+
  ggtitle("Methods working on all the data")

k <- 4
idx <- which(data.m$Method%in%Method[[k-2]])
data.m.temp <- data.m[idx,]
p[[k]] <- ggplot(data.m.temp,aes(x=Sample_size,y=log(1000*value),fill=Method))+
  geom_bar(stat="identity",position=position_dodge(),
           alpha=0.8)+
  theme_Publication()+
  scale_x_discrete(name = "Sample size", limits=c("5000","25000","50000",
                                                  "1e+05"))+
  scale_y_continuous(name="Log (computing time)",limits = c(0,14.5))+
  theme(legend.position="none")+
  discrete_scale("fill","Publication",manual_pal(values = c("#7fc97f","#004953","#662506","#a6cee3","#fb9a99","#984ea3","#ffff33")))+
  ggtitle("Methods only working on complete data")




png(file = "Computation_time.png",width = 16, height = 8, units = "in", res = 300)
grid.arrange(p[[1]],p[[2]],
             p[[3]],p[[4]],nrow=2
)
dev.off()

