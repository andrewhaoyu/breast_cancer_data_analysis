#setwd('/spin1/users/zhangh24/breast_cancer_data_analysis/')
load(paste0("./discovery_SNP/sensitivity_analysis/result/sensitivity_result.Rdata"))
library(data.table)
data1 <- fread("./data/iCOGS_euro_v10_10232017.csv",header=T)
data2 <- fread("./data/Onco_euro_v10_10232017.csv",header=T)
all.countries <- unique(c(data1$StudyCountry,data2$StudyCountry))
level1 <- levels(as.factor(all.countries))
n.coun <- 20
all.countries[n.coun+1] <- "Summary"

discover.casecase <- read.csv("./data/discovery_snp_casecase.csv",header = T)
new.order <- c(1,9,15,13,2,6,5,8,17,4,10,18,
                 3,14,7,19)
discover.casecase <- discover.casecase[,]

library(forestplot)
plot.list <- list()

for(i1 in 1:length(new.order)){
  print(i1)
  data <- sensitivity.result[[new.order[i1]]]
  data <- rbind(data,discover.casecase[new.order[i1],2:ncol(discover.casecase)])
  all.countries <- factor(all.countries,levels=rev(c(level1,"Summary")))
  levels(all.countries)
  
  
  
  p <- c(1,4,7,10)
  odds = as.vector(as.matrix(data[,p]))
  low = c(as.matrix(data[,p+1]))
  high <- c(as.matrix(data[,p+2]))
  tumor <- c(rep("ER",n.coun+1),rep("PR",n.coun+1),
             rep("HER2",n.coun+1),rep("Grade",n.coun+1))
  tumor <- factor(tumor,levels = c("ER","PR","HER2","Grade"))
  Type <- c(rep("ind",n.coun),"sum")
  Type <- rep(Type,4)
  #Size <- c(rep(1,n.coun),1)
  new.data <- data.frame(odds,
                         low,
                         high,
                         Type,
                         label = all.countries,
                         tumor)
  
  library(ggplot2)
  fp <- ggplot(new.data,aes(x=label,y=odds,ymin=low,ymax=high,shape=Type,colour=Type))+
    geom_pointrange()+
    scale_colour_manual(values=c("#386cb0","#fdb462"))+
    #geom_line(yintercept=1,lty=2)+
    coord_flip()+
    theme_bw()+
    geom_hline(yintercept=1,size=1,lty=2)+
    theme(legend.position="none")+
    ylab("Odds Ratio (95% CI)")+
    #scale_y_continuous(limits=c(0.90,1.1))+
    xlab("Removed Countries")+
    facet_grid(.~tumor)+
    ggtitle(paste0("Forest plot of ", discover.casecase[new.order[i1],1])
    )+
    theme(plot.title = element_text(hjust=0.5))
  plot.list[[i1]] <- fp
  
  
}
i1 = 3
plot.list[[i1]]
library(gridExtra)
width = 36
height = 25
png("./discovery_SNP/sensitivity_analysis/result/sensitivity_plot1.png",width=width,
    height= height, units="cm", res= 600)
grid.arrange(plot.list[[1]],plot.list[[2]],plot.list[[3]],plot.list[[4]],plot.list[[5]],plot.list[[6]],nrow=3,ncol=2)
dev.off()
png("./discovery_SNP/sensitivity_analysis/result/sensitivity_plot2.png",width=width,
    height= height, units="cm", res= 600)
grid.arrange(plot.list[[7]],plot.list[[8]],plot.list[[9]],plot.list[[10]],plot.list[[11]],plot.list[[12]],nrow=3,ncol=2)
dev.off()
png("./discovery_SNP/sensitivity_analysis/result/sensitivity_plot3.png",width=width,
    height= (height*2/3), units="cm", res= 600)
grid.arrange(plot.list[[13]],plot.list[[14]],plot.list[[15]],plot.list[[16]],nrow=2,ncol=2)
dev.off()


dev.off()
png("./discovery_SNP/sensitivity_analysis/result/sensitivity_plot6.png",width=width,
    height= height/3, units="cm", res= 600)
plot.list[[16]]
dev.off()
# png("./discovery_SNP/sensitivity_analysis/result/sensitivity_plot7.png",width=(width),
#     height= (height/3), units="cm", res= 600)
# plot.list[[19]]
# dev.off()
