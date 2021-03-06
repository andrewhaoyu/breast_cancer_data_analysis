setwd("/Users/zhangh24/GoogleDrive/breast_cancer_data_analysis")
odds.ratio <- read.csv("./risk_prediction/Nasim_new_prs/result/odds.ratio.csv")
odds.ratio.overall <- read.csv("./risk_prediction/Nasim_prs/result/odds.ratio.overall.csv")
odds.ratio.er <- read.csv("./risk_prediction/Nasim_prs/result/odds.ratio.er.csv")
quantile.thres <- factor(c("<1%","1%-5%","5%-10%",
                           "10%-20%","20%-40%","40%-60%","60%-80%","80%-90%","90%-95%","95%-99%",">99%"))

quantile.thres <- factor(quantile.thres,levels(quantile.thres)[c(1,3,7,4,5,6,8,9,10,11,2)])

n <- nrow(odds.ratio)
odds.ratio <- odds.ratio[,-1,drop=F]
method <- rep("Intrinsic subtyeps PRS")
odds.ratio <- cbind(quantile.thres,method,odds.ratio)

odds.ratio.overall <- odds.ratio.overall[,-1,drop=F]
method <- rep("Overall PRS",n)
odds.ratio.overall <- cbind(quantile.thres,method,odds.ratio.overall)
odds.ratio.er <- odds.ratio.er[,-1,drop=F]
method <- rep("ER specific PRS",n)
odds.ratio.er <- cbind(quantile.thres,method,odds.ratio.er)



select.names <- c("Luminal A-like ",
    "Luminal B-like",  
                       "Luminal B/HER2-negative-like",
                       " HER2-enriched-like",
                       "TN")
names.subtypes <-  c("Luminal_A",
                      "Luminal_B",
                      "Luminal_B_HER2Neg",
                      "HER2_Enriched",
                      "TN")



data <- rbind(odds.ratio,odds.ratio.overall,
              odds.ratio.er)
scale_colour_Publication <- function(...){
  library(scales)
  discrete_scale("colour","Publication",manual_pal(values = c("#386cb0","#fb9a99","chartreuse4","#FFD042","#7fc97f","#004953","#662506","#a6cee3","#984ea3","#EF7E3D","#c0392b")), ...)
  
}

library(ggplot2)
p <- list()
for(k in 1:length(select.names)){
  data.temp <- data[,c(1:2,k+2)]
  colnames(data.temp) <- c("quantile","method","or")
print({
  png(file = paste0("./risk_prediction/Nasim_new_prs/result/",names.subtypes[k],".png"),width=10,height =6,unit = "in",res=300)  
  ggplot(data.temp,aes(x=quantile,
                       y= or,group=method,color=method))+
    geom_line(size=0.9)+
    theme_Publication()+
    scale_colour_Publication()+
    xlab(NULL)+
    ylab("Odds ratio")+
    ggtitle(paste0(select.names[k]))+
    scale_y_continuous(limits=(c(0,6)),
                       breaks= c(0,1,2,4,6))+
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
}
  )
dev.off()
}

