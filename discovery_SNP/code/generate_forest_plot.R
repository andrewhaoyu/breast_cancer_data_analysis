#goal generate forest plot for discovery SNPs
library(data.table)
setwd('/Users/zhangh24/GoogleDrive/breast_cancer_data_analysis')
library(data.table)
data <- as.data.frame(fread("./discovery_SNP/result/discovery_snp_analysis.csv"))
n.sub <- nrow(data)
n.var <- c(rep("ER",n.sub),
           rep("PR",n.sub),
           rep("HER2",n.sub),
           rep("Grade",n.sub),
           rep("Luminal A-like",n.sub),
           rep("Luminal B,HER2-negative-like",n.sub),
           rep("Luminal B-like",n.sub),
             rep("HER2 enriched-like",n.sub),
             rep("TN",n.sub),
           rep("BRCA1 carriers",n.sub)
           )
snps <- rep(data$rs_id,n.cat)
n.cat <- 10
OR <- rep(0,length(n.var))
low <- rep(0,length(n.var))
high <- rep(0,length(n.var))
temp <- 1

  for(j in c(7,9,11,13,21,23,25,27,29,31)){
    for(i in 1:n.sub){
      tempstr <- data[i,j]
    temp.split <- strsplit(tempstr,"\\(")
    OR[temp] <- as.numeric(temp.split[[1]][1])
    temp.split2 <- strsplit(temp.split[[1]][2],"-")
    low[temp] <-as.numeric(temp.split2[[1]][1])
    high[temp] <- as.numeric(strsplit(temp.split2[[1]][2],")")[[1]][1])
    temp <- temp+1
  }
}
library(dplyr)
#######prepare data for forest plot
data.plot <- data.frame(snps,n.var,OR,low,high,stringsAsFactors = F)
colnames(data.plot) <- c("SNP","Var","OR","low","high")
case.char <- c("ER","PR","HER2","Grade")
casecon.car <- c("Luminal A-like","Luminal B,HER2-negative-like",
                 "Luminal B-like",
                 "HER2 enriched-like",
                 "TN",
                 "BRCA1 carriers" )
########plot intrinsic subtyeps and case-case plor for each snp
for(i in 1:n.sub){
  snp.temp <- data$rs_id[i]
  png(paste0("./discovery_SNP/result/forest_plot/casecase_",snp.temp,".png"),height=15,width = 17,res=300,units="cm")
  data.plot.snp1 = data.plot%>% filter(SNP==snp.temp) %>% filter(Var%in%case.char)
  temp.case.char <- factor(data.plot.snp1$Var)
  data.plot.snp1$Var <- factor(temp.case.char,levels(temp.case.char)[c(2,3,4,1)])
  print(
    ggplot(data.plot.snp1,aes(x=Var,y=OR,ymin=low,ymax=high,colour=Var))+
      geom_pointrange()+
      #scale_colour_manual(values=c("#386cb0","#fdb462"))+
      #geom_line(yintercept=1,lty=2)+
      coord_flip()+
      theme_bw()+
      geom_hline(yintercept=1,size=1,lty=2)+
      theme(legend.position="none")+
      #ylab("Case-case odds ratio")+
      #xlab("Tumor characteristics")+
      #ylim(0.90, 1.65)+
      
      #scale_y_continuous(breaks=c(0,1,2.5,5,8))+
      #xlab("SNP")+
      # facet_grid(.~method)+
      ggtitle(paste0("Case-case OR of ",snp.temp)
      )+
      theme(plot.title = element_text(hjust=0.5,face="bold",size = 20),
            axis.text=element_text(face="bold",size=20),
            axis.title.x = element_text(face="bold",size=19),
            axis.title.y = element_text(face="bold"))+
      theme(axis.title.x=element_blank(),
            axis.text.y=element_blank(),
            axis.title.y=element_blank())
  )
  dev.off()
  
  png(paste0("./discovery_SNP/result/forest_plot/casecontrol_",snp.temp,".png"),height=15,width = 17,res=300,units="cm")
  data.plot.snp2 = data.plot%>% filter(SNP==snp.temp) %>% filter(Var%in%casecon.car)
  temp.case.char <- factor(data.plot.snp2$Var)
  data.plot.snp2$Var <- factor(temp.case.char,levels(temp.case.char)[c(1,6,2,4,5,3)])
  print(
    ggplot(data.plot.snp2,aes(x=Var,y=OR,ymin=low,ymax=high,colour=Var))+
      geom_pointrange()+
      #scale_colour_manual(values=c("#386cb0","#fdb462"))+
      #geom_line(yintercept=1,lty=2)+
      coord_flip()+
      theme_bw()+
      geom_hline(yintercept=1,size=1,lty=2)+
      theme(legend.position="none")+
      #ylab("Case-control odds ratio")+
      #xlab("Intrinsic subtypes")+
      #ylim(0.65, 1.65)+
      #scale_y_continuous(breaks=c(0,1,2.5,5,8))+
      #xlab("SNP")+
      # facet_grid(.~method)+
      ggtitle(paste0("Case-control OR of ",snp.temp)
      )+
      theme(plot.title = element_text(hjust=0.5,face="bold",size = 20),
            axis.text=element_text(face="bold",size=20),
            axis.title.x = element_text(face="bold",size=10),
            axis.title.y = element_text(face="bold",size=10))+
      theme(axis.title.x=element_blank(),
            axis.text.y=element_blank(),
            axis.title.y=element_blank())
  )
  dev.off()
  
}
#adjust for specific snp rs7924772
i <- 26
snp.temp <- data$rs_id[i]
png(paste0("./discovery_SNP/result/forest_plot/casecontrol_",snp.temp,".png"),height=15,width = 17,res=300,units="cm")
data.plot.snp2 = data.plot%>% filter(SNP==snp.temp) %>% filter(Var%in%casecon.car)
temp.case.char <- factor(data.plot.snp2$Var)
data.plot.snp2$Var <- factor(temp.case.char,levels(temp.case.char)[c(1,6,2,4,5,3)])
print(
  ggplot(data.plot.snp2,aes(x=Var,y=OR,ymin=low,ymax=high,colour=Var))+
    geom_pointrange()+
    #scale_colour_manual(values=c("#386cb0","#fdb462"))+
    #geom_line(yintercept=1,lty=2)+
    coord_flip()+
    theme_bw()+
    geom_hline(yintercept=1,size=1,lty=2)+
    theme(legend.position="none")+
    #ylab("Case-control odds ratio")+
    #xlab("Intrinsic subtypes")+
    ylim(0.89, 1.07)+
    #scale_y_continuous(breaks=c(0,1,2.5,5,8))+
    #xlab("SNP")+
    # facet_grid(.~method)+
    ggtitle(paste0("Case-control OR of ",snp.temp)
    )+
    theme(plot.title = element_text(hjust=0.5,face="bold",size = 20),
          axis.text=element_text(face="bold",size=20),
          axis.title.x = element_text(face="bold",size=10),
          axis.title.y = element_text(face="bold",size=10))+
    theme(axis.title.x=element_blank(),
          axis.text.y=element_blank(),
          axis.title.y=element_blank())
)
dev.off()


i <- 28
snp.temp <- data$rs_id[i]
png(paste0("./discovery_SNP/result/forest_plot/casecontrol_",snp.temp,".png"),height=15,width = 17,res=300,units="cm")
data.plot.snp2 = data.plot%>% filter(SNP==snp.temp) %>% filter(Var%in%casecon.car)
temp.case.char <- factor(data.plot.snp2$Var)
data.plot.snp2$Var <- factor(temp.case.char,levels(temp.case.char)[c(1,6,2,4,5,3)])
print(
  ggplot(data.plot.snp2,aes(x=Var,y=OR,ymin=low,ymax=high,colour=Var))+
    geom_pointrange()+
    #scale_colour_manual(values=c("#386cb0","#fdb462"))+
    #geom_line(yintercept=1,lty=2)+
    coord_flip()+
    theme_bw()+
    geom_hline(yintercept=1,size=1,lty=2)+
    theme(legend.position="none")+
    #ylab("Case-control odds ratio")+
    #xlab("Intrinsic subtypes")+
    ylim(0.90, 1.07)+
    #scale_y_continuous(breaks=c(0,1,2.5,5,8))+
    #xlab("SNP")+
    # facet_grid(.~method)+
    ggtitle(paste0("Case-control OR of ",snp.temp)
    )+
    theme(plot.title = element_text(hjust=0.5,face="bold",size = 20),
          axis.text=element_text(face="bold",size=20),
          axis.title.x = element_text(face="bold",size=10),
          axis.title.y = element_text(face="bold",size=10))+
    theme(axis.title.x=element_blank(),
          axis.text.y=element_blank(),
          axis.title.y=element_blank())
)
dev.off()
i <- 29
snp.temp <- data$rs_id[i]
png(paste0("./discovery_SNP/result/forest_plot/casecase_",snp.temp,".png"),height=15,width = 17,res=300,units="cm")
data.plot.snp1 = data.plot%>% filter(SNP==snp.temp) %>% filter(Var%in%case.char)
temp.case.char <- factor(data.plot.snp1$Var)
data.plot.snp1$Var <- factor(temp.case.char,levels(temp.case.char)[c(2,3,4,1)])
print(
  ggplot(data.plot.snp1,aes(x=Var,y=OR,ymin=low,ymax=high,colour=Var))+
    geom_pointrange()+
    #scale_colour_manual(values=c("#386cb0","#fdb462"))+
    #geom_line(yintercept=1,lty=2)+
    coord_flip()+
    theme_bw()+
    geom_hline(yintercept=1,size=1,lty=2)+
    theme(legend.position="none")+
    #ylab("Case-case odds ratio")+
    #xlab("Tumor characteristics")+
    ylim(0.89, 1.02)+
    
    #scale_y_continuous(breaks=c(0,1,2.5,5,8))+
    #xlab("SNP")+
    # facet_grid(.~method)+
    ggtitle(paste0("Case-case OR of ",snp.temp)
    )+
    theme(plot.title = element_text(hjust=0.5,face="bold",size = 20),
          axis.text=element_text(face="bold",size=20),
          axis.title.x = element_text(face="bold",size=19),
          axis.title.y = element_text(face="bold"))+
    theme(axis.title.x=element_blank(),
          axis.text.y=element_blank(),
          axis.title.y=element_blank())
)
dev.off()