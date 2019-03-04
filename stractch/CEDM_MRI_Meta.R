#Goal CEDM MRI Sensitivity analysis
sensi <- read.csv("./data/CEDM_MRI_sensitivity.csv")
Meta95 <- function(p_vec,sample_vec){
  n <- length(p_vec)
  var_vec <- rep(0,n)
  for(i in 1:n){
    p <- p_vec[i]
    if(p==1){
      p= 0.99
    }
    var_vec[i] <- p*(1-p)/sample_vec[i]
  }
  meta_var <- sum(1/var_vec)^-1
  meta_p <- meta_var*sum(p_vec/var_vec)
  meta_p_high <- meta_p+1.96*sqrt(meta_var)
  meta_p_low <- meta_p-1.96*sqrt(meta_var)
  if(meta_p_high>=1){
    meta_p_high <- 1
  }
  return(c(meta_p,meta_p_low,meta_p_high))
}

CI95 <- function(p,sample){
  if(p==1){
    p=0.99
  }
  low <- p-1.96*sqrt(p*(1-p)/sample)
  high <- p+1.96*sqrt(p*(1-p)/sample)
  if(high>=1){
    high = 1
  }
  return(c(low,high))
}
n <- nrow(sensi)
MRI.95 <- matrix(0,nrow(sensi),2)
for(i in 1:n){
  MRI.95[i,] <- CI95(sensi[i,2],sensi[i,4])
  
}
colnames(MRI.95) <- c("low","high")
MRI.meta <- Meta95(sensi[,2],sensi[,4])
CESM.95 <- matrix(0,nrow(sensi),2)
for(i in 1:n){
  CESM.95[i,] <- CI95(sensi[i,3],sensi[i,4])
  
}
colnames(CESM.95) <- c("low","high")
CESM.meta <- Meta95(sensi[,3],sensi[,4])
label <- rep(c(as.character(sensi[,1]),"Summary"),2)
Sensitivity <- c(sensi[,2],MRI.meta[1],sensi[,3],CESM.meta[1])
low <- c(MRI.95[,1],MRI.meta[2],CESM.95[,1],CESM.meta[2])
high <- c(MRI.95[,2],MRI.meta[3],CESM.95[,2],CESM.meta[3])
Type <- rep(c(rep("ind",n),"sum"),2)
method <- factor(c(rep("MRI",n+1),rep("CESM",n+1)),
                 levels = c("MRI","CESM"))
#method <- as.factor(method)
levels(method) <- c("MRI","CESM")
label <- factor(label,
                levels = rev(c(as.character(sensi[,1]),"Summary")))


new.data <- data.frame(Sensitivity,low,high,Type,label,method)
library(ggplot2)
ggplot(new.data,aes(x=label,y=Sensitivity,ymin=low,ymax=high,shape=Type,colour=Type))+
  geom_pointrange()+
  scale_colour_manual(values=c("#386cb0","#fdb462"))+
  #geom_line(yintercept=1,lty=2)+
  coord_flip()+
  theme_bw()+
  geom_hline(yintercept=1,size=1,lty=2)+
  theme(legend.position="none")+
  ylab("Sensitivity")+
  #scale_y_continuous(limits=c(0.90,1.1))+
  xlab("Studies")+
  facet_grid(.~method)+
  ggtitle(paste0("Forest plot of sensitivity")
  )+
  theme(plot.title = element_text(hjust=0.5,face="bold"),
        axis.text=element_text(face="bold"),
        axis.title.x = element_text(face="bold"),
        axis.title.y = element_text(face="bold"))










sensi <- read.csv("./data/CEDM_MRI_specificity.csv")
Meta95 <- function(p_vec,sample_vec){
  n <- length(p_vec)
  var_vec <- rep(0,n)
  for(i in 1:n){
    p <- p_vec[i]
    if(p==1){
      p= 0.99
    }
    var_vec[i] <- p*(1-p)/sample_vec[i]
  }
  meta_var <- sum(1/var_vec)^-1
  meta_p <- meta_var*sum(p_vec/var_vec)
  meta_p_high <- meta_p+1.96*sqrt(meta_var)
  meta_p_low <- meta_p-1.96*sqrt(meta_var)
  if(meta_p_high>=1){
    meta_p_high <- 1
  }
  return(c(meta_p,meta_p_low,meta_p_high))
}

CI95 <- function(p,sample){
  if(p==1){
    p=0.99
  }
  low <- p-1.96*sqrt(p*(1-p)/sample)
  high <- p+1.96*sqrt(p*(1-p)/sample)
  if(high>=1){
    high = 1
  }
  return(c(low,high))
}
n <- nrow(sensi)
MRI.95 <- matrix(0,nrow(sensi),2)
for(i in 1:n){
  MRI.95[i,] <- CI95(sensi[i,2],sensi[i,4])
  
}
colnames(MRI.95) <- c("low","high")
MRI.meta <- Meta95(sensi[,2],sensi[,4])
CESM.95 <- matrix(0,nrow(sensi),2)
for(i in 1:n){
  CESM.95[i,] <- CI95(sensi[i,3],sensi[i,4])
  
}
colnames(CESM.95) <- c("low","high")
CESM.meta <- Meta95(sensi[,3],sensi[,4])
label <- rep(c(as.character(sensi[,1]),"Summary"),2)
Sensitivity <- c(sensi[,2],MRI.meta[1],sensi[,3],CESM.meta[1])
low <- c(MRI.95[,1],MRI.meta[2],CESM.95[,1],CESM.meta[2])
high <- c(MRI.95[,2],MRI.meta[3],CESM.95[,2],CESM.meta[3])
Type <- rep(c(rep("ind",n),"sum"),2)
method <- factor(c(rep("MRI",n+1),rep("CESM",n+1)),
                 levels = c("MRI","CESM"))
#method <- as.factor(method)
levels(method) <- c("MRI","CESM")
label <- factor(label,
                levels = rev(c(as.character(sensi[,1]),"Summary")))


new.data <- data.frame(Sensitivity,low,high,Type,label,method)
library(ggplot2)
ggplot(new.data,aes(x=label,y=Sensitivity,ymin=low,ymax=high,shape=Type,colour=Type))+
  geom_pointrange()+
  scale_colour_manual(values=c("#386cb0","#fdb462"))+
  #geom_line(yintercept=1,lty=2)+
  coord_flip()+
  theme_bw()+
  geom_hline(yintercept=1,size=1,lty=2)+
  theme(legend.position="none")+
  ylab("Specificity")+
  #scale_y_continuous(limits=c(0.90,1.1))+
  xlab("Studies")+
  facet_grid(.~method)+
  ggtitle(paste0("Forest plot of specificity")
  )+
  theme(plot.title = element_text(hjust=0.5,face="bold"),
        axis.text=element_text(face="bold"),
        axis.title.x = element_text(face="bold"),
        axis.title.y = element_text(face="bold"))
