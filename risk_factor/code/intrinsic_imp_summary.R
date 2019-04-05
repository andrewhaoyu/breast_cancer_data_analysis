#goal summary the results from intrinsic multiple imputation
setwd('/Users/zhangh24/GoogleDrive/breast_cancer_data_analysis')
load(paste0("./risk_factor/result/intrinsic_imp/intrinsic_imp_merge.Rdata"))

coef.1 <- result[[1]]
covar.1 <- result[[2]]
#interested in first 17 variable
n.in <- 16
colnames(coef.1)[1:n.in] <- c(paste0("agemenarche_cat",c(1,2,3)),
                              paste0("parity_cat",c(1,2,3,4)),
                              paste0("mensagelast_cat",c(1,2)),
                              paste0("agefftp_cat",c(2,3,4)),
                              paste0("breastmos_cat",c(2,3,4,5)))
# paste0("lastchildage_cat",c(2,3,4,9)))



#reorganize the covariance matrix to match coef data.frame
#columns are the variables
n.var <- ncol(coef.1)
#rows are the 5 intrinsic subtypes plus 9 
n.sub <- 5
var.1 <- matrix(diag(covar.1),nrow=n.sub)

result <- GenerateP(coef.1,var.1)
p <- result[[1]]
low.95 <- result[[2]]
high.95 <- result[[3]]

coef.1.n <- coef.1[,1:n.in]
var.1.n <- var.1[,1:n.in]
p.n <- p[,1:n.in]
low.95.n <- low.95[,1:n.in]
high.95.n <- high.95[,1:n.in]

coef.l <- as.vector(t(coef.1.n))
var.l <- as.vector(t(var.1.n))
p.l <- as.vector(t(p.n))
low.95.l <- as.vector(t(low.95.n))
high.95.l <- as.vector(t(high.95.n))
variable.l <- rep(colnames(coef.1.n),n.sub)
subtypes.l <- c(rep("Luminal A-like",ncol(coef.1.n)),
                rep("Luminal B,HER2-negative-like",ncol(coef.1.n)),
                rep("Luminal B-like",ncol(coef.1.n)),
                rep("HER2 enriched-like",ncol(coef.1.n)),
                rep("TN",ncol(coef.1.n)))
                #rep("missing",ncol(coef.1.n)))
#create the variable category for colour
variable.cat <- substr(variable.l,1,nchar(variable.l)-1)

new.data <- data.frame(variable.l,subtypes.l,
                       coef.l,
                       var.l,
                       p.l,
                       low.95.l,
                       high.95.l,
                       exp(coef.l),
                       exp(low.95.l),
                       exp(high.95.l),
                       variable.cat,
                       stringsAsFactors = F)
colnames(new.data) <- c("variable",
                        "subtypes",
                        "coef",
                        "variance",
                        "p-value",
                        "low95.coef",
                        "high95.coef",
                        "OR",
                        "ORlow95",
                        "ORhigh95",
                        "variable_cat")

subtypes <- c("Luminal A-like","Luminal B,HER2-negative-like",
              "Luminal B-like",
              "HER2 enriched-like",
              "TN")
new.data.c <- new.data
write.csv(new.data.c,file = "./risk_factor/result/intrinsic_imp.csv")
for(i in 1:length(subtypes)){
  subtype.temp = subtypes[i]
  png(paste0("./risk_factor/result/intrinsic_imp/intrinic_imp_result_",subtype.temp,".png"),height=20,width = 15,res=300,units="cm")
  new.data.plot <- new.data.c %>% filter(subtypes==subtype.temp)
  print(
    ggplot(new.data.plot,aes(x=variable,y=OR,ymin=ORlow95,ymax=ORhigh95,colour=variable_cat))+
      geom_pointrange()+
      #scale_colour_manual(values=c("#386cb0","#fdb462"))+
      #geom_line(yintercept=1,lty=2)+
      coord_flip()+
      theme_bw()+
      geom_hline(yintercept=1,size=1,lty=2)+
      theme(legend.position="none")+
      ylab("Odds ratio")+
      #scale_y_continuous(breaks=c(0,1,2.5,5,8))+
      #xlab("SNP")+
      # facet_grid(.~method)+
      ggtitle(paste0("Forest plot for odds ratio of ",subtype.temp)
      )+
      theme(plot.title = element_text(hjust=0.5,face="bold"),
            axis.text=element_text(face="bold"),
            axis.title.x = element_text(face="bold"),
            axis.title.y = element_text(face="bold"))
  )
  dev.off()
  
}



GenerateP <- function(coef,sigma){
  z <- coef/sqrt(sigma)
  p <- 2*pnorm(-abs(z),lower.tail = T)
  low.95 <- coef-1.96*sqrt(sigma)
  high.95 <- coef+1.96*sqrt(sigma)
  return(list(p,low.95,high.95))
}
