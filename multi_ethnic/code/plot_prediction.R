#goal: plot prediction results
setwd('/Users/zhangh24/GoogleDrive/breast_cancer_data_analysis/')
library(ggplot2)
p.thr <- c(10^-8,5E-8,10^-7,5E-7,10^-6,5E-6,10^-5,5E-5,10^-4,5E-4,10^-3,5E-3,10^-2,0.1,0.3,0.5)
log10P <- -log10(p.thr)
log10P <- round(log10P,1)
n.cut <- length(p.thr)
prop.result <- matrix(0,3,n.cut)
n.snp.result <- matrix(0,3,n.cut)
vad.r2 <- rep(0,3)
pop <- c("European",
         "African",
         "Latino")
for(i in 1:3){
  load(paste0("./multi_ethnic/result/LDP_summary_",i)) 
  n.snp.result[i,] <- round(LDP.result[[1]],0)
  prop.result[i,] <- round(LDP.result[[2]],2)
  vad.r2[i] <- mean(LDP.result[[5]])
  r2.test <- LDP.result[[3]]
  data <- data.frame(p.thr,r2.test,prop.result[i,])
  colnames(data) <- c("Pthr","R2","Prop")
  
  
  png(paste0("./multi_ethnic/result/prediction_",pop[i],".png"),width = 24, height=12,unit = "cm",res = 300)
 print({
   ggplot(data)+geom_line(aes(-log10(Pthr),R2))+
     fte_theme()+
     xlab("-log10(P-value)")+
     ylab("Adjusted R2")+
     ggtitle(paste0("Prediction result for ",pop[i]," population"))+
     scale_x_discrete(limit= log10P)
 }) 
 dev.off()
 png(paste0("./multi_ethnic/result/proportion_",pop[i],".png"),width = 24, height=12,unit = "cm",res = 300)
 print({
   ggplot(data)+geom_line(aes(-log10(Pthr),Prop))+
     fte_theme()+
     xlab("-log10(P-value)")+
     ylab("Proportion of causal SNPs in the selected SNPs")+
     ggtitle(paste0("Prediction result for ",pop[i]," population"))+
     scale_x_discrete(limit= log10P)
 }) 
 dev.off() 
  
}
write.csv(n.snp.result,"./multi_ethnic/result/number_of_selected_SNP.csv")
write.csv(prop.result,"./multi_ethnic/result/causal_proportion.csv")

