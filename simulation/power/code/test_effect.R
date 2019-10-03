if(s==1){
  #theta_test <- c(0.05,0,0,0,0)
  theta_test <- c(0.08,0,0,0,0)
  #theta_test <- c(0.25,0,0,0,0)
}else if(s==2){
  #theta_test <- c(0,0.05,0,0,0)
  theta_test <- c(0,0.08,0,0,0)
}else{
  #theta_test <- c(c(0,0.05),rnorm(3,0,0.02))
  theta_test <- c(c(0,0.08),rnorm(3,0,0.02))
}

    temp.simu <- SimulateDataPower(theta_intercept,theta_test,theta_covar,n)
    y.pheno.mis <- temp.simu[[1]]
    G <- temp.simu[[2]]
    x_covar <- temp.simu[[3]]
model <- glm(y.pheno.mis[,1]~G+x_covar,family=binomial())    
summary(model)








ggplot(data.m,aes(x=Sample_size,y=value,fill=variable))+
  geom_bar(stat="identity",position=position_dodge(),
           alpha=0.8)+
  scale_fill_Publication()+
  theme_Publication()+
  scale_x_discrete(name = "Sample size", limits=c("5000","25000","50000",
                                                  "1e+05"))+
  scale_y_continuous(name="Power",limits = c(0,1))+
  ggtitle(Heter_pattern_temp)+
  theme(legend.position="none")







ggplot(data.m,aes(x=Sample_size,y=value,fill=variable))+
  geom_bar(stat="identity",position=position_dodge(),alpha=0.8)+
  scale_fill_Publication()+
  theme_Publication()+
  scale_x_discrete(name="Sample_size",limits=c("5000","25000","50000","1e+05"))+
  scale_y_continuous(name="Power",limits=c(0,1))+
  ggtitle("No Heterogeneity")+
  theme(legend.position="none")






ggplot(data.m,aes(x=Sample_size,y=value,fill=variable))+
  geom_bar(stat="identity",position=position_dodge())+
  scale_fill_Publication()+
  theme_Publication()+
  scale_x_discrete(name="Sample Size",limits=c("5000","25000",
                                                "50000","1e+05"))+
  scale_y_continuous(name="Power",limits=c(0,1))+
  ggtitle("No heterogeneity")+
  theme(legend.position="right")


ggplot(data.m,aes(x=Sample_size,y=value,fill=variable))+
  geom_bar(stat="identity",position=position_dodge())+
  theme_Publication()+
  scale_fill_brewer(palette="Greys")+
  scale_x_discrete(name="Sample Size",limits=c("5000","25000",
                                               "50000","1e+05"))+
  scale_y_continuous(name="Power",limits=c(0,1))+
  ggtitle("No heterogeneity")+
  theme(legend.position="right")
  