n.signals <- c(1,2,3,4)
n.locus <- c(36,13,3,4)
data <- data.frame(n.signals,n.locus)
colnames(data) <- c("signal","locus")
p <- ggplot(data,aes(signal,locus))+geom_bar(stat="identity",fill="#DD8888")+
  xlab("N signals")+ylab("N locus")+ggtitle("Independent signals per region")+theme_bw()

