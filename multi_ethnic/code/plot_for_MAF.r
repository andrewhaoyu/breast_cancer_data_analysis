#plot the MAF freq between EUR and AFR in KG pruned
load("/Users/zhangh24/GoogleDrive/breast_cancer_data_analysis/multi_ethnic/result/pruned_MAF.Rdata")
plot(pruned.snp.clean$MAF.EUR,pruned.snp.clean$all.snp.AFR.MAF.AFR)
hist(pruned.snp.clean$MAF.EUR)
hist(pruned.snp.clean$all.snp.AFR.MAF.AFR)

colnames(pruned.snp.clean)[7] <- "MAF.AFR"
library(ggplot2)
library()
ggplot(pruned.snp.clean,aes(x=MAF.EUR))+
  geom_bar(stat="bin",fill="dodgerblue4")+
  theme_minimal()+
  xlab("MAF of SNPs in European population")
ggplot(pruned.snp.clean,aes(x=MAF.AFR))+
  geom_bar(stat="bin",fill="#c0392b")+
  theme_minimal()+
  xlab("MAF of SNPs in African population")
