#prepare the data for creating ideogram using Phenogram
setwd("/Users/zhangh24/GoogleDrive/breast_cancer_data_analysis")
data <- read.csv("./data/212_known_discovery_snp_paper_order.csv")
cate <- rep(0,nrow(data))
cate[1:178] <- c("Known variants")
cate[179:200] <- c("a")
cate[201:208] <- c("b")
cate[209:210] <- c("c")
data <- data.frame(data,cate,stringsAsFactors = F)
colnames(data) <- c("snp","chr","pos","phenotype")
write.table(data,file = paste0("./discovery_SNP/ideogram_plot/result/ideogram_data.txt"),col.names = T,quote=F,row.names = F,sep = "\t")

