#goal: signed-rank test for gahna

data <- read.csv("./data/GWAS_Gahna.csv",stringsAsFactors = F)
idx <- which(data[,ncol(data)]==1)
data.n <- data[idx,]
library(dplyr)
library(MASS)
wilcox.test(log(as.numeric(data.n$BCAC_RA_OverallOR)),log(as.numeric(data.n$Ghana_RA_Overall_OR)),paired=T)


  wilcox.test(log(as.numeric(data.n$BCAC_RA_Erpos_OR)),log(as.numeric(data.n$Ghana_RA_Erpos_OR)),paired=T)

  
  wilcox.test(log(as.numeric(data.n$BCAC_RA_Eneg_OR)),log(as.numeric(data.n$Ghana_RA_Eneg_OR)),paired=T)
  