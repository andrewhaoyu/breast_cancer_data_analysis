setwd('/spin1/users/zhangh24/breast_cancer_data_analysis/')
library(data.table)
data1 <- fread("./data/iCOGS_euro_v10_10232017.csv",header=T)
data2 <- fread("./data/Onco_euro_v10_10232017.csv",header=T)
all.countries <- unique(c(data1$StudyCountry,data2$StudyCountry))
n.coun <- length(all.countries)
n.snp <- 19
places <- 5
sensitivity.result <- list()
for(i1 in 1:n.snp){
  print(i1)
ER.low <- ER.high <- ER.odds <- rep(0,n.coun)

PR.low <- PR.high <- PR.odds <- rep(0,n.coun)

HER2.low <- HER2.high <- HER2.odds <- rep(0,n.coun)

grade.low <- grade.high <- grade.odds <- rep(0,n.coun)


    for(i2 in 1:20){
    
    load(paste0("./discovery_SNP/sensitivity_analysis/result/sensitivity_analysis",i1,"_",i2,".Rdata"))
      logodds <- result[[1]]
      sigma <- result[[2]]
      var.logodds <- diag(sigma)
      logodds.low <- logodds-1.96*sqrt(var.logodds)
      logodds.high <- logodds+1.96*sqrt(var.logodds)
      odds <- exp(logodds)
      odds.low <- exp(logodds.low)
      odds.high <- exp(logodds.high)
      odds <- round(odds,places)
      odds.low <- round(odds.low,places)
      odds.high <- round(odds.high,places)
      
      ER.low[i2] <- odds.low[2]
      PR.low[i2] <- odds.low[3]
      HER2.low[i2] <- odds.low[4]
      grade.low[i2] <- odds.low[5]
      ER.odds[i2] <-odds[2]
      PR.odds[i2] <- odds[3]
      HER2.odds[i2] <- odds[4]
      grade.odds[i2] <- odds[5]
      ER.high[i2] <- odds.high[2]
      PR.high[i2] <- odds.high[3]
      HER2.high[i2] <- odds.high[4]
      grade.high[i2] <- odds.high[5]
      
    
      result <- data.frame(ER.odds,ER.low,ER.high,
                           PR.odds,PR.low,PR.high,
                           HER2.odds,HER2.low,HER2.high,
                           grade.odds,grade.low,grade.high) 
    
    
    }

sensitivity.result[[i1]] <- result


}


save(sensitivity.result,file=paste0("./discovery_SNP/sensitivity_analysis/result/sensitivity_result.Rdata"))
