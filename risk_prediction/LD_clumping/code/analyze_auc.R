#-------------------------------------------------------------------
# Update Date: 11/24/2018
# Create Date: 11/24/2018
# Goal: analyze auc
# Author: Haoyu Zhang
#-------------------------------------------------------------------
#---------------------------------------#---------------------------------------
load("./risk_prediction/LD_clumping/result/auc.result.rdata")
head(auc.result)
auc.result %>% filter(subtypes=="Luminal_B")
names.subtypes <- c("Luminal_A",
                    "Luminal_B",
                    "Luminal_B_HER2Neg",
                    "HER2Enriched",
                    "TripleNeg")
method.names <- c("standard","two-stage","eb")

total <- length(names.subtypes)*length(method.names)
method.out <- rep("c",total)
subtypes.out <- rep("c",total)
auc.out <- rep(0,total)
ind <- 1
for(i in 1:length(names.subtypes)){
  for(j in 1:length(method.names)){
    auc.out[ind] <- auc.result %>% filter(subtypes==names.subtypes[i]&
                            method==method.names[j]) %>%
      select(auc) %>% 
      max
    method.out[ind] <- method.names[j]
    subtypes.out[ind] <- names.subtypes[i]
    ind <- ind + 1
  }
}
result.temp <- data.frame(method.out,subtypes.out,auc.out)
colnames(result.temp) <- c("method","subtypes","auc")
auc.table <- NULL

auc.table <- NULL
  

result.temp %>% filter(
                         subtypes==names.subtypes[j]) %>% 
select(auc)
