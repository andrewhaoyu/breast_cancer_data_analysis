setwd("/spin1/users/zhangh24/breast_cancer_data_analysis/")
load(paste0("./genetic_correlation/ONCO/result/result.clean.completeglm.Rdata"))
load(paste0("./genetic_correlation/result/ICOG.result.clean.completeglm.Rdata"))
colnames(ICOG.result.clean)[11:40] <- paste0("ICOG",c(1:30))
colnames(ONCO.result.clean)[11:40] <- paste0("ONCO",c(1:30))
icog.onco.merge <- merge(ICOG.result.clean,
                         ONCO.result.clean,
                         by.x="chr.pos",
                         by.y="chr.pos",
                         all.x=F,
                         all.y=F)
save(icog.onco.merge,file=paste0("./genetic_correlation/ICOG/result/icog.onco.merge.completeglm.Rdata"))
