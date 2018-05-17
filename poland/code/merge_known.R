heter.result.com <- NULL
for(i1 in 1:178){
  load(paste0("./poland/result/known_snps/heter_result_",i1,".Rdata"))  
  heter.result.com <- rbind(heter.result.com,heter.result)
}

write.csv(heter.result.com,file=paste0("./poland/result/known_snps/known_snps_result.csv"),row.names=F,quote=F)
