heter.result.all <- NULL
for(i1 in 1:2){
  load(paste0("./poland/result/whole_genome/heter_result_",i1,".Rdata"))
  heter.result.all <- rbind(heter.result.all,heter.result)
}
