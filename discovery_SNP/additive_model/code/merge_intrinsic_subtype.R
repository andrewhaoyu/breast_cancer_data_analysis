result <- NULL
for(i1 in 1:28){
  load(paste0("./discovery_SNP/additive_model/result/intrinsic_subtype_",i1,".Rdata"))
  result <- rbind(result,test.result.second.wald)
}
SNP <- c(colnames(icog.julie),colnames(discovery.snp.icog)[1:18])
final.result <- cbind(SNP,result)
write.csv(final.result,file= "./discovery_SNP/additive_model/result/intrinsic_subtype_final.csv")
