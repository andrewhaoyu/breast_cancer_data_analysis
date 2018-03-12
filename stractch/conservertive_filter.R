library(data.table)
discovery_snp <- as.data.frame(fread("./data/discovery_snps_annotated_clean.csv",header=T))
fine_mapping <- as.data.frame(fread("./data/fine_mapping_annotated_clean.csv"))

dis <- 10^6
check.data <- NULL
for(i in 11:29){
  pos <- discovery_snp$position[i]
  CHR <- discovery_snp$CHR[i]
  for(j in 1:178){
    pos.known <- fine_mapping$position[j]
    CHR.known <- fine_mapping$CHR[j]
    if(CHR==CHR.known&pos>=(pos.known-dis)&pos<=(pos.known+dis)){
      print(c(i,j))
      temp1 <- discovery_snp[i,c(1,3,2)]
      colnames(temp1) <- c("SNP","CHR","Position")
      temp2 <- fine_mapping[j,c(1,3,4)]
      colnames(temp2) <- c("SNP","CHR","Position")
      result <- rbind(temp1,temp2)
      check.data <- rbind(check.data,result)      
    }
      
  }
}
write.csv(check.data,file="./data/check_SNPs.csv",row.names = F,col.names = T)
