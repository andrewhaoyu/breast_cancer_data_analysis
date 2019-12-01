load("/spin1/users/zhangh24/breast_cancer_data_analysis/whole_genome_age/ICOG/Intrinsic_subtypes/result/meta_result_shared_1p_final.Rdata")

library(dat.talbe)
library(dplyr)


colnames(meta_result_shared_1p)[21:45] <- 
  paste0("cov",c(1:25))


#load CIMBA BCAC meta-analysis data
#chromosome 23 data are not meta-analyzed yet
CIMBA.result <- as.data.frame(fread("/spin1/users/zhangh24/breast_cancer_data_analysis/data/brca1_bcac_tn_meta.txt"))



CIMBA.BCAC <- merge(meta_result_shared_1p,CIMBA.result,by.x = "var_name",by.y = "MarkerName")
idx.am <- which(CIMBA.BCAC$a1!=toupper(CIMBA.BCAC$Allele2))
#reverse the allele that not aligned in bcac and CIMBA
temp1 = CIMBA.BCAC$a1[idx.am] 
temp2 = CIMBA.BCAC$a0[idx.am]
CIMBA.BCAC$a1[idx.am] = temp2
CIMBA.BCAC$a0[idx.am] = temp1

CIMBA.BCAC$Effect[idx.am] <- -CIMBA.BCAC$Effect[idx.am]
CIMBA.BCAC$Freq1[idx.am] <- 1-CIMBA.BCAC$Freq1[idx.am]
CIMBA.BCAC$MinFreq[idx.am] <- 1-CIMBA.BCAC$MinFreq[idx.am]
CIMBA.BCAC$MaxFreq[idx.am] <- 1-CIMBA.BCAC$MaxFreq[idx.am]

CIMBA.BCAC_update <- CIMBA.BCAC %>% 
  mutate(var_effect = StdErr^2) %>% 
  select(var_name,"5",cov25,Effect,var_effect)

#BCAC TN and CIMBA BRCA1 meta analysis
meta_effect <- rep(0,nrow(CIMBA.BCAC_update))
meta_std <- rep(0,nrow(CIMBA.BCAC_update))
meta_p <- rep(0,nrow(CIMBA.BCAC_update))
for(i in 1:nrow(CIMBA.BCAC_update)){
  if(i%%1000==0){
    print(i)  
  }
  
  result_temp <- MetaFixedPfunction_temp(CIMBA.BCAC_update[i,2:5],1)
  meta_effect[i] <- result_temp[[1]]
  meta_std[i] <- sqrt(result_temp[[2]])
  meta_p[i] <- as.numeric(result_temp[[3]][2])
  
}



CIMBA.BCAC_update$meta_effect <- meta_effect
CIMBA.BCAC_update$meta_std <- meta_std
CIMBA.BCAC_update$meta_p <- meta_p
#Some alleles change the order due to BCAC and CIMBA meta-analysis.
#we need to change them back
colnames(CIMBA.BCAC_update)[1] <- "MarkerName"
CIMBA.BCAC_update$a1_update <- CIMBA.BCAC$a1
CIMBA.BCAC_update$a0_update <- CIMBA.BCAC$a0
CIMBA.BCAC_update$Freq1_update <- CIMBA.BCAC$Freq1
CIMBA.BCAC_update$MinFreq_update <- CIMBA.BCAC$MinFreq
CIMBA.BCAC_update$MaxFreq_update <- CIMBA.BCAC$MaxFreq





#update CIMBA data
CIMBA.result.update <- left_join(CIMBA.result,CIMBA.BCAC_update,by ="MarkerName")
idx <- which(!is.na(CIMBA.result.update$Effect.y))
#use the update information allele1,allele2,freq1,minfreq,maxfreq,effect.x, stderr,p-value
CIMBA.result.update[idx,c(2,3,4,6,7,8,9,10)] <-
  CIMBA.result.update[idx,c(26,25,28,27,29,20,23,24)]
head(CIMBA.result.update)
#drop the uncessary data
CIMBA.result.update <- CIMBA.result.update[,c(1:17)]
write.table(CIMBA.result.update,file = "/spin1/users/zhangh24/breast_cancer_data_analysis/data/CIMBA_BCAC_meta_analysis_083019.txt",quote=F,row.names=F)



