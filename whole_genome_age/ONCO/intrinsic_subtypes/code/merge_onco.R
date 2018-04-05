
setwd("/spin1/users/zhangh24/breast_cancer_data_analysis/")
filedir <- './whole_genome_age/ONCO/intrinsic_subtypes/result/'
files <- dir(filedir,pattern="intrinsic_subytpe_onco")
total <- 567*5
missingid <- matrix(0,total,2)
temp <- 0
for(i1 in 1:567){
  print(i1)
  for(i2 in 1:5){
    text <- paste0("intrinsic_subytpe_onco",i1,"_",i2)
    if((text%in%files)==F){
      temp <- temp+1
      missingid[temp,] <- c(i1,i2)
    }
  }
}
missingid <- missingid[1:temp,]
length(unique(missingid[,1]))
onco.unique.resubmit <- unique(missingid[,1])
save(onco.unique.resubmit,file="./whole_genome_age/ONCO/intrinsic_subtypes/result/onco.unique.resubmit.Rdata")
submit <- rep("c",length(onco.unique.resubmit)*15)
temp <- 1
for(i in 1:length(onco.unique.resubmit)){
  for(j in 1:15){
    submit[temp] <- paste0("Rscript /spin1/users/zhangh24/breast_cancer_data_analysis/whole_genome_age/ONCO/intrinsic_subtypes/code/intrinsic_subtype_onco.R ",onco.unique.resubmit[i]," ",j)
    temp <- temp+1
  }
  
}
write.table(submit,file="/spin1/users/zhangh24/breast_cancer_data_analysis/whole_genome_age/ONCO/intrinsic_subtypes/code/onco_resubmit.sh",
            row.names=F,quote=F,col.names=F)
