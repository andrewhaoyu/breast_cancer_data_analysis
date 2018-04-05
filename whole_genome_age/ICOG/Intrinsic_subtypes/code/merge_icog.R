setwd("/spin1/users/zhangh24/breast_cancer_data_analysis/")
filedir <- './whole_genome_age/ICOG/Intrinsic_subtypes/result/'
files <- dir(filedir,pattern="intrinsic_subytpe_icog")
total <- 564*5
missingid <- matrix(0,total,2)
temp <- 0
for(i1 in 1:564){
  print(i1)
  for(i2 in 1:5){
    text <- paste0("intrinsic_subytpe_icog",i1,"_",i2)
    if((text%in%files)==F){
      temp <- temp+1
      missingid[temp,] <- c(i1,i2)
    }
  }
}
missingid <- missingid[1:temp,]
icog.unique.resubmit <- unique(missingid[,1])
