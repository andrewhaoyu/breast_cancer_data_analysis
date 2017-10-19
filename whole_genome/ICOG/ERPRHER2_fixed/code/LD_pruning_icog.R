#load("/spin1/users/zhangh24/breast_cancer_data_analysis/whole_genome/ICOG/ERPRHER2_fixed/result/extract_result.Rdata")

setwd("/spin1/users/zhangh24/breast_cancer_data_analysis/")
load("/spin1/users/zhangh24/breast_cancer_data_analysis/whole_genome/ONCO/ERPRHER2_fixed/result/extract_result.Rdata")

data2 <- read.csv("./data/ONCO_pruning.csv",header=T)
x.test.all.mis2 <- data2[,c(27:205)]

idx.control <- which(data2$Behaviour1==0)

x.test.all.mis2.control <- x.test.all.mis2[idx.control,]

load("/spin1/users/zhangh24/breast_cancer_data_analysis/whole_genome/ICOG/ERPRHER2_fixed/result/extract_list.Rdata")
# extract.list <- extract.list[-c(1700,1701,1702,1703,1705,1706,1707),]
# save(extract.list,file="/spin1/users/zhangh24/breast_cancer_data_analysis/whole_genome/ICOG/ERPRHER2_fixed/result/extract_list.Rdata")
#try <- duplicated(extract.list)

extract.result.onco <- extract.result[[2]]
idx.match <- match(extract.list$SNP.ONCO,extract.result[[1]])
extract.result.onco <- extract.result.onco[,idx.match]
extract.result.onco.control <- extract.result.onco[idx.control,]


idx.known.ld.flag <- NULL

for(i in 1:ncol(extract.result.onco.control)){
  print(i)
  temp <- extract.result.onco.control[,i]
  
  for(j in 1:179){
    check.ld <- cor(x.test.all.mis2.control[,j],temp)^2
    if(check.ld>=0.1){
      idx.known.ld.flag <- c(idx.known.ld.flag,i)
      break
    }
      
  }
  
  
  
}

extract.list.ld <- extract.list[idx.known.ld.flag,]

save(extract.list.ld,file="/spin1/users/zhangh24/breast_cancer_data_analysis/whole_genome/ICOG/ERPRHER2_fixed/result/extract.list.ld")
