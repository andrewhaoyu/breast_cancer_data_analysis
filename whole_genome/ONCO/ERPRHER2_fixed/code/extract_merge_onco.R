load("/spin1/users/zhangh24/breast_cancer_data_analysis/whole_genome/ICOG/ERPRHER2_fixed/result/extract_list.Rdata")
setwd("/spin1/users/zhangh24/breast_cancer_data_analysis/")
pheno.file <- "./data/pheno.onco"
load(pheno.file)
n.sub = nrow(pheno)
y.pheno.mis2 <- cbind(pheno$Behaviour1,pheno$PR_status1,pheno$ER_status1,pheno$HER2_status1)
colnames(y.pheno.mis2) = c("Behaviour1","PR_status1",
                           "ER_status1","HER2_status1")
n <- nrow(y.pheno.mis2)
n.snps <- nrow(extract.list)
extract.result <- matrix(0,n,n.snps)
snp.id <- rep(0,n.snps)
total <- 0
for(i in 491:567){
  print(i)
  load(paste0("/spin1/users/zhangh24/breast_cancer_data_analysis/whole_genome/ONCO/ERPRHER2_fixed/result/ERPRHER2_fixed_onco_extract",i,".Rdata"))
  if(is.null(result[[1]])==0){
    temp <- length(result[[1]])
    snp.id[(1:temp)+total] <- result[[1]]
    extract.result[,(1:temp)+total] <- result[[2]]
    total <- temp+total
  }
  
}