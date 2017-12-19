load("/spin1/users/zhangh24/breast_cancer_data_analysis/discovery_SNP/conditional_analysis/result/conditional.results.5th.Rdata")

load("/spin1/users/zhangh24/breast_cancer_data_analysis/discovery_SNP/conditional_analysis/result/all.conditional.snps.Rdata")

n.5th <- nrow(conditional.results.5th)

idx.5th <- rep(0,n.5th)

for(i in 1:n.5th){
  
  idx <- which(all.conditional.snps$SNP.ONCO==conditional.results.5th$SNP.ONCO[i])
  idx.5th[i] <- idx    
  
  
}


load("/spin1/users/zhangh24/breast_cancer_data_analysis/discovery_SNP/conditional_analysis/result/conditional.snp.list.icog.Rdata")
load("/spin1/users/zhangh24/breast_cancer_data_analysis/discovery_SNP/conditional_analysis/result/conditional.snp.list.onco.Rdata")
load("/spin1/users/zhangh24/breast_cancer_data_analysis/discovery_SNP/conditional_analysis/result/all.conditional.snps.Rdata")

n.icog <- length(conditional.snp.list.icog[[1]])
n.onco <- length(conditional.snp.list.onco[[1]])

idx.match.icog <- match(all.conditional.snps$SNP.ICOGS,conditional.snp.list.icog[[1]])

test.icog <- conditional.snp.list.icog[[1]][idx.match.icog]
test.icog.value <- conditional.snp.list.icog[[2]][,idx.match.icog]

conditional.snp.list.icog.clean <- list(test.icog,
                                        test.icog.value)
idx.match.onco <- match(all.conditional.snps$SNP.ONCO,conditional.snp.list.onco[[1]])

test.onco <- conditional.snp.list.onco[[1]][idx.match.onco]
test.onco.value <- conditional.snp.list.onco[[2]][,idx.match.onco]
conditional.snp.list.onco.clean <- list(test.onco,
                                        test.onco.value)


icog.id.5th <- conditional.snp.list.icog.clean[[1]][idx.5th]
onco.id.5th <-conditional.snp.list.onco.clean[[1]][idx.5th]

test <- cbind(icog.id.5th,onco.id.5th)
all.equal(icog.id.5th,conditional.results.5th[,1])
all.equal(onco.id.5th,conditional.results.5th[,2])



snpvalue.icog.5th <- conditional.snp.list.icog.clean[[2]][,idx.5th]
snpvalue.onco.5th <- conditional.snp.list.onco.clean[[2]][,idx.5th]

icog.5th <- list(icog.id.5th,
                 snpvalue.icog.5th)
onco.5th <- list(onco.id.5th,
                 snpvalue.onco.5th)
save(icog.5th,file="/spin1/users/zhangh24/breast_cancer_data_analysis/discovery_SNP/conditional_analysis/result/icog.5th.Rdata")
save(onco.5th,file="/spin1/users/zhangh24/breast_cancer_data_analysis/discovery_SNP/conditional_analysis/result/onco.5th.Rdata")
save(conditional.results.5th,file="/spin1/users/zhangh24/breast_cancer_data_analysis/discovery_SNP/conditional_analysis/result/conditional.results.5th")

