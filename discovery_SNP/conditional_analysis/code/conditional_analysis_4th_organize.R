load("/spin1/users/zhangh24/breast_cancer_data_analysis/discovery_SNP/conditional_analysis/result/conditional.results.4th.Rdata")

load("/spin1/users/zhangh24/breast_cancer_data_analysis/discovery_SNP/conditional_analysis/result/all.conditional.snps.Rdata")

n.4th <- nrow(conditional.results.4th)

idx.4th <- rep(0,n.4th)

for(i in 1:n.4th){
  
  idx <- which(all.conditional.snps$SNP.ONCO==conditional.results.4th$SNP.ONCO[i])
  idx.4th[i] <- idx    
  
  
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


icog.id.4th <- conditional.snp.list.icog.clean[[1]][idx.4th]
onco.id.4th <-conditional.snp.list.onco.clean[[1]][idx.4th]

test <- cbind(icog.id.4th,onco.id.4th)
all.equal(icog.id.4th,conditional.results.4th[,1])
all.equal(onco.id.4th,conditional.results.4th[,2])



snpvalue.icog.4th <- conditional.snp.list.icog.clean[[2]][,idx.4th]
snpvalue.onco.4th <- conditional.snp.list.onco.clean[[2]][,idx.4th]

icog.4th <- list(icog.id.4th,
                 snpvalue.icog.4th)
onco.4th <- list(onco.id.4th,
                 snpvalue.onco.4th)
save(icog.4th,file="/spin1/users/zhangh24/breast_cancer_data_analysis/discovery_SNP/conditional_analysis/result/icog.4th.Rdata")
save(onco.4th,file="/spin1/users/zhangh24/breast_cancer_data_analysis/discovery_SNP/conditional_analysis/result/onco.4th.Rdata")
save(conditional.results.4th,file="/spin1/users/zhangh24/breast_cancer_data_analysis/discovery_SNP/conditional_analysis/result/conditional.results.4th")

