load("/spin1/users/zhangh24/breast_cancer_data_analysis/discovery_SNP/conditional_analysis/result/conditional.results.3rd.Rdata")

load("/spin1/users/zhangh24/breast_cancer_data_analysis/discovery_SNP/conditional_analysis/result/all.conditional.snps.Rdata")

n.3rd <- nrow(conditional.results.3rd)

idx.3rd <- rep(0,n.3rd)

for(i in 1:n.3rd){
  
  idx <- which(all.conditional.snps$SNP.ONCO==conditional.results.3rd$SNP.ONCO[i])
  idx.3rd[i] <- idx    
  
  
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


icog.id.3rd <- conditional.snp.list.icog.clean[[1]][idx.3rd]
onco.id.3rd <-conditional.snp.list.onco.clean[[1]][idx.3rd]

test <- cbind(icog.id.3rd,onco.id.3rd)
all.equal(icog.id.3rd,conditional.results.3rd[,1])
all.equal(onco.id.3rd,conditional.results.3rd[,2])



snpvalue.icog.3rd <- conditional.snp.list.icog.clean[[2]][,idx.3rd]
snpvalue.onco.3rd <- conditional.snp.list.onco.clean[[2]][,idx.3rd]

icog.3rd <- list(icog.id.3rd,
                 snpvalue.icog.3rd)
onco.3rd <- list(onco.id.3rd,
                 snpvalue.onco.3rd)
save(icog.3rd,file="/spin1/users/zhangh24/breast_cancer_data_analysis/discovery_SNP/conditional_analysis/result/icog.3rd.Rdata")
save(onco.3rd,file="/spin1/users/zhangh24/breast_cancer_data_analysis/discovery_SNP/conditional_analysis/result/onco.3rd.Rdata")
save(conditional.results.3rd,file="/spin1/users/zhangh24/breast_cancer_data_analysis/discovery_SNP/conditional_analysis/result/conditional.results.3rd")

