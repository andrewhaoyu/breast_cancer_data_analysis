load("/spin1/users/zhangh24/breast_cancer_data_analysis/discovery_SNP/conditional_analysis/result/conditional.results.2nd.Rdata")

load("/spin1/users/zhangh24/breast_cancer_data_analysis/discovery_SNP/conditional_analysis/result/all.conditional.snps.Rdata")

n.2nd <- nrow(conditional.results.2nd)

idx.2nd <- rep(0,n.2nd)

for(i in 1:n.2nd){

    idx <- which(all.conditional.snps$SNP.ONCO==conditional.results.2nd$SNP.ONCO[i])
    idx.2nd[i] <- idx    

  
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


icog.id.2nd <- conditional.snp.list.icog.clean[[1]][idx.2nd]
onco.id.2nd <-conditional.snp.list.onco.clean[[1]][idx.2nd]

test <- cbind(icog.id.2nd,onco.id.2nd)
all.equal(icog.id.2nd,conditional.results.2nd[,1])
all.equal(onco.id.2nd,conditional.results.2nd[,2])



snpvalue.icog.2nd <- conditional.snp.list.icog.clean[[2]][,idx.2nd]
snpvalue.onco.2nd <- conditional.snp.list.onco.clean[[2]][,idx.2nd]

icog.2nd <- list(icog.id.2nd,
                   snpvalue.icog.2nd)
onco.2nd <- list(onco.id.2nd,
                   snpvalue.onco.2nd)
save(icog.2nd,file="/spin1/users/zhangh24/breast_cancer_data_analysis/discovery_SNP/conditional_analysis/result/icog.2nd.Rdata")
save(onco.2nd,file="/spin1/users/zhangh24/breast_cancer_data_analysis/discovery_SNP/conditional_analysis/result/onco.2nd.Rdata")
save(conditional.results.2nd,file="/spin1/users/zhangh24/breast_cancer_data_analysis/discovery_SNP/conditional_analysis/result/conditional.results.2nd")

