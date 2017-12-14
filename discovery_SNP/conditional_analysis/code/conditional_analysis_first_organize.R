load("/spin1/users/zhangh24/breast_cancer_data_analysis/discovery_SNP/conditional_analysis/result/conditional.results.first.Rdata")

load("/spin1/users/zhangh24/breast_cancer_data_analysis/discovery_SNP/conditional_analysis/result/all.conditional.snps.Rdata")

n.first <- nrow(conditional.results.first)

idx.first <- rep(0,n.first)

for(i in 1:n.first){
  if(i!=56){
    print(i)
    idx <- which(all.conditional.snps$SNP.ICOGS==conditional.results.first$SNP.ICOGS[i]&
                   all.conditional.snps$SNP.ONCO==conditional.results.first$SNP.ONCO[i])
    idx.first[i] <- idx    
  }else{
    idx <- which(all.conditional.snps$SNP.ONCO==conditional.results.first$SNP.ONCO[i])
    idx.first[i] <- idx    
  }
  
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


icog.id.first <- conditional.snp.list.icog.clean[[1]][idx.first]
onco.id.first <-conditional.snp.list.onco.clean[[1]][idx.first]

test <- cbind(icog.id.first,onco.id.first)
all.equal(icog.id.first,conditional.results.first[,1])
all.equal(onco.id.first,conditional.results.first[,2])



snpvalue.icog.first <- conditional.snp.list.icog.clean[[2]][,idx.first]
snpvalue.onco.first <- conditional.snp.list.onco.clean[[2]][,idx.first]

icog.first <- list(icog.id.first,
                   snpvalue.icog.first)
onco.first <- list(onco.id.first,
                   snpvalue.onco.first)
save(icog.first,file="/spin1/users/zhangh24/breast_cancer_data_analysis/discovery_SNP/conditional_analysis/result/icog.first.Rdata")
save(onco.first,file="/spin1/users/zhangh24/breast_cancer_data_analysis/discovery_SNP/conditional_analysis/result/onco.first.Rdata")
save(conditional.results.first,file="/spin1/users/zhangh24/breast_cancer_data_analysis/discovery_SNP/conditional_analysis/result/conditional.results.first")

