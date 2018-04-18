args = commandArgs(trailingOnly = T)
i1 = as.numeric(args[[1]])
library(bc2)
library(bcutility)
load("/spin1/users/zhangh24/breast_cancer_data_analysis/discovery_SNP/conditional_analysis/result/conditional.check.data.Rdata")

y.pheno.mis1 = conditional.check.data$y.pheno.mis1
y.pheno.mis2 = conditional.check.data$y.pheno.mis2
x.covar.mis1 = conditional.check.data$x.covar.mis1
x.covar.mis2 = conditional.check.data$x.covar.mis2
x.test.all1 = conditional.check.data[[2*i1-1]]
x.test.all2 = conditional.check.data[[2*i1]]


z.standard = GenerateZstandard(y.pheno.mis1)
z.standard <- GenerateZstandard(y.pheno.mis1)
z.random.support <- cbind(1,z.standard[,1])
z.random.test <- z.standard[,2:4]

total.gene <- ncol(x.test.all1)-1
p.value <- rep(0,total.gene)

for(i in 1:total.gene){
  gene1 = x.test.all1[,i+1]
  x.covar1 = cbind(x.covar.mis1,x.test.all1[,1:i])
  gene2 = x.test.all2[,i+1]
  x.covar2 = cbind(x.covar.mis2,x.test.all2[,1:i])
  p.value[i] <- two_data_two_stage_random(y.pheno.mis1,
                                        gene1,
                                        x.covar1,
                                        y.pheno.mis2,
                                        gene2,
                                        x.covar2)
}
save(p.value, file = paste0("/spin1/users/zhangh24/breast_cancer_data_analysis/discovery_SNP/conditional_analysis/result/conditional.check.result",i1,".Rdata"))
