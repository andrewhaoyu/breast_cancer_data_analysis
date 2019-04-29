#simulate phenotypes data for AFR and EUR
#MAF based on 1000KG
#sample size EUR n =100000
#sample size AFR n = 15000
#heritability for EUR 0.8
#heritability for AFR 0.8
#5000 causal SNPs for each population
#4000 shared causal SNPs
#Genetic correlation for the shared SNPs is 0.6
#1000 SNPs for each as independent causal
# load("/spin1/users/zhangh24/KG.vcf/MAF_result/pruned_MAF.Rdata")
# set.seed(123)
# n.snp <- nrow(pruned.snp.clean)
# pruned.snp.permu <- pruned.snp.clean[sample(c(1:n.snp),n.snp),]
# save(pruned.snp.permu,file = "/spin1/users/zhangh24/KG.vcf/MAF_result/pruned_MAF_permu.Rdata")
set.seed(123)
load("/spin1/users/zhangh24/KG.vcf/MAF_result/pruned_MAF_permu.Rdata")
MAF.EUR <- pruned.snp.permu$MAF.EUR
MAF.AFR <- pruned.snp.permu$MAF.AFR
n.EUR <- 100000
n.AFR <- 15000
n.shared <- 4000
n.nonshare <- 1000
n.cau <- n.shared+n.nonshare
her.EUR <- 0.8
her.AFR <- 0.8
gr <- 0.6
library(mvtnorm)


Sigma <- matrix(c(her.EUR/n.cau,
                  gr*sqrt(her.EUR/n.cau*her.AFR/n.cau),
                  gr*sqrt(her.EUR/n.cau*her.AFR/n.cau),
                  her.AFR/n.cau),2,2)


beta <- rmvnorm(n.shared,c(0,0),
                sigma=Sigma)
#since the SNP MAF data is permutated
#for EUR and AFR, we will use 1:4000 as shared SNPs
#use 4000:5000 as nonshared SNPs for EUR
#5000:6000 as nonshared SNPs for AFR
G_EUR_shared <- matrix(rbinom(n.shared*n.EUR,2,MAF.EUR[1:n.shared]),n.EUR,n.shared,byrow = T)
G_AFR_shared <- matrix(rbinom)
#for easyness of coding, we will random simulate 2*nonshared SNPs
#but put second half for AFR as nonshare
two_nonshare = 2*n.nonshare
G_EUR_nonshare <- matrix(rbinom(two_nonshare*n.EUR,2,MAF.EUR[(n.shared+1):(n.shared+two_nonshare)]),n.EUR,two_nonshare,byrow = T)
G_EUR_effect <- cbind(G_EUR_shared,G_EUR_nonshare)


