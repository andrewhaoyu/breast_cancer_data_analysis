#plot permuate genotype results
library(qqman)
load("./known_SNPs/known_SNPs_analysis_G_revised/batch_effect_test/result/permutate_genoytp_result.rdata")

p <- result[,12,drop=F]
n <- length(p)
SNP <- paste0("rs",c(1:n))
CHR <- rep(1,n)
BP <- c(1,n)
gwas_result <- data.frame(SNP,CHR,BP,p)
colnames(gwas_result)[4] <- "P"
png(filename = "./known_SNPs/known_SNPs_analysis_G_revised/batch_effect_test/result/permuate_genotype.png",width = 8,height = 6, unit = "in",
    res = 300)
qq(gwas_result$P)
dev.off()


#plot all phenotype results
load("./known_SNPs/known_SNPs_analysis_G_revised/batch_effect_test/result/permutate_phenotype_result.rdata")

p <- c(result[,4],
       result[,6],
       result[,8],
       result[,10])
n <- length(p)
SNP <- paste0("rs",c(1:n))
CHR <- rep(1,n)
BP <- c(1,n)
gwas_result <- data.frame(SNP,CHR,BP,p)
png(filename = "./known_SNPs/known_SNPs_analysis_G_revised/batch_effect_test/result/permuate_phenotypes.png",width = 8,height = 6, unit = "in",
    res = 300)
qq(gwas_result$p)
dev.off()
