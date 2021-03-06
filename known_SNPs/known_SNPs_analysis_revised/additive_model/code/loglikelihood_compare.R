data <- read.csv("./known_SNPs/known_SNPs_analysis_revised/additive_model/result/loglikelihood_difference.csv",header=T)

log.difference <- data[,1]-data[,2]
sum(log.difference>=0)
sum(log.difference<=0)
hist(log.difference,breaks = 20,main = "log-likelihood difference between additive model and intrinsic subtypes")
t.test(log.difference)
