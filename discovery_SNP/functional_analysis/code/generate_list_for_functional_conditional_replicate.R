setwd("/spin1/users/zhangh24/breast_cancer_data_analysis/")
load(paste0("./discovery_SNP/functional_analysis/result/functional_snp_conditional_replicate.Rdata"))
dim(functional_snp_conditional)
idx.na <- which(is.na(functional_snp_conditional$SNP.ICOG))

functional_snp_conditional <- functional_snp_conditional[-idx.na,]
total <- 0
n <- nrow(functional_snp_conditional)
p.value <- rep(0,n)
for(i1 in 1:100){
  print(i1)
  load(paste0("./discovery_SNP/functional_analysis/result/ICOG/p.value.result_replicate",i1,".Rdata"))
  temp <- length(p.value.result)
  p.value[total+(1:temp)] <- p.value.result
  total <- total + temp
}
# total
# n
functional_snp_conditional$p.value <- p.value

FunctionalFilter <- function(chr.temp,pos.temp,data){
  idx.temp <- which(data$CHR==chr.temp&
                      data$position==pos.temp)
  p.temp <- data$p.value[idx.temp]
  top.SNP <- data[idx.temp,12]
  idx <- which(data$CHR==chr.temp&
                 data$position>=(pos.temp-500000)&
                 data$position<=(pos.temp+500000)&
                 data$p.value<=(p.temp*100))
  lead.snp <- rep(top.SNP,length(idx))
  return(cbind(data[idx,],lead.snp))
}

############# 2 SNPs are used in the replication
############# first is rs6860806(pos 131640536 chr 5)
library(data.table)
idx.top <- which(functional_snp_conditional$CHR==5&
                   functional_snp_conditional$position==131640536)
idx.chr <- which(functional_snp_conditional$CHR==5)
functional_snp_conditional_chr <- functional_snp_conditional[idx.chr,]
p.top <- functional_snp_conditional$p.value[idx.top]
ccv1 <- FunctionalFilter(5,131640536,functional_snp_conditional_chr)



############# second is rs17743054 (pos 42900892 chr 18)
library(data.table)
idx.top <- which(functional_snp_conditional$CHR==18&
                   functional_snp_conditional$position==42900892)
idx.chr <- which(functional_snp_conditional$CHR==18)
functional_snp_conditional_chr <- functional_snp_conditional[idx.chr,]
p.top <- functional_snp_conditional$p.value[idx.top]
ccv2 <- FunctionalFilter(18,42900892,functional_snp_conditional_chr)
#try <- which(functional_snp_conditional_chr$p.value<=(1.82*10^-4))

ccv <- rbind(ccv1,ccv2)


write.csv(ccv,file= "./discovery_SNP/functional_analysis/result/functional_analysis_list_replicate_haoyu.csv",row.names=F,quote=F)

# SNP.list <- strsplit(functional.result$rs_id,split= ":")
# SNP <- rep("c",nrow(functional.result))
# for(i in 1:nrow(functional.result)){
#   SNP[i] <- SNP.list[[i]][1]
# }