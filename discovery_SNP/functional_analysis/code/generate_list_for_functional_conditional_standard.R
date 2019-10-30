setwd("/data/zhangh24/breast_cancer_data_analysis/")
load(paste0("./discovery_SNP/functional_analysis/result/functional_snp_conditional_standard.Rdata"))
functional_snp_conditional <- functional_snp_conditional[1:5480,]
total <- 0
n <- nrow(functional_snp_conditional)
p.value <- rep(0,n)
for(i1 in 1:500){
  print(i1)
  load(paste0("./discovery_SNP/functional_analysis/result/ICOG/p.value.result_standard",i1,".Rdata"))
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

############# 2 SNPs detected by conditoinal standard model
############# first is rs11264454:156153043:A:G (pos 156153043 chr 1)
library(data.table)
idx.top <- which(functional_snp_conditional$CHR==1&
               functional_snp_conditional$position==156153043)
idx.chr <- which(functional_snp_conditional$CHR==1)
functional_snp_conditional_chr <- functional_snp_conditional[idx.chr,]
p.top <- functional_snp_conditional$p.value[idx.top]
ccv1 <- FunctionalFilter(1,156153043,functional_snp_conditional_chr)



############# second is rs141930488 (pos 51248274 chr 5)
library(data.table)
idx.top <- which(functional_snp_conditional$CHR==5&
                   functional_snp_conditional$position==51248274)
idx.chr <- which(functional_snp_conditional$CHR==5)
functional_snp_conditional_chr <- functional_snp_conditional[idx.chr,]
p.top <- functional_snp_conditional$p.value[idx.top]
ccv2 <- FunctionalFilter(5,51248274,functional_snp_conditional_chr)
#try <- which(functional_snp_conditional_chr$p.value<=(1.82*10^-4))

ccv <- rbind(ccv1,ccv2)


write.csv(ccv,file= "./discovery_SNP/functional_analysis/result/functional_analysis_list_standard_haoyu.csv",row.names=F,quote=F)

# SNP.list <- strsplit(functional.result$rs_id,split= ":")
# SNP <- rep("c",nrow(functional.result))
# for(i in 1:nrow(functional.result)){
#   SNP[i] <- SNP.list[[i]][1]
# }
