setwd("/spin1/users/zhangh24/breast_cancer_data_analysis/")
load(paste0("./discovery_SNP/functional_analysis/result/functional_snp_conditional.Rdata"))
total <- 0
n <- nrow(functional_snp_conditional)
p.value <- rep(0,n)
for(i1 in 1:1000){
  print(i1)
  load(paste0("./discovery_SNP/functional_analysis/result/ICOG/p.value.result",i1,".Rdata"))
  temp <- length(p.value.result)
  p.value[total+(1:temp)] <- p.value.result
  total <- total + temp
}
# total
# n
functional_snp_conditional$p.value <- p.value

load("/spin1/users/zhangh24/breast_cancer_data_analysis/whole_genome_age/ICOG/ERPRHER2GRADE_casecase/result/meta_result_shared_1p.Rdata")
meta_result_shared_1p_MTOP <- meta_result_shared_1p
idx <- which(meta_result_shared_1p$CHR==22&
               meta_result_shared_1p$position==30592808)
meta_result_shared_1p[idx,]
############# 5 SNPs only detected by mixed effect two-stage model
############# 2 SNPs only detected by fixed effect two-stage model
############# 4 SNPs only detected by conditoinal mixed effet two-stage model
library(data.table)
MTOP_SNP <- fread("./data/Discovery_SNP_MTOP.csv",header=T)
functional.result <- NULL

FunctionalFilter <- function(chr.temp,pos.temp,data){
  idx.temp <- which(data$CHR==chr.temp&
                      data$position==pos.temp)
  p.temp <- data$p.value[idx.temp]
  
  idx <- which(data$CHR==chr.temp&
                 data$position>=(pos.temp-50000)&
                 data$position<=(pos.temp+50000)&
                 data$p.value<=(p.temp*100))
  return(data[idx,])
}


for(i in 1:nrow(MTOP_SNP)){
  print(i)
  chr.temp = as.numeric(MTOP_SNP[i,3])
  pos.temp = as.numeric(MTOP_SNP[i,2])
  result.temp <- FunctionalFilter(chr.temp,pos.temp,meta_result_shared_1p)
  print(nrow(result.temp))
  functional.result <- rbind(functional.result,result.temp)
  
}

############# 2 SNPs only detected by fixed effect two-stage model
library(data.table)
FTOP_SNP <- fread("./data/Discovery_SNP_FTOP.csv",header=T)
#functional.result <- NULL
load("/spin1/users/zhangh24/breast_cancer_data_analysis/whole_genome_age/ICOG/ERPRHER2GRADE_fixed_baseline/result/meta_result_shared_1p.Rdata")
for(i in 1:nrow(FTOP_SNP)){
  print(i)
  chr.temp = as.numeric(FTOP_SNP[i,3])
  pos.temp = as.numeric(FTOP_SNP[i,2])
  result.temp <- FunctionalFilter(chr.temp,pos.temp,meta_result_shared_1p)
  print(nrow(result.temp))
  functional.result <- rbind(functional.result,result.temp)
  
}

############# 4 SNPs only detected by conditoinal mixed effet two-stage model

MTOP_con_SNP <- fread("./data/Discovery_SNP_conditional_MTOP.csv",header=T)
#functional.result <- NULL

for(i in 1:nrow(MTOP_con_SNP)){
  print(i)
  chr.temp = as.numeric(MTOP_con_SNP[i,3])
  pos.temp = as.numeric(MTOP_con_SNP[i,2])
  result.temp <- FunctionalFilter(chr.temp,pos.temp,functional_snp_conditional)
  print(nrow(result.temp))
  functional.result <- rbind(functional.result,result.temp)
  
}
