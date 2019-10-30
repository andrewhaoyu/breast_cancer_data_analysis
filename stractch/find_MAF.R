setwd("/data/zhangh24/breast_cancer_data_analysis/")
load("/data/zhangh24/breast_cancer_data_analysis/whole_genome/ICOG/ERPRHER2GRADE_casecase/result/meta_result_shared_1p.Rdata")

discovery_icog_id <- read.table(paste0("./discovery_SNP/result/extract_id_icog_discovery.txt"),header=T)

idx <- which(meta_result_shared_1p$SNP.ICOGS==as.character(discovery_icog_id[21,1]))
meta_result_shared_1p[idx,]

setwd('/data/zhangh24/breast_cancer_data_analysis/')
load(paste0("./discovery_SNP/CIMBA_BCAC_meta_analysis/result/CIMBA.BCAC.combine.Rdata"))
idx <- which(CIMBA.BCAC.combine$SNP.ICOGS==as.character(discovery_icog_id[11,1]))
CIMBA.BCAC.combine[idx,]




metaF <- function(x1,s1,x2,s2){
  s.result <- sqrt(((s1^2)^-1+(s2^2)^-1)^-1)
  x.result <- ((s1^2)^-1+(s2^2)^-1)^-1*((s1^2)^-1*x1+(s2^2)^-1*x2)
  return(list(x.result,
              s.result))
}




z <- -0.3790843/ 0.06676339


metaF(-0.400307,sqrt(0.007654879),-0.3495,0.1033)



2*pnorm(abs(z),lower.tail = F)


