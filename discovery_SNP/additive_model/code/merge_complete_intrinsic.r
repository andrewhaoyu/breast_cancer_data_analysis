setwd('/data/zhangh24/breast_cancer_data_analysis/')
result <- NULL
for(i1 in 1:35){
  load(paste0("./discovery_SNP/additive_model/result/intrinsic_subtype_logodds_complete_data_analysis_",i1,".Rdata"))

    result <- rbind(result,test.result.second.wald)

  
}
#SNP <- c(colnames(icog.julie),colnames(discovery.snp.icog)[1:18])

##################discovery snp were ordered based on the order they are extracted
discovery_snp <- read.csv("/data/zhangh24/breast_cancer_data_analysis/data/discovery_snp_summary_new.csv",header=T)
#SNP <- discovery_snp$SNP.ICOGS

################match the discovery snps to the order in the paper
discovery_snp_paper_order <- read.csv("./data/discovery_snp_paper_order.csv",header=T)
chr.pos.paper <- paste0(discovery_snp_paper_order$CHR,":",discovery_snp_paper_order$position)


chr.pos <- paste0(discovery_snp$CHR.x,":",discovery_snp$position)

idx.match <- match(chr.pos.paper,
                   chr.pos)
discovery_snp_new <- discovery_snp[idx.match,]
result <- result[idx.match,]
final_result <- cbind(discovery_snp_paper_order[,1],result)

write.csv(final_result,file= "./discovery_SNP/additive_model/result/complete_data_intrinsic.csv")
